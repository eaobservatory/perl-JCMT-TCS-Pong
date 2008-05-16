/*                      o c s P t c s C o m p u t e P o n g G r i d

 *  Function name:
      ocsPtcsComputePongGrid

 *  Function:
        Compute the vertices (i.e., reflection points) of the PONG pattern

 *  Description:
        Compute the vertices (i.e., reflection points) of the PONG pattern

 *  Language:
      C

 *  Call:
      (void) =  ocsPtcsComputePongGrid(&scanParam, &pongScan, &status)

 *  Parameters:   (">" input, "!" modified, "W" workspace, "<" output)
      (>) scanParam (ScanParamPtr)   Scan parameters from configuration file
      (!) pongScan  (PongScanPtr)    Computed pong scan parameters
      (!) status    (StatusType *)   Modified status.

 *  External functions used:
      none

 *  External values used:
      none

 *  Prior requirements:
        As above.

 *  Support: Russell Kackley, JAC

 *

 *  History:
      17-May-2007 - RDK - Original version

 */

/* Added during import to XS - so that we do not need to be dependent
   on the whole of TCS and DRAMA - TJ */
#include <stdlib.h>
typedef int StatusType;
typedef int SdsIdType;
typedef int systype;
typedef int VtType;
/* we do not really care about error status */
#define PTCS__MALLOCERR -1
#define PTCS__GERROR -2
/* 2pi */
#define D2PI 6.2831853071795864769252867665590057683943387987502
#define StatusOkP(_value_)  (*(_value_) == STATUS__OK)

#include <math.h>
/* Commented out during import 
#include "DitsFix.h"
#include "slamac.h"
#include "jit.h"
#include "ptcs_err.h"
*/
#include "ocsPtcsScan.h"
#include "ocsPtcsPong.h"



void ocsPtcsComputePongGrid( ScanParamPtr scanParam, PongScanPtr pongScan, StatusType *status ) {

  double height = scanParam->height, width = scanParam->width;
  double spacing = scanParam->dy, velocity = scanParam->velocity;
  double vert_spacing;
  int x_numvert, y_numvert;
  int *most, *least;
  int i;
  double grid_space;
  int x_ngridseg, y_ngridseg, n_seg;
  double *xgrid, *ygrid;
  double *xvertex, *yvertex;
  int x_init, y_init;
  int x_off, y_off;
  int steps;
  double v2;

  if (!StatusOkP(status)) return;

  /* Determine number of vertices (reflection points) along each side
     of the box which satisfies the common-factors criterion and the
     requested size / spacing. The algorithm to determine the vertices
     is based on pseudocode provided in SC2/ANA/S210/008 "The Impact
     of Scanning Pattern Strategies on Uniform Sky Coverage of Large
     Maps". */

  /* Find out how far apart the vertices are along the axes */
  vert_spacing = ( 2.0 * spacing ) / sqrt ( 2.0 );

  /* Determine how many vertices (minimum) there must be along each
     axis to cover the required area */
  x_numvert = ceil ( width / vert_spacing );
  y_numvert = ceil ( height / vert_spacing );

  /* Determine which is lower and check to make sure that one is even
     while the other is odd */
  if ( x_numvert >= y_numvert ) {
    most = &x_numvert;
    least = &y_numvert;
  } else {
    most = &y_numvert;
    least = &x_numvert;
  }

  /* If both are odd or both are even, increment the lesser of the
     two, and update which is least */
  if ( ( x_numvert % 2 ) == (y_numvert % 2 ) ) {
    *least += 1;
    if ( x_numvert >= y_numvert ) {
      most = &x_numvert;
      least = &y_numvert;
    } else {
      most = &y_numvert;
      least = &x_numvert;
    }
  }

  /* Check for common factors between the two, and adjust as necessary
     until x_numvert and y_numvert do not share any common factors */
  for (i = 3; i <= *least; i++) {
    if ( ( ( *least % i ) == 0 ) && ( ( *most % i ) == 0 ) ) {
      /* Found a common factor, increment most and start over */
      *most += 2;
      i = 3;
    }
  }

  /* At the end of this step, we are left with two values, x_numvert
     and y_numvert, of which one is even and one is odd, and which
     share no common factors. */

  /* The entire pattern is defined on a grid with points spaced half
     the distance between adjacent vertices as calculated above along
     the same side. Calculate spacing between these grid points and
     the length along each side in units of grid_space sized
     segments. */
  grid_space = vert_spacing / 2.;
  x_ngridseg = x_numvert * 2;
  y_ngridseg = y_numvert * 2;

  /* Calculate the total number of straight line segments in the
     pattern, and allocate an array to store the x and y tanplane
     offsets of the endpoints in order */
  n_seg = x_ngridseg + y_ngridseg;
/*   if (jitDebugGet() & PTCS_DEBUG_LEVEL3) */
/*     MsgOut(status, "x_numvert %d y_numvert %d n_seg %d", x_numvert, y_numvert, n_seg); */

  xgrid = calloc(x_ngridseg + 1, sizeof(double));
  if (xgrid == NULL) {
    *status = PTCS__MALLOCERR;
/*     ErsRep(0, status, "Failed to calloc space for xgrid"); */
    return;
  }
  ygrid = calloc(y_ngridseg + 1, sizeof(double));
  if (ygrid == NULL) {
    *status = PTCS__MALLOCERR;
/*     ErsRep(0, status, "Failed to calloc space for ygrid"); */
    return;
  }

  for (i = 0; i <= x_ngridseg; i++) {
    xgrid[i] = (i - x_ngridseg/2) * grid_space;
  }
  for (i = 0; i <= y_ngridseg; i++) {
    ygrid[i] = (i - y_ngridseg/2) * grid_space;
  }
 
  xvertex = calloc(n_seg + 1, sizeof(double));
  if (xvertex == NULL) {
    *status = PTCS__MALLOCERR;
/*     ErsRep(0, status, "Failed to calloc space for xvertex"); */
    return;
  }
  yvertex = calloc(n_seg + 1, sizeof(double));
  if (yvertex == NULL) {
    *status = PTCS__MALLOCERR;
/*     ErsRep(0, status, "Failed to calloc space for yvertex"); */
    return;
  }

  /* Initialization of the scan pattern */
  x_init = 0;                /* starting grid coordinates */
  y_init = 1;

  x_off = x_init;            /* current grid offsets */
  y_off = y_init;

  xvertex[0] = xgrid[x_init]; /* starting tanplane offsets of the scan */
  yvertex[0] = ygrid[y_init];

  /* Loop over line segments and calculate list of endpoints for each
     reflection in order */
  steps = 0;
  for (i = 1; i <= n_seg; i++) {
    int x_refl, y_refl;

    /* increment steps to the next boundary reflection (whichever
       comes first - along the sides or top/bottom of the
       rectangle) */

    if( (x_ngridseg - x_off) <= (y_ngridseg - y_off) ) {
      steps += x_ngridseg-x_off;    /* side reflection next */
    } else {
      steps += y_ngridseg-y_off;    /* top/bottom reflection next */
    }

    /* Number of steps divided by number of grid segments along each
       side tells us how many reflections we've been through. */
    x_refl = floor( (steps+x_init)/x_ngridseg );
    y_refl = floor( (steps+y_init)/y_ngridseg );

    /* The remainder is the grid offset */
    x_off = (steps+x_init) % x_ngridseg;
    y_off = (steps+y_init) % y_ngridseg;

    /* If the number of reflections is even, the offset is in the
       positive direction from the start of the lookup table. If it is
       odd, the offset is in the negative direction from the end of
       the lookup table */
    if( x_refl % 2 == 0 ) {
      xvertex[i] = xgrid[x_off];               /* even reflections */
    } else {
      xvertex[i] = xgrid[x_ngridseg-x_off];    /* odd reflections */
    }

    if( y_refl % 2 == 0 ) {
      yvertex[i] = ygrid[y_off];               /* even reflections */
    } else {
      yvertex[i] = ygrid[y_ngridseg-y_off];    /* odd reflections */
    }
  }

  if (xgrid != NULL) free(xgrid);
  if (ygrid != NULL) free(ygrid);

  pongScan->n_seg = n_seg;
  pongScan->grid_space = grid_space;
  pongScan->x_ngridseg = x_ngridseg;
  pongScan->y_ngridseg = y_ngridseg;
  pongScan->xvertex = xvertex;
  pongScan->yvertex = yvertex;

  /* Compute the coefficients for the CURVY_PONG path */
  pongScan->alpha_x = x_numvert * vert_spacing / 2.0;
  pongScan->alpha_y = y_numvert * vert_spacing / 2.0;
  v2 = velocity * sqrt(2.0) / 2.0;
  pongScan->beta_x = x_numvert * vert_spacing * 2.0 / v2;
  pongScan->beta_y = y_numvert * vert_spacing * 2.0 / v2;

  /* Compute the time required to complete one complete circuit of the
     CURVY_PONG path */
  pongScan->period = x_numvert * y_numvert * vert_spacing * 2.0 / v2;
}

/*                      o c s P t c s C o m p u t e P o n g P a t h

 *  Function name:
      ocsPtcsComputePongPath

 *  Function:
        Compute the details of the SQUARE_PONG or ROUNDED_PONG path

 *  Description:
        Compute the details of the SQUARE_PONG or ROUNDED_PONG path

 *  Language:
      C

 *  Call:
      (void) =  ocsPtcsComputePongPath(&scanParam, &pongScan, &status)

 *  Parameters:   (">" input, "!" modified, "W" workspace, "<" output)
      (>) scanParam (ScanParamPtr)   Scan parameters from configuration file
      (!) pongScan  (PongScanPtr)    Computed pong scan parameters
      (!) status    (StatusType *)   Modified status.

 *  External functions used:
      none

 *  External values used:
      none

 *  Prior requirements:
        As above.

 *  Support: Russell Kackley, JAC

 *

 *  History:
      17-May-2007 - RDK - Original version

 */

#define ARC_LENGTH_NOM 45 /* 45 deg arc at start and end of straight-line segment for ROUNDED_PONG */;
#define EPSILON 1.e-6     /* EPSILON is used to determine if a path is basically zero-length */

void ocsPtcsComputePongPath( ScanParamPtr scanParam, PongScanPtr pongScan, StatusType *status ) {

  double spacing = scanParam->dy, velocity = scanParam->velocity;
  double curve_radius = spacing / 2;
  int pathSegmentIndex = -1;
  double pathLength = 0.0;
  PongSegmentPtr path_data;
  double *xvertex = pongScan->xvertex;
  double *yvertex = pongScan->yvertex;
  double grid_space = pongScan->grid_space;
  int x_ngridseg = pongScan->x_ngridseg;
  int y_ngridseg = pongScan->y_ngridseg;
  double xstart_nom, ystart_nom;
  double xend_nom, yend_nom;
  double dx_nom, dy_nom, length_nom;
  double mx, my;
  double xstart, ystart;
  double xend, yend;
  double dx, dy, length;
  double arc_length = ARC_LENGTH_NOM;
  double sqrt2 = sqrt(2.);
  double posx_boundary = grid_space * x_ngridseg / 2;
  double negx_boundary = -posx_boundary;
  double posy_boundary = grid_space * y_ngridseg / 2;
  double negy_boundary = -posy_boundary;
  double segLength;
  int i;

  if (!StatusOkP(status)) return;
/*   if (jitDebugGet() & PTCS_DEBUG_LEVEL3) MsgOut(status, "grid_space %f x_ngridseg %d y_ngridseg %d",grid_space,x_ngridseg,y_ngridseg); */
  path_data = calloc(3 * pongScan->n_seg, sizeof(PongSegmentType));
  if (path_data == NULL) {
    *status = PTCS__MALLOCERR;
/*     ErsRep(0, status, "Failed to calloc space for path_data"); */
    return;
  }

  /* TSTART for the first path segment is 0.0 */
  path_data[0].tstart = 0.0;
  
  for (i = 0; i < pongScan->n_seg; i++) {
    double xstart_arc_center = 0., ystart_arc_center = 0.;
    double start_theta_init = 0., start_theta_fin=0., start_rot_dir = 0.;
    double xend_arc_center = 0., yend_arc_center = 0.;
    double end_theta_init = 0., end_theta_fin=0., end_rot_dir = 0.;

    /* Get the nominal start and end points for this row */
    xstart_nom = xvertex[i];
    ystart_nom = yvertex[i];

    xend_nom = xvertex[i+1];
    yend_nom = yvertex[i+1];
/*     if (jitDebugGet() & PTCS_DEBUG_LEVEL3) MsgOut(status, "xstart_nom %f ystart_nom %f xend_nom %f yend_nom %f ",xstart_nom,ystart_nom,xend_nom,yend_nom); */

    /* Compute the length and direction vector (mx, my) of the row */
    dx_nom = xend_nom - xstart_nom;
    dy_nom = yend_nom - ystart_nom;
    length_nom = sqrt(dx_nom * dx_nom + dy_nom * dy_nom);
    mx = dx_nom / length_nom;
    my = dy_nom / length_nom;

    /* For SQUARE_PONG, we use the nominal start and end points for
       the actual path segment. For ROUNDED_PONG, we modify the path
       segment by shortening the straight part of the path and putting
       arcs at the beginning and end of the path segment. */
    if (strcmp(pongScan->type, "SQUARE") == 0) {
      xstart = xstart_nom;
      ystart = ystart_nom;
      xend = xend_nom;
      yend = yend_nom;
      dx = xend - xstart;
      dy = yend - ystart;
      length = sqrt(dx * dx + dy * dy);
    } else {

    /* For ROUNDED_PONG, the actual start and end points are computed
       from the nominal start and end points, modified by the arcs
       placed at the start and end to round off the sharp corners. */
      xstart = xstart_nom + mx * curve_radius;
      ystart = ystart_nom + my * curve_radius;
      xend = xend_nom - mx * curve_radius;
      yend = yend_nom - my * curve_radius;
      dx = xend - xstart;
      dy = yend - ystart;
      length = sqrt(dx * dx + dy * dy);
      
      /* If the straight-line segment has a non-zero length, where
	 zero is defined as less than EPSILON, use the nominal arc
	 length. If the straight-line segment is basically zero, then
	 we don't even use a straight-line segment and combine the
	 start and end arcs into a single arc of twice the nominal arc
	 length */
      if (length > EPSILON) {
	arc_length = ARC_LENGTH_NOM;
      } else {
	arc_length = 2.0 * ARC_LENGTH_NOM;
      }

      /* For ROUNDED_PONG, the arcs at the start and end of the
	 straight segment are the most complicated parts of the
	 path. We need to figure out the angular position at which the
	 arc starts and the direction of motion while in the
	 arc. rot_dir == 1 => CCW; rot_dir == -1 => CW. These values
	 depend on whether the start/end points are on the +X, -X, +Y,
	 or -Y boundaries, and the direction of the path
	 leaving/approaching the boundary. */
      if (fabs(xstart_nom - posx_boundary) < EPSILON) {
	xstart_arc_center = xstart_nom - curve_radius * sqrt2;
	ystart_arc_center = ystart_nom;
	start_theta_init = 0.;
	start_rot_dir = (dy_nom > 0)? 1.: -1.;
      } else if (fabs(xstart_nom - negx_boundary) < EPSILON) {
	xstart_arc_center = xstart_nom + curve_radius * sqrt2;
	ystart_arc_center = ystart_nom;
	start_theta_init = 180.;
	start_rot_dir = (dy_nom > 0)? -1.: 1.;
      } else if (fabs(ystart_nom - posy_boundary) < EPSILON) {
	xstart_arc_center = xstart_nom;
	ystart_arc_center = ystart_nom - curve_radius * sqrt2;
	start_theta_init = 90.;
	start_rot_dir = (dx_nom > 0)? -1.: 1.;
      } else if (fabs(ystart_nom - negy_boundary) < EPSILON) {
	xstart_arc_center = xstart_nom;
	ystart_arc_center = ystart_nom + curve_radius * sqrt2;
	start_theta_init = 270.;
	start_rot_dir = (dx_nom > 0)? 1.: -1.;
      } else {
	*status = PTCS__GERROR;
/* 	ErsRep(0, status, "Unexpected error in ocsPtcsComputeRoundedPongPath: xstart_nom or ystart_nom do not match expected values xstart_nom %f expected %f diff1 %g diff2 %g ystart_nom %f expected %f ",xstart_nom,grid_space * x_ngridseg / 2,xstart_nom-grid_space * x_ngridseg / 2,xstart_nom+grid_space * x_ngridseg / 2,ystart_nom,grid_space * y_ngridseg / 2); */
      }
      start_theta_fin = start_theta_init + start_rot_dir * arc_length;     
      
      if (fabs(xend_nom - posx_boundary) < EPSILON) {
	xend_arc_center = xend_nom - curve_radius * sqrt2;
	yend_arc_center = yend_nom;
	end_rot_dir = (dy_nom > 0)? 1.: -1.;
	end_theta_init = 0. - end_rot_dir * arc_length;
      } else if (fabs(xend_nom - negx_boundary) < EPSILON) {
	xend_arc_center = xend_nom + curve_radius * sqrt2;
	yend_arc_center = yend_nom;
	end_rot_dir = (dy_nom > 0)? -1.: 1.;
	end_theta_init = 180. - end_rot_dir * arc_length;
      } else if (fabs(yend_nom - posy_boundary) < EPSILON) {
	xend_arc_center = xend_nom;
	yend_arc_center = yend_nom - curve_radius * sqrt2;
	end_rot_dir = (dx_nom > 0)? -1.: 1.;
	end_theta_init = 90. - end_rot_dir * arc_length;
      } else if (fabs(yend_nom - negy_boundary) < EPSILON) {
	xend_arc_center = xend_nom;
	yend_arc_center = yend_nom + curve_radius * sqrt2;
	end_rot_dir = (dx_nom > 0)? 1.: -1.;
	end_theta_init = 270. - end_rot_dir * arc_length;
      } else {
	*status = PTCS__GERROR;
/* 	ErsRep(0, status, "Unexpected error in ocsPtcsComputeRoundedPongPath: xend_nom or yend_nom do not match expected values xend_nom %f expected %f yend_nom %f expected %f ",xend_nom,grid_space * x_ngridseg / 2,yend_nom,grid_space * y_ngridseg / 2); */
      }
      end_theta_fin = end_theta_init + end_rot_dir * arc_length;
    }

    /* For ROUNDED_PONG, put an arc at the start of the segment */
    if (strcmp(pongScan->type, "ROUNDED") == 0) {
      pathSegmentIndex++;
      path_data[pathSegmentIndex].type = arc_segment;
      path_data[pathSegmentIndex].arc.xcenter = xstart_arc_center;
      path_data[pathSegmentIndex].arc.ycenter = ystart_arc_center;
      path_data[pathSegmentIndex].arc.radius = curve_radius;
      path_data[pathSegmentIndex].arc.theta_init = start_theta_init;
      path_data[pathSegmentIndex].arc.theta_fin = start_theta_fin;
      path_data[pathSegmentIndex].arc.rot_dir = start_rot_dir;
      segLength = D2PI * curve_radius * arc_length / 360.;
      pathLength += segLength;
      /* TSTART for path segments after the first is the TEND from the
	 previous segment */
      if (pathSegmentIndex > 0) {
	path_data[pathSegmentIndex].tstart = path_data[pathSegmentIndex - 1].tend;
      }
      path_data[pathSegmentIndex].tend = path_data[pathSegmentIndex].tstart + segLength / velocity;
    }

    /* The straight part of the segment. Skip any segments that are
       basically zero-length. */
    if (length > EPSILON) {
      pathSegmentIndex++;
      path_data[pathSegmentIndex].type = straight_segment;
      path_data[pathSegmentIndex].straight.xstart = xstart;
      path_data[pathSegmentIndex].straight.ystart = ystart;
      path_data[pathSegmentIndex].straight.xend = xend;
      path_data[pathSegmentIndex].straight.yend = yend;
      path_data[pathSegmentIndex].straight.mx = mx;
      path_data[pathSegmentIndex].straight.my = my;
      pathLength += length;
      /* TSTART for path segments after the first is the TEND from the
	 previous segment */
      if (pathSegmentIndex > 0) {
	path_data[pathSegmentIndex].tstart = path_data[pathSegmentIndex - 1].tend;
      }
      path_data[pathSegmentIndex].tend = path_data[pathSegmentIndex].tstart + length / velocity;
    }

    /*  For ROUNDED_PONG where the straight-line segment is non-zero,
	put an arc at the end of the segment */
    if (strcmp(pongScan->type, "ROUNDED") == 0 && length > EPSILON) {
      pathSegmentIndex++;
      path_data[pathSegmentIndex].type = arc_segment;
      path_data[pathSegmentIndex].arc.xcenter = xend_arc_center;
      path_data[pathSegmentIndex].arc.ycenter = yend_arc_center;
      path_data[pathSegmentIndex].arc.radius = curve_radius;
      path_data[pathSegmentIndex].arc.theta_init = end_theta_init;
      path_data[pathSegmentIndex].arc.theta_fin = end_theta_fin;
      path_data[pathSegmentIndex].arc.rot_dir = end_rot_dir;
      segLength = D2PI * curve_radius * arc_length / 360.;
      pathLength += segLength;
      path_data[pathSegmentIndex].tstart = path_data[pathSegmentIndex - 1].tend;
      path_data[pathSegmentIndex].tend = path_data[pathSegmentIndex].tstart + segLength / velocity;
    }
  }

  /* Save the number of path segments and the pointer to the path_data
     array */
  pongScan->numPathSegments = pathSegmentIndex + 1;
  pongScan->path_data = path_data;

  /* Compute the time required to complete one complete circuit of the
     SQUARE or ROUNDED PONG path */
  pongScan->period = pathLength / velocity;
}

