#ifndef _OCS_PTCS_TYPES_INCLUDE_
#define _OCS_PTCS_TYPES_INCLUDE_ 1

#include "ptcs.h"
#define PTCS_PATTERNLEN 50

typedef enum motiontype { no_motion, once, group, discrete, continuous } motiontype;

typedef struct OffsetPointType {
    char system[PTCS_SYSLEN];
    char type[PTCS_TYPELEN];    
    double dc1;
    double dc2;
} OffsetPointType;

typedef struct ScanParamType {
    double height;
    double width;
    double optimalPa;
    double velocity;
    double dy;
    char system[PTCS_SYSLEN];
    char type[PTCS_TYPELEN];
    long reversal;
    OffsetPointType offset;
    char pattern[PTCS_PATTERNLEN];
    long nterms;
} ScanParamType; 
typedef ScanParamType * ScanParamPtr;

/* sequence step structure */

typedef struct StepType {
  unsigned long number;      /* RTS step number */
  double tai;                /* TAI time at the midpoint of the integration step */
  char source[PTCS_NAMELEN]; /* Source name  */
  long index;                /* Index into observing area component  */
  char beam[PTCS_NAMELEN];   /* Beam name  */
  double stepTime;           /* RTS step time for this sequence  */
  int stspl_index;           /* Index location at which this step should be placed in the STATE array */
  int publish;               /* publish = 1 => After placing this step into STATE, publish the STATE parameter */
  struct StepType *next;     /* Pointer to the next step in the step queue */
} StepType;

typedef StepType * StepPtr;

/* spherical coordinate structure */

typedef struct SphCoordType {
    double c1;
    double c2;
} SphCoordType;

/* position structure */

typedef struct PosnType {
    char system[PTCS_SYSLEN];
    double fpAngle;
    double rotAngle;
    SphCoordType actual;
    SphCoordType demand;
    SphCoordType base;
} PosnType;

/* time-stamped position structure */

typedef struct TSPosnType {
    double tai;
    double airmass;
    char instAp[PTCS_NAMELEN];
    PosnType tracking;
    PosnType azel;
} TSPosnType;

typedef TSPosnType * TSPosnPtr;

typedef enum pong_segment_type { straight_segment, arc_segment } pong_segment_type;
typedef enum pong_type { straight_pong, curvy_pong, rounded_pong } pong_type;

typedef struct PongStraightSegment {
  double xstart;
  double ystart;
  double xend;
  double yend;
  double mx;
  double my;
} PongStraightSegmentType;
typedef PongStraightSegmentType * PongStraightSegmentPtr;

typedef struct PongArcSegment {
  double xcenter;
  double ycenter;
  double radius;
  double theta_init;
  double theta_fin;
  double rot_dir;
} PongArcSegmentType;
typedef PongArcSegmentType * PongArcSegmentPtr;

typedef struct PongSegment {
  pong_segment_type type;
  PongStraightSegmentType straight;
  PongArcSegmentType arc;
  double tstart;
  double tend;
} PongSegmentType;
typedef PongSegmentType * PongSegmentPtr;

typedef struct PongScan {
  char type[PTCS_NAMELEN];
  int n_seg;
  double grid_space;
  int x_ngridseg;
  int y_ngridseg;
  double *xvertex;
  double *yvertex;
  PongSegmentPtr path_data;
  int numPathSegments;
  double alpha_x;
  double alpha_y;
  double beta_x;
  double beta_y;
  double period;
} PongScanType;
typedef PongScanType * PongScanPtr;

typedef struct ScanRow {
  long index;
  double sx;
  double sy;
  double ex;
  double ey;
  double mx;
  double my;
  long min;
} ScanRowType;
typedef ScanRowType * ScanRowPtr;

typedef enum scanPatternType { no_scan, raster, discrete_boustrophedon, continuous_boustrophedon, pong, skydip } scanPatternType;

typedef enum skydipType { no_skydip, discrete_skydip, continuous_skydip } skydipType;

typedef enum skydipDirectionType { forward, reverse } skydipDirectionType;

#endif
