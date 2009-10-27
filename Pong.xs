/*        -*- C -*-

   Pong patterns

						timj@jach.hawaii.edu

   Copyright (C) 2008 Science and Technology Facilities Council.
   All Rights Reserved.

*/

#include "EXTERN.h"   /* std perl include */
#include "perl.h"     /* std perl include */
#include "XSUB.h"     /* XSUB include */

#include "ocsPtcsTypes.h"

MODULE = JCMT::TCS::Pong     PACKAGE = JCMT::TCS::Pong

void
_get_pong_dur( height, width, dy, velocity, type)
  double height
  double width
  double dy
  double velocity
  char* type
 PREINIT:
  ScanParamType scanParam;
  PongScanType pongScan;
  int status = 0;
 PPCODE:
  scanParam.height = height;
  scanParam.width = width;
  scanParam.dy = dy;
  scanParam.velocity = velocity;
  strcpy(scanParam.pattern, "PONG");
  strcpy(scanParam.type, type );

  /* Compute the vertices (i.e., reflection points) for the PONG
     pattern. This also computes the period of the CURVY_PONG
     pattern. */
  ocsPtcsComputePongGrid( &scanParam, &pongScan, &status);
  
  /* For SQUARE and ROUNDED PONG, compute the details of the path. The
     period is re-computed for the SQUARE and ROUNDED_PONG patterns
     because they are slightly different than the CURVY_PONG
     pattern. The period is returned in pongScan.period. */
  if (strcmp(pongScan.type, "SQUARE") == 0 || strcmp(pongScan.type, "ROUNDED") == 0) {
    ocsPtcsComputePongPath( &scanParam, &pongScan, &status);
  }

  /* We no longer need the vertex information, so free the memory
     now */
  if (pongScan.xvertex != NULL) {
    free(pongScan.xvertex);
    pongScan.xvertex = NULL;
  }
  if (pongScan.yvertex != NULL) {
    free(pongScan.yvertex);
    pongScan.yvertex = NULL;
  }

  /* The time required to traverse the entire pattern is in
     pongScan.period */
  XPUSHs( sv_2mortal(newSVnv(pongScan.period)));
