#ifndef _OCS_PTCS_PONG_INCLUDE_
#define _OCS_PTCS_PONG_INCLUDE_ 1

#include "ocsPtcsTypes.h"

void ocsPtcsPongPath( VtType vt, double tai, double taioffz, double *udiffa, double *udiffb, StatusType *status );
void ocsPtcsComputePongGrid( ScanParamPtr scanParam, PongScanPtr pongScan, StatusType *status );
void ocsPtcsComputePongPath( ScanParamPtr scanParam, PongScanPtr pongScan, StatusType *status );

#endif
