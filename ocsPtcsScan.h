#ifndef _OCS_PTCS_SCAN_INCLUDE_
#define _OCS_PTCS_SCAN_INCLUDE_ 1

#include "ocsPtcsTypes.h"

void ocsPtcsConfigPrepareScan(StatusType *status);
void ocsPtcsScanGetBeta(double *beta, StatusType *status);
void ocsPtcsComputeFirstScanRow(long index, SdsIdType *retArg, StatusType *status);
void ocsPtcsComputeScanRow(long index, SdsIdType *retArg, StatusType *status);
void ocsPtcsComputeScanRowFinish(double sx, double sy, double ex, double ey, double mx, double my, SdsIdType *retArg, StatusType *status);
void ocsPtcsComputeRemainingRows(long *num, StatusType *status);
void ocsPtcsComputeScanGetSystem(char *name, systype *sys, char *category, StatusType *status);

#endif
