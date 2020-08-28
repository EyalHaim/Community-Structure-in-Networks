#ifndef SPPROJECT_MODULARITYMAXIMIZATION_H
#define SPPROJECT_MODULARITYMAXIMIZATION_H
#include "spmat.h"


void ModularityMaximization(spmat *Bg, double *originalS);
void InitVectorsForModularityMaximization(int **unMoved, double **sumsOfCols,int **indices, int subSizeGroup);
void FreeModularityMaximization(int *unMoved, double *sumsOfCols,int *indices);




#endif /*SPPROJECT_MODULARITYMAXIMIZATION_H*/
