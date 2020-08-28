#ifndef SPPROJECT_DIVISIONUTILLIS_H
#define SPPROJECT_DIVISIONUTILLIS_H
#include "list.h"
#include "spmat.h"
#define EPSILON 0.00001
#define IS_POSITIVE(X) ((X) > 0.00001)

spmat *BuildSubMatrix(spmat *Amatrix, ELEMENT *g_ListNodes);
void CreateBgRowSumsArray(spmat *A_SubMatrix,int *rowptr, int M, int sumOfSubDegArray, int subGroupSize, double *bg_rowSumsArray);
int CreateSubDegArray(spmat * A_SubMatrix, ELEMENT *g_ListNodes, int *subDegArrayPtr, int *deg_array);
int CalculateNNZofSubMatrix(ELEMENT *g_CurrentNode, int *rowptr);
void CreateZeroMatrix(spmat *A_SubMatrix);
double CalculateSumOfKiVi(double *vector, int *deg_array, int M, int subGroupSize);
void PowerIteration(spmat *Bhat, double* initialVec, double *eigenVec, double *eigenVal, double norma);
void RandomizedVec(double *initialVec,int size);
int CheckEpsilon(double *vecbOld, double *vecbNew, int size);
double CalculateEigenValue(spmat *Bhat, double *eigenVec, double norma, double *tempVec);
double VecWithVecMult(const double *vec1, const double *vec2, int size);
double CalculateNorma(spmat *Bhat);
void CalculateSvector(double *S_vec, double *vector, int size, int isOnesVector);

#endif /*SPPROJECT_DIVISIONUTILLIS_H */
