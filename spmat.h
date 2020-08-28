#ifndef _SPMAT_H
#define _SPMAT_H
#include "list.h"

typedef struct _spmat {
    int		 n;
    int      M;
    int 	 *colind;
    int      *rowptr;
    int      lastInsertedIndex;
    int      *deg_array;    /* Array of Ki values */
    double   *Bg_rowSumsArray; /* the i's value is the sum of row i in B[g] */

} spmat;

spmat* AllocateSparseMatrix(int size, int nnz);
void AddRow(spmat *A, const int *neighbors, int rowSize, int i);
void AddRowList(spmat *A, ELEMENT *neighbors, int i);
void SparseWithVectorMult(spmat *A, double *v, double *result);
void BhatWithVectorMult(spmat *Bhat, double *vector, double *result,int isNormalized, int isRegularMult, double norma);
void freeSparseMatrix(spmat *A);


#endif
