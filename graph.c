#include <stdio.h>
#include <stdlib.h>
#include "graph.h"
#include <math.h>
#include "errorHandler.h"
#include <string.h>

/*Used for building a sparse matrix - add a row to the matrix */
void AddRowToMatrix(mat *A, const int *neighbors, int rowSize, int i){

    int  j;
    int *colind, *rowptr, *colindCurrIndex;

    colind = A -> colind;
    rowptr = A -> rowptr;
    colindCurrIndex = colind + (A -> lastInsertedIndex) + 1;

    for(j=0;j<rowSize;j++){
        *colindCurrIndex = *neighbors;
        colindCurrIndex++;
        neighbors++;
    }

    if(i==0)
        *rowptr = 0;

    *(rowptr+i+1) = *(rowptr + i) + rowSize;
    A->lastInsertedIndex += rowSize;
}

/*Multiplication of B/Bhat/ShiftedBhat with a vector*/
void BWithVectorMult(mat *Bhat, double *vector, double *result,int isNormalized, int isShiftedMat, int isBhat, double norma){
    /* isShiftedMat -  1 if we multiply shifted mat, 0 if we multiply regular mat*/
    /* isNormalized -  1 if we divide the result vec by his normalized value, 0 - otherwise */
    /* isBhat       -  1 if we multiply Bhat, 0 if we multiply B */

    int subGroupSize, i, *deg_array;
    double *resultCurrIndex, *FigCurrIndex, sumNormalized, sumOfKiVi;

    resultCurrIndex = result;
    FigCurrIndex = Bhat -> Bg_rowSumsArray;
    subGroupSize = Bhat ->n;
    deg_array = Bhat ->deg_array;

    SparseWithVectorMult(Bhat, vector, result);
    sumNormalized = 0.0;
    sumOfKiVi = CalculateSumOfKiVi(vector, Bhat->deg_array, Bhat->M, Bhat->n);          /* Will be used for KiK/m mult */

    for(i = 0; i < subGroupSize; i++){
        *resultCurrIndex -= (*deg_array)*sumOfKiVi;                                     /* Minus the KiKv matrix part */
        if(isBhat)
            *resultCurrIndex -= (*vector)*(*FigCurrIndex);                              /* Minus the delta*Fig part if Bhat*/

        if(isShiftedMat)                                                                /* Plus the norma if relevant*/
            *resultCurrIndex += (*vector) * norma;

        if(isNormalized)                                                                /* Will eventually be used to calc normalized vec*/
            sumNormalized += (*resultCurrIndex) * (*resultCurrIndex);

        resultCurrIndex++;
        FigCurrIndex++;
        vector++;
        deg_array++;
    }

    if(isNormalized){
        resultCurrIndex = result;
        sumNormalized = sqrt(sumNormalized);

        if(!IS_POSITIVE(sumNormalized))
            DivideByZeroError();

        for(i=0; i < subGroupSize; i++){                                                /*Divide the vector by his normalized value */
            *resultCurrIndex /= sumNormalized;
            resultCurrIndex++;
        }
    }
}

/*Multiply sparse matrix by given vector*/
void SparseWithVectorMult(mat *A, double *v, double *result){

    int len, startIndex, finishIndex, i, j, *colind, *rowptr;
    double tmpSum, *resultVec;

    colind = A->colind; rowptr = A->rowptr; len = A->n;

    if(rowptr == NULL){                                            /*0's matrix case - returns 0's vector*/
        memset(result, 0, len*sizeof(double));
        return;
    }

    resultVec = result;
    for (i = 0;i<len; i++){                                        /*Iterate over the matrix rows*/
        startIndex = *rowptr;
        rowptr++;
        finishIndex = *rowptr;
        tmpSum=0.0;

        for(j =startIndex; j<finishIndex;j++){
            tmpSum += *(v + *colind);
            colind++;
        }

        *resultVec = tmpSum;
        resultVec++;
    }
}

/* Calculate the sum of Vi*Ki/M for efficient Bhat multiply */
double CalculateSumOfKiVi(double *vector, int *deg_array, int M, int subGroupSize){

    int i;
    double sum;

    sum = 0.0;
    for(i = 0; i < subGroupSize; i++){
        sum+= (*vector)*(((double)(*deg_array))/M);
        vector++;
        deg_array++;
    }
    return sum;
}

/*Calculates the norma of given Bhat*/
double CalculateNorma(mat *Bhat, double *helpVec){

    int size, i;
    double *vector, colSum, maxSum, *vectorCurrIndex;

    size = Bhat -> n;
    maxSum = 0.0;
    vector = (double *)calloc(size, sizeof(double));                         /*Initiate 0's vector*/
    EnsureMallocSucceeded(vector);
    vectorCurrIndex = vector;

    for(i = 0; i < size; i++){                                               /*Iterate over each column and sums its absolute values*/
        *vectorCurrIndex = 1.0;
        colSum = CalculateSumOfAbsValuesOfBhatCol(Bhat, vector, helpVec, i);

        if(IS_POSITIVE((colSum) - (maxSum)))
            maxSum = colSum;

        *vectorCurrIndex = 0.0;
        vectorCurrIndex++;
    }

    free(vector);
    return maxSum;
}

/*Calculate the absolute sum of specific column in Bhat*/
double CalculateSumOfAbsValuesOfBhatCol(mat *Bhat, double *vector, double *result, int col){
    int subGroupSize, i, *deg_array, *rowptr, *colind, startIndex, finishIndex;
    double *resultCurrIndex, *FigCurrIndex, sum, sumOfKiVi;

    resultCurrIndex = result;
    FigCurrIndex = Bhat -> Bg_rowSumsArray; subGroupSize = Bhat ->n;
    deg_array = Bhat ->deg_array; colind = Bhat->colind; rowptr = Bhat->rowptr;

    rowptr += col;
    startIndex = *rowptr;
    rowptr++;
    finishIndex = *rowptr;
    colind += startIndex;

    sum = 0.0;
    sumOfKiVi = CalculateSumOfKiVi(vector, Bhat->deg_array, Bhat->M, Bhat->n);        /*Will be used for KiKj/m mult*/

    for(i = 0; i < subGroupSize; i++){
        *resultCurrIndex = 0;
        if(startIndex<finishIndex && i==(*colind)){                                   /*Extract the relevant values from A[g] */
            *resultCurrIndex = 1;
            startIndex++;
            colind++;
        }
        *resultCurrIndex -= (*deg_array)*sumOfKiVi;                                   /* Minus the KiKv matrix part */
        *resultCurrIndex -= (*vector)*(*FigCurrIndex);                                /* Minus the delta*Fig part */
        sum += fabs(*resultCurrIndex);

        resultCurrIndex++;
        FigCurrIndex++;
        vector++;
        deg_array++;
    }

    return sum;
}

/*Allocate the resources used in mat structure*/
mat* AllocateMatrix(int n, int nnz){
    mat *A;

    A = (mat *)malloc(sizeof(mat));
    EnsureMallocSucceeded(A);
    A -> colind = (int *)malloc(sizeof(int)*nnz);
    EnsureMallocSucceeded((A->colind));
    A -> rowptr = (int *)malloc(sizeof(int)*(n+1));
    EnsureMallocSucceeded((A->rowptr));
    A-> lastInsertedIndex = -1;
    A -> n = n;

    return A;
}

/*Free the resources used in mat structure*/
void FreeMatrix(struct _mat *A, int isOriginal){

    if(A->colind != NULL)
        free(A->colind);
    if(A->rowptr != NULL)
        free(A->rowptr);
    if(A->Bg_rowSumsArray != NULL)
        free(A->Bg_rowSumsArray);
    if(isOriginal && A->deg_array != NULL)
        free(A->deg_array);

    free(A);
}


