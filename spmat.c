#include <stdio.h>
#include <stdlib.h>
#include "spmat.h"
#include <math.h>
#include "errorHandler.h"
#include <string.h>
#include "divisionUtils.h"

void AddRow(spmat *A, const int *neighbors, int rowSize, int i){          /*Used for building the sparse matrix*/
    /* adds a row to the sparse matrix */
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

void AddRowList(spmat *B, ELEMENT *neighbors, int i){                  /*Used for building the SubMatrix of A (Neighbors are ELEMENT List and not INT List)*/
    /* adds a row to the sparse matrix */
    /*if NULL row ptr i +1 = row ptr i => */
    int  j, rowSize;
    int *colind, *rowptr, *colindCurrIndex;

    colind = B -> colind;
    rowptr = B -> rowptr;
    if(neighbors != NULL){
        rowSize = neighbors ->length;
        colindCurrIndex = colind + (B -> lastInsertedIndex) + 1;

        for(j=0;j<rowSize;j++){
            *colindCurrIndex = neighbors -> value;
            colindCurrIndex++;
            neighbors = neighbors -> next;
        }

        if(i==0){*rowptr = 0;}
        *(rowptr+i+1) = *(rowptr + i) + rowSize;
        B->lastInsertedIndex += rowSize;

    } else{
        if(i==0){*rowptr = 0;}
        *(rowptr+i+1) = *(rowptr + i);
    }
}

void SparseWithVectorMult(spmat *A, double *v, double *result){
    /*multiply sparse matrix by given vector*/
    int len, startIndex, finishIndex, i, j, *colind, *rowptr;
    double tmpSum, *resultVec;

    colind = A->colind; rowptr = A->rowptr; len = A->n;

    if(rowptr == NULL){                                            /*Zero Matrix*/
        memset(result, 0, len*sizeof(double));
        return;
    }

    resultVec = result;
    for (i = 0;i<len; i++){
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

void BhatWithVectorMult(spmat *Bhat, double *vector, double *result,int isNormalized, int isRegularMult, double norma){
    /* isRegularMult - 1 if we multiply regular mat, 0 if we multiply shifted mat */
    /* isNormalized -  1 if we want to divide the result vec by his normalized value, 0 - otherwise */

    int subGroupSize, i, *deg_array;
    double sumOfKiVi ,*resultCurrIndex, *FigCurrIndex, sumNormalized;

    resultCurrIndex = result;
    FigCurrIndex = Bhat -> Bg_rowSumsArray;
    subGroupSize = Bhat ->n;
    deg_array = Bhat ->deg_array;
    SparseWithVectorMult(Bhat, vector, result);
    sumNormalized = 0.0;
    sumOfKiVi = CalculateSumOfKiVi(vector, Bhat->deg_array, Bhat->M, Bhat->n);          /*will be used for KiK/m mult */

    for(i = 0; i < subGroupSize; i++){
        *resultCurrIndex -= (*deg_array)*sumOfKiVi;                                     /* Minus the KiKv matrix part */
        *resultCurrIndex -= (*vector)*(*FigCurrIndex);                                  /* Minus the delta*Fig part */

        if(!isRegularMult){                                                             /* Plus the norma if relevant*/
            *resultCurrIndex += (*vector) * norma;
        }

        if(isNormalized){                                                               /*Will eventually be used to calc normalized vec*/
            sumNormalized += (*resultCurrIndex) * (*resultCurrIndex);
        }

        resultCurrIndex++;
        FigCurrIndex++;
        vector++;
        deg_array++;
    }

    if(isNormalized){
        resultCurrIndex = result;
        sumNormalized = sqrt(sumNormalized);
        for(i=0; i < subGroupSize; i++){                                                 /*Divide the vector by his normalized value */
            *resultCurrIndex /= sumNormalized;
            resultCurrIndex++;
        }
    }
}


void  freeSparseMatrix(struct _spmat *A){
    /*free the resources of the structure */
    if(A->colind != NULL)
         free(A->colind);
    if(A->rowptr != NULL)
        free(A->rowptr);
    if(A->deg_array != NULL)
        free(A->deg_array);
    if(A->Bg_rowSumsArray != NULL)
        free(A->Bg_rowSumsArray);
    free(A);
}

spmat* AllocateSparseMatrix(int n, int nnz){
    spmat *A;
    A = (spmat *)malloc(sizeof(spmat));
    EnsureMallocSucceeded(A);
    A -> colind = (int *)malloc(sizeof(int)*nnz);
    EnsureMallocSucceeded((A->colind));
    A -> rowptr = (int *)malloc(sizeof(int)*(n+1));
    EnsureMallocSucceeded((A->rowptr));
    A-> lastInsertedIndex = -1;
    A -> n = n;
    return A;
}

