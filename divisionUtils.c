#include "divisionUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include "errorHandler.h"
#include <time.h>
#include <math.h>


void CreateBgRowSumsArray(spmat *A_SubMatrix,int *rowptr, int M, int sumOfSubDegArray, int subGroupSize, double *rowSumsArray){
    /*Returns the Fig values for given sub matrix */
    int i, *subDegArray;
    subDegArray = A_SubMatrix -> deg_array;
    i = 0;
    A_SubMatrix -> Bg_rowSumsArray = rowSumsArray;
    for(i=0; i< subGroupSize; i++){
        *rowSumsArray = *(rowptr+1)-*(rowptr) - (double)(*subDegArray)*sumOfSubDegArray/M;
        rowSumsArray++;
        rowptr++;
        subDegArray++;
    }
}

int CreateSubDegArray(spmat *A_SubMatrix, ELEMENT *g_ListNodes, int *subDegArrayPtr, int *deg_array){
    /*Create deg array for given sub matrix */
    int sum;
    sum = 0;

    A_SubMatrix -> deg_array = subDegArrayPtr;
    while(g_ListNodes != NULL){
        *subDegArrayPtr = deg_array[g_ListNodes->value];
        sum += *subDegArrayPtr;
        g_ListNodes = g_ListNodes -> next;
        subDegArrayPtr++;
    }
    return sum;
}

int CalculateNNZofSubMatrix(ELEMENT *g_CurrentNode, int *rowptr){
    /* Calculate NNZ of given sub matrix by rowptr array*/

    int nnz, i;
    nnz = 0, i=0;

    while (g_CurrentNode != NULL) {
        if (i == g_CurrentNode->value) {
            nnz += *(rowptr + 1) - *(rowptr);
            g_CurrentNode = g_CurrentNode->next;
        }

        rowptr++;
        i++;
    }

    return nnz;
}

void CreateZeroMatrix(spmat *A_SubMatrix) {
    /* Creates zero matrix representation */

    free(A_SubMatrix -> colind);
    free(A_SubMatrix -> rowptr);
    A_SubMatrix -> rowptr = NULL;
    A_SubMatrix -> colind = NULL;
}

/*This function will build Submatrix of A according to the nodes mentioned in group g (Linked list).*/
spmat *BuildSubMatrix(spmat *Amatrix, ELEMENT *g_ListNodes) {
    int subGroupSize, i, nnz, *colind,row, col, isFirstAdded ,emptyMatAssertion, *rowptr, newCol;
    int *deg_array,*start, *startIndex, *end, *currDegArrayPtr, sumOfSubDegArray;
    double *Bg_rowSumsArray;
    ELEMENT *neighborsHead, *neighborsCurrent,  *g_CurrentNode, *g_tmpPointer;
    spmat * A_SubMatrix;
    colind = Amatrix->colind; rowptr = Amatrix->rowptr;
    deg_array = Amatrix -> deg_array; subGroupSize = g_ListNodes->length, g_CurrentNode = g_ListNodes;
    emptyMatAssertion = 0, i=0;                                                                                /*i will contain the current row in A_subMatrix*/
    nnz = CalculateNNZofSubMatrix(g_CurrentNode, rowptr);
    A_SubMatrix = AllocateSparseMatrix(subGroupSize, nnz);

    currDegArrayPtr = (int *) malloc(sizeof(int) * subGroupSize);
    EnsureMallocSucceeded(currDegArrayPtr);
    sumOfSubDegArray = CreateSubDegArray(A_SubMatrix, g_ListNodes, currDegArrayPtr, deg_array);  /*TODO SumSubDegArray = nnz!    */   /*sumOfSubDegArray will contain*/

    g_CurrentNode = g_ListNodes;
    rowptr = Amatrix->rowptr;
    while (g_CurrentNode != NULL) {                                                                              /*For each node in g, calculate his neighbors for the sub matrix Bg*/
        row = g_CurrentNode->value;
        startIndex = rowptr + row;
        start = colind + (*startIndex), end = colind + (*(startIndex + 1));                                      /*Start\end will contain the bounderis of currentNode in colind array = neighbors of node "row" in original graph*/
        g_tmpPointer = g_ListNodes;
        isFirstAdded = 1;
        neighborsHead = NULL, neighborsCurrent = NULL;                                                           /*neighborsHead will point to the head of the list of neighbors of currentNode according to g group*/
        newCol = 0;
        while (g_tmpPointer != NULL && start < end) {                                                            /*iterate over g to filter "g_currentNode" neighbors*/
            col = g_tmpPointer->value;
            if(col>*start)
                start++;
            else if (col<*start){
                g_tmpPointer = g_tmpPointer ->next;
                   newCol++;}
            else {
                if (!isFirstAdded) {
                    neighborsCurrent->next = (ELEMENT *) malloc(sizeof(ELEMENT));
                    EnsureMallocSucceeded((neighborsCurrent->next));
                    neighborsHead->length++, neighborsCurrent = neighborsCurrent->next;
                    neighborsCurrent->value = newCol;
                } else {
                    neighborsHead = (ELEMENT *) malloc(sizeof(ELEMENT));
                    EnsureMallocSucceeded(neighborsHead);
                    neighborsHead->value = newCol;
                    neighborsCurrent = neighborsHead, neighborsHead->length = 1;
                    isFirstAdded = 0;
                }
                start++;
                newCol++;
                g_tmpPointer = g_tmpPointer->next;
            }
        }

        AddRowList(A_SubMatrix, neighborsHead, i);
        if (neighborsHead != NULL) {
            neighborsCurrent->next = NULL;
            emptyMatAssertion += neighborsHead -> length;
            freeELEMENTList(neighborsHead);
        }


        g_CurrentNode = g_CurrentNode->next;
        i++;
    }

    Bg_rowSumsArray = (double *)malloc(sizeof(double )*subGroupSize);
    EnsureMallocSucceeded(Bg_rowSumsArray);
    CreateBgRowSumsArray(A_SubMatrix,A_SubMatrix->rowptr,Amatrix -> M, sumOfSubDegArray,subGroupSize, Bg_rowSumsArray);

    if(emptyMatAssertion == 0)
        CreateZeroMatrix(A_SubMatrix);

    A_SubMatrix -> M = Amatrix -> M;
    return A_SubMatrix;
}

double CalculateSumOfKiVi(double *vector, int *deg_array, int M, int subGroupSize){
    /* Calculate the sum of Vi*Ki for efficient Bhat multiply */

    int i;
    double sum;
    for(i=0;i<subGroupSize;i++){
        sum+= (*vector)*(*deg_array)/M;
        vector++;
        deg_array++;
    }

    return sum;
}

void PowerIteration(spmat *Bhat,double *oldVec, double *eigenVec, double *eigenVal, double norma){
    /*Update eigenVec and EigenVal of given sub matrix x*/

    int i, len, counter; /*TODO delete counter*/
    double *oldVecCurrIndex, *eigenVecCurrIndex, *tmp, *tmpVec;
    len = Bhat -> n;
    oldVecCurrIndex = oldVec, eigenVecCurrIndex = eigenVec;
    RandomizedVec(eigenVecCurrIndex, len);

    for(i=0; i<len; i++){                             /*Initiate oldVec*/
        *oldVecCurrIndex = *eigenVecCurrIndex;
        oldVecCurrIndex++; eigenVecCurrIndex++;
    }
    counter =0;
    do{                                                /*Power Iteration*/
        counter ++;
        tmp = oldVec;
        oldVec = eigenVec;
        eigenVec = tmp;
        BhatWithVectorMult(Bhat, oldVec, eigenVec, 1, 0, norma);
    } while(CheckEpsilon(oldVec, eigenVec, len));

    tmpVec = (double *)malloc(sizeof(double )*len);
    EnsureMallocSucceeded(tmpVec);
    *eigenVal = CalculateEigenValue(Bhat, eigenVec, norma, tmpVec);
    *eigenVal -= norma;
    free(tmpVec);
}

double CalculateEigenValue(spmat *Bhat, double *eigenVec, double norma, double *tempVec){
    double eigenValue;

    BhatWithVectorMult(Bhat, eigenVec, tempVec, 0, 0, norma);
    eigenValue = VecWithVecMult(eigenVec, tempVec, Bhat -> n);
    eigenValue /= VecWithVecMult(eigenVec, eigenVec, Bhat -> n);
    return eigenValue;
}

double VecWithVecMult(const double *vec1,const double *vec2, int len){
    int i;
    double sum;

    sum = 0;
    for(i=0; i< len; i++){
        sum += (*vec1)*(*vec2);
    }
    return sum;
}

void RandomizedVec(double *initialVec,int size){
    /* Insert random values to pre allocated array */

    int i;
    srand(time(NULL));

    for(i=0; i<size; i++)
    {
        *initialVec = rand();
        initialVec++;
    }
}

int CheckEpsilon(double *vecbOld, double *vecbNew, int size){
    /*Check if 2 vectors are close enough to stop the power iteration algorithm*/

    int i;
    for (i = 0; i<size; i++){
        if (fabs(*vecbNew - *vecbOld) >= EPSILON) {
            return 1;
        }
        vecbOld++;
        vecbNew++;
    }

    return 0;
}

double CalculateNorma(spmat *Bhat){
    int size, i, j;
    double *vector, *resultVec, colSum, maxSum, *vectorCurrIndex, *resultCurrIndex;
    size = Bhat -> n;
    maxSum = 0.0;
    vector = (double *)calloc(size, sizeof(double));                        /*Initiate 0's vector*/
    EnsureMallocSucceeded(vector);
    resultVec = (double  *)malloc(size*sizeof(double ));
    EnsureMallocSucceeded(resultVec);
    vectorCurrIndex = vector;

    for(i = 0; i< size; i++){                                               /*iterate over each column and sums its absolute values*/
        *vectorCurrIndex = 1.0;
        BhatWithVectorMult(Bhat, vector, resultVec, 0, 1, 0);
        colSum = 0.0;
        resultCurrIndex = resultVec;
        for(j = 0; j< size; j++){
            colSum += fabs(*resultCurrIndex);
            resultCurrIndex++;
        }

        if(colSum > maxSum)
            maxSum = colSum;

        *vectorCurrIndex = 0.0;
        vectorCurrIndex++;
    }

    free(vector);
    free(resultVec);
    return maxSum;
}

void CalculateSvector(double *S_vec, double *vector, int size, int isOnesVector){
    /*isOnesVector - 1 to initiate all vector with 1's, 0 to initiate vector by eigenVec*/
    int i;
    double value;
    for(i = 0; i< size; i++){
           value = *vector;
           if(IS_POSITIVE(value) || isOnesVector){
               *S_vec = 1;
           }
           else{
               *S_vec = -1;
           }
           vector++, S_vec++;
    }
}






