#include "division.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

/*Creates the community structure of the graph*/
GROUPS *FindOptimalDivision(mat *Amatrix, int *numOfClusters) {

    GROUPS *headP;
    GROUPS *headO;
    InitGroups(Amatrix->n, &headO, &headP);
    *numOfClusters = FindDivision(Amatrix, headP, &headO);
    return headO;
}

/*Initiate the algorithm with a trivial division*/
void InitGroups(int size, GROUPS **headO, GROUPS **headP){

    int i;
    ELEMENT *tempHead, *temp, *node;
    (*headP) = (GROUPS*)malloc(sizeof(GROUPS));
    EnsureMallocSucceeded((*headP));
    (*headP) -> next = NULL;

    /*Create list with all the nodes*/
    for(i=0;i<size;i++){
        node = (ELEMENT *)malloc(sizeof(ELEMENT));
        EnsureMallocSucceeded(node);
        node -> value = i;
        node -> next = NULL;
        if(i>0)
            temp -> next = node;
        else
            tempHead = node;
        temp = node;
    }

    tempHead -> length = size;
    (*headP) -> head = tempHead;
    (*headO) = NULL;
}

/*Manages "Algorithm 3" - The main algorithm*/
int FindDivision(mat *Amatrix, GROUPS *headP, GROUPS **headO){

    int numberOfGroups, isFirst, *submatDegArray, *unMoved, *indices, n, isOrigMat;
    double *helpVec, *SMultB, norma;
    ELEMENT  *g1, *g2, *g;
    GROUPS  *node1, *node2, *toFree, *headOCurr;
    mat *Bhat;

    numberOfGroups = 0, isFirst = 1, n= Amatrix->n, isOrigMat = 1;
    AllocateFindDivisionVecs(&submatDegArray, &helpVec, &SMultB, &indices, &unMoved, n);

    while (headP != NULL){                                              /*At the end of this loop, O pointer will contain the optimal Division*/
        g = headP -> head;                                              /*pull group g from P*/
        toFree = headP;
        headP = headP -> next;
        free(toFree);
        g1 = NULL, g2 = NULL;

        Bhat = BuildBhatSubMatrix(Amatrix, g, submatDegArray);          /*Bhat will be represented with collection of elements in a struct*/
        if(isOrigMat){                                                  /*Calculate only the norma of the first Bhat */
            norma = CalculateNorma(Bhat, helpVec);
            Bhat -> norma = norma, Amatrix ->norma = norma;
            isOrigMat = 0;}

        DivideG(Bhat, g, &g1, &g2, helpVec, SMultB, unMoved, indices);
        node1 = (GROUPS*)malloc(sizeof(GROUPS));
        EnsureMallocSucceeded(node1);
        if(g1 == NULL || g2 == NULL){                                   /*If g1 or g2 of size 0 add g to GROUP O */
            node1->head = g; node1->next = NULL;
            if(isFirst){
                *headO = node1, headOCurr = *headO;
                isFirst = 0;}
            else{
                headOCurr -> next = node1;
                headOCurr = headOCurr->next;}
            numberOfGroups++;}
        else {
            node2 = (GROUPS*)malloc(sizeof(GROUPS));
            EnsureMallocSucceeded(node2);
            node1->head = g1; node1 ->next = NULL;
            node2 -> head = g2; node2 -> next = NULL;
            if(g1->length == 1){                                        /*Add g1 to group O if size == 1*/
                if(isFirst){
                    *headO = node1, headOCurr = *headO;
                    isFirst = 0;}
                else {
                    headOCurr->next = node1;
                    headOCurr = headOCurr->next;}
                numberOfGroups ++;}
            else{                                                      /*Add g1 to group P if size > 1 */
                node1->next = headP;
                headP = node1;}
            if(g2->length == 1){                                       /*Add g2 to group O if size == 1 */
                if(isFirst){
                    *headO = node2, headOCurr = *headO;
                    isFirst = 0;}
                else {
                    headOCurr->next = node2;
                    headOCurr = headOCurr->next;}
                numberOfGroups ++;}
            else{                                                      /*Add g2 to P if size > 1 */
                node2->next = headP;
                headP = node2;}
        }
        FreeMatrix(Bhat,0);
    }

    FreeFindDivisionVecs(&submatDegArray, &helpVec, &SMultB, &indices, &unMoved);
    return numberOfGroups;
}

/*Manages "Algorithm 2" - Divide a group into 2 groups*/
void DivideG(mat *Bhat, ELEMENT *g, ELEMENT **g1Ptr, ELEMENT **g2Ptr, double *helpVec, double *SMultB, int *unMoved, int *indices) {
    int subGroupSize;
    double *eigenVec, *initialVec, eigenVal, norma,Q0;

    subGroupSize = Bhat -> n;
    eigenVec = (double *)malloc(sizeof(double ) *subGroupSize);             /*Pre allocated vec for Power Iteration*/
    EnsureMallocSucceeded(eigenVec);
    initialVec = (double *) malloc(sizeof(double)*subGroupSize);            /*Pre allocated vec for Power Iteration*/
    EnsureMallocSucceeded(initialVec);

    norma = Bhat ->norma;
    PowerIteration(Bhat, initialVec, eigenVec, &eigenVal, norma, helpVec);

    if(!IS_POSITIVE(eigenVal)){
        CalculateSvector(eigenVec, subGroupSize, 1);
    }
    else{
        CalculateSvector(eigenVec, subGroupSize, 0);
        BWithVectorMult(Bhat,eigenVec,helpVec, 0, 0,1, 0);
        Q0 = VecWithVecMult(eigenVec, helpVec, subGroupSize);
        if(!IS_POSITIVE(Q0))
            CalculateSvector(eigenVec, subGroupSize, 1);

        ModularityMaximization(Bhat, eigenVec, SMultB, unMoved, indices);
    }

    UpdateListsByDivision(g, g1Ptr, g2Ptr, eigenVec);

    free(eigenVec);
    free(initialVec);
}

/*Allocate vectors used in FindDivision function */
void AllocateFindDivisionVecs(int **submatDegArray, double **helpVec, double **SMultB, int **indices, int **unMoved, int n){
    *helpVec = (double *)malloc(sizeof(double)*n);                 /*Pre-allocated vector that used for future calculations */
    EnsureMallocSucceeded(helpVec);
    *submatDegArray= (int *)malloc(sizeof(int) * n);               /*Pre-allocated vec for building sub matrix*/
    EnsureMallocSucceeded(submatDegArray);
    *SMultB = (double *)malloc(sizeof(double)*n);                  /*Pre allocated vec for ModularityMaximization*/
    EnsureMallocSucceeded(SMultB);
    *unMoved = (int *)malloc(sizeof(int)*n);                       /*Pre allocated vec for ModularityMaximization*/
    EnsureMallocSucceeded(unMoved);
    *indices = (int *)malloc(sizeof(int)*n);                       /*Pre allocated vec for ModularityMaximization*/
    EnsureMallocSucceeded(indices);
}

/*Free the dynamic memory that is allocated in FindDivsion function */
void FreeFindDivisionVecs(int **submatDegArray, double **helpVec, double **SMultB, int **indices, int **unMoved){
    free(*submatDegArray);
    free(*helpVec);
    free(*SMultB);
    free(*indices);
    free(*unMoved);
}

/*Update eigenVec and EigenVal of given sub matrix x*/
void PowerIteration(mat *Bhat,double *oldVec, double *eigenVec, double *eigenVal, double norma, double *helpVec) {

    int i, len;
    long int counter, limit;
    double *oldVecCurrIndex, *eigenVecCurrIndex, *tmp;
    len = Bhat->n;
    limit = 6000*len;
    limit += 500000;
    RandomizedVec(eigenVec, len);
    oldVecCurrIndex = oldVec, eigenVecCurrIndex = eigenVec;

    for (i = 0; i < len; i++) {                                  /*Initiate oldVec*/
        *oldVecCurrIndex = *eigenVecCurrIndex;
        oldVecCurrIndex++;
        eigenVecCurrIndex++;
    }
    counter = 0;
    do {                                                         /*Power Iteration*/
        counter++;
        tmp = oldVec;
        oldVec = eigenVec;
        eigenVec = tmp;
        BWithVectorMult(Bhat, oldVec, eigenVec, 1, 1,1, norma);
    } while (CheckEpsilon(oldVec, eigenVec, len) && counter < limit);

    if(counter < limit) {
        if (counter % 2 == 1) {                                 /*Update eigenVec with the newest values*/
            oldVecCurrIndex = oldVec;
            eigenVecCurrIndex = eigenVec;
            for (i = 0; i < len; i++) {
                *eigenVecCurrIndex = *oldVecCurrIndex;
                eigenVecCurrIndex++;
                oldVecCurrIndex++;
            }
        }

        *eigenVal = CalculateEigenValue(Bhat, eigenVec, norma, helpVec);
        *eigenVal -= norma;
    }
    else{
        InfinityLoopError();
    }
}

/* Insert random values to pre allocated array */
void RandomizedVec(double *initialVec,int size){

    int i;
    for(i = 0; i < size; i++)
    {
        *initialVec = rand();
        initialVec++;
    }
}

/*Check if 2 vectors are close enough to stop the Power Iteration algorithm*/
int CheckEpsilon(double *vecbOld, double *vecbNew, int size){

    int i;

    for (i = 0; i<size; i++){
        if (fabs((*vecbNew) - (*vecbOld)) >= EPSILON)
            return 1;

        vecbOld++;
        vecbNew++;
    }

    return 0;
}

/*Calculate the eigenvalue of shifted Bhat*/
double CalculateEigenValue(mat *Bhat, double *eigenVec, double norma, double *tempVec){

    double eigenValue, tmp;

    BWithVectorMult(Bhat, eigenVec, tempVec, 0, 1,1, norma);
    eigenValue = VecWithVecMult(eigenVec, tempVec, Bhat -> n);
    tmp = VecWithVecMult(eigenVec, eigenVec, Bhat -> n);

    if(tmp > -EPSILON && tmp < EPSILON)
        DivideByZeroError();

    eigenValue = eigenValue/tmp;

    return eigenValue;
}

/*isOnesVector - 1 to initiate all vector with 1's, 0 to initiate vector by eigenVec*/
void CalculateSvector(double *s, int size, int isOnesVector){

    int i;
    double value;

    for(i = 0; i < size; i++){
        value = *s;
        if(isOnesVector || IS_POSITIVE(value)){
            *s = 1.0;
        }
        else {
            *s = -1.0;
        }

        s++;
    }
}

/*Calculates the multiplication result of 2 vectors*/
double VecWithVecMult(double *vec1, double *vec2, int len){

    int i;
    double sum;

    sum = 0;
    for(i = 0; i < len; i++){
        sum += (*vec1)*(*vec2);
        vec1++, vec2++;
    }

    return sum;
}

/*Maximize the modularity of given division*/
void ModularityMaximization(mat *Bhat, double *originalS, double *SMultB, int *unMoved, int *indices){

    int  *unMovedCurrIndex, *deg_array, *deg_arrayCurrIndex, *indicesCurrIndex, *rowptr, *colind;
    int subGroupSize, outerIndex, innerIndex, maxScoreIndex, maxImproveIndex, M, i;
    double maxScore, currScore, deltaQ, maxImprove, currImprove, prevImprove, Sj, maxIndexDegDivideM;
    double *Fig, *FigCurrIndex, *SMultBCurrIndex, *originalSCurrIndex;

    subGroupSize = Bhat -> n, M = Bhat -> M, deg_array = Bhat -> deg_array, Fig = Bhat -> Bg_rowSumsArray;

    do{
        memset(unMoved, 0, sizeof(int)*subGroupSize);
        indicesCurrIndex = indices;
        maxImprove = -HUGE_VAL;
        BWithVectorMult(Bhat, originalS, SMultB, 0,0,0,0);

        for(outerIndex = 0; outerIndex < subGroupSize; outerIndex++){
            unMovedCurrIndex = unMoved, deg_arrayCurrIndex = deg_array, FigCurrIndex = Fig, SMultBCurrIndex = SMultB;
            originalSCurrIndex = originalS;
            maxScore = -HUGE_VAL;

            for(innerIndex=0; innerIndex <subGroupSize; innerIndex++){
                if(!(*unMovedCurrIndex)){                                                                           /*Check if current index is unmoved*/
                    currScore = 0.0 - 4.0*(*originalSCurrIndex)*(*SMultBCurrIndex);
                    currScore -= 4.0*(((double)((*deg_arrayCurrIndex)*(*deg_arrayCurrIndex)))/(double)M);           /*Calculation of curr score*/

                    if(IS_POSITIVE((currScore) - (maxScore))) {
                        maxScore = currScore, maxScoreIndex = innerIndex;
                    }
                }
                unMovedCurrIndex++, deg_arrayCurrIndex++, FigCurrIndex++;
                originalSCurrIndex++, SMultBCurrIndex++;
            }

            Sj = originalS[maxScoreIndex]* -1;                                                                      /*"j" -> maxScoreIndex */
            originalS[maxScoreIndex] = Sj;
            colind = Bhat ->colind, rowptr = Bhat -> rowptr;
            maxIndexDegDivideM = (double)(deg_array[maxScoreIndex])/(double)M;
            SMultBCurrIndex = SMultB, deg_arrayCurrIndex = deg_array;

            UpdateSMultBArray(SMultBCurrIndex, deg_arrayCurrIndex, rowptr, colind, subGroupSize, maxScoreIndex, Sj,  /*Update SMultB array according to Sj*/
                              maxIndexDegDivideM);

            *indicesCurrIndex = maxScoreIndex;
            indicesCurrIndex++;

            currImprove = maxScore;
            if(outerIndex > 0)
                currImprove = prevImprove + maxScore;

            prevImprove = currImprove;
            if(IS_POSITIVE((currImprove) - (maxImprove)))
                maxImprove = currImprove, maxImproveIndex = outerIndex;

            unMoved[maxScoreIndex] = 1;                                                                              /*Move relevant index*/
        }

        indicesCurrIndex = indices + (maxImproveIndex+1);
        for(i = maxImproveIndex + 1; i < subGroupSize; i++){
            originalS[(*indicesCurrIndex)] *= -1;
            indicesCurrIndex++;}

        deltaQ = (maxImproveIndex == (subGroupSize - 1) ? 0 : maxImprove);

    } while(IS_POSITIVE(deltaQ));
}

/*Update SMultB array according to the moved vertex*/
void UpdateSMultBArray(double *SMultBCurrIndex, int *deg_arrayCurrIndex, int *rowptr, int *colind,
                       int subGroupSize, int maxScoreIndex, double Sj, double maxIndexDegDivideM) {

    int i, startIndex, finishIndex;
    double Ag_value;

    rowptr += maxScoreIndex;
    startIndex = *rowptr;
    finishIndex = *(rowptr + 1);
    colind += startIndex;

    for(i = 0; i < subGroupSize; i++) {
        Ag_value = 0.0;
        if (startIndex < finishIndex && i == (*colind)) {                                                        /*Extract the values of A[g] relevant row*/
            Ag_value = 1.0;
            startIndex++;
            colind++;
        }
        *SMultBCurrIndex += 2 * (Ag_value - ((double)(*deg_arrayCurrIndex) * maxIndexDegDivideM)) * (double)Sj;  /*Calculation to update SMultB after Sj has moved*/
        deg_arrayCurrIndex++, SMultBCurrIndex++;
    }
}

/*Creates a list of nodes which represents the division*/
void UpdateListsByDivision(ELEMENT  *g, ELEMENT **g1,ELEMENT **g2, double *s){
    /*Divide g into g1,g2 by the division vector s*/

    int isFirst_g1, isFirst_g2;
    ELEMENT *g2CurrIndex, *g1CurrIndex, *g2Head, *g1Head, *nodeToAdd;
    double *sCurrIndex;

    nodeToAdd = g, g1Head = *g1, g2Head = *g2;
    sCurrIndex = s;
    isFirst_g1 = 1, isFirst_g2 =1;

    while(nodeToAdd != NULL){
        if(*sCurrIndex == 1)
            nodeToAdd = BuildListNode(&isFirst_g1, &g1Head, &g1CurrIndex, nodeToAdd);
        else
            nodeToAdd = BuildListNode(&isFirst_g2, &g2Head, &g2CurrIndex, nodeToAdd);

        sCurrIndex++;
    }
    *g1 = g1Head;
    *g2 = g2Head;

    if(isFirst_g1)
        *g1 = NULL;
    if (isFirst_g2)
        *g2 = NULL;
}

/*Sub function used in UpdateListsByDivision to build node*/
ELEMENT *BuildListNode(int *isFirst, ELEMENT **listHead, ELEMENT **listCurrentNode, ELEMENT *newNode){
    ELEMENT *nextIndex;

    if(*isFirst){
        *listHead = newNode;
        (*listHead) -> length = 1;
        *listCurrentNode = *listHead;
        nextIndex = (*listCurrentNode) -> next;
        (*listCurrentNode) -> next = NULL;
        *isFirst = 0;
    }
    else{
        (*listCurrentNode) -> next = newNode;
        (*listHead) -> length++;
        (*listCurrentNode) =  (*listCurrentNode)->next;
        nextIndex = (*listCurrentNode) -> next;
        (*listCurrentNode) -> next = NULL;
    }

    return nextIndex;
}

/*Creates all the relevant elements required to represent B_g[hat]*/
mat *BuildBhatSubMatrix(mat *Amatrix, ELEMENT *g_ListNodes, int *submatDegArray) {

    int subGroupSize, subMatCurrRow, nnz,origMatCurrRow, col, emptyMatAssertion, length, newCol;
    int *deg_array, *start, *startIndex, *end, *neighbors, *neighborsCurrIndex, *colind, *rowptr;
    double *Bg_rowSumsArray;
    ELEMENT *g_CurrentNode, *g_tmpPointer;
    mat *subMatrix;

    colind = Amatrix->colind; rowptr = Amatrix->rowptr;
    deg_array = Amatrix -> deg_array; subGroupSize = g_ListNodes->length;
    g_CurrentNode = g_ListNodes, emptyMatAssertion = 0, subMatCurrRow=0;

    neighbors = (int *)malloc(sizeof(int)*subGroupSize);
    nnz = CalculateNNZofSubMatrix(g_CurrentNode, rowptr);
    subMatrix = AllocateMatrix(subGroupSize, nnz);

    CreateSubDegArray(subMatrix, g_ListNodes, submatDegArray, deg_array);               /*Creates the deg array for the relevant sub matrix */

    g_CurrentNode = g_ListNodes;
    rowptr = Amatrix->rowptr;
    while (g_CurrentNode != NULL) {                                                     /*For each node in g, calculate his neighbors for the sub matrix A[g]*/
        origMatCurrRow = g_CurrentNode->value;
        startIndex = rowptr + origMatCurrRow;
        start = colind + (*startIndex), end = colind + (*(startIndex + 1));             /*Start\end will contain the boundaries of currentNode in colind array = neighbors of node origMatCurrRow in original graph*/
        g_tmpPointer = g_ListNodes;
        newCol = 0, length = 0;
        neighborsCurrIndex = neighbors;

        while (g_tmpPointer != NULL && start < end) {                                   /*Iterate over g to filter "g_currentNode" neighbors*/
            col = g_tmpPointer->value;
            if(col>*start)
                start++;
            else if (col<*start){
                g_tmpPointer = g_tmpPointer ->next;
                newCol++;}
            else {
                *neighborsCurrIndex = newCol;
                length++;
                start++, newCol++;
                g_tmpPointer = g_tmpPointer->next;
                neighborsCurrIndex++;
            }
        }
        AddRowToMatrix(subMatrix, neighbors, length, subMatCurrRow);                    /*Add row to the sub matrix */
        emptyMatAssertion += length;

        g_CurrentNode = g_CurrentNode->next;
        subMatCurrRow++;
    }

    Bg_rowSumsArray = (double *)malloc(sizeof(double )*subGroupSize);
    EnsureMallocSucceeded(Bg_rowSumsArray);
    CreateBgRowSumsArray(subMatrix, Amatrix -> M, nnz, subGroupSize, Bg_rowSumsArray);       /*Creates the vector Fig  */
    if(emptyMatAssertion == 0)
        CreateZeroMatrix(subMatrix);

    subMatrix -> M = Amatrix -> M; subMatrix->norma = Amatrix->norma;
    free(neighbors);

    return subMatrix;
}

/*Computes the Fig vector (Bg Matrix row sums array) for given sub matrix */
void CreateBgRowSumsArray(mat *subMatrix, int M, int sumOfSubDegArray, int subGroupSize, double *Bg_rowSumsArray){

    int i, *subDegArray, *rowptr;

    rowptr = subMatrix->rowptr;
    subDegArray = subMatrix -> deg_array;
    subMatrix -> Bg_rowSumsArray = Bg_rowSumsArray;

    for(i=0; i< subGroupSize; i++){
        *Bg_rowSumsArray = *(rowptr+1)-*(rowptr) - (double)(*subDegArray)*(((double)sumOfSubDegArray)/M);
        Bg_rowSumsArray++;
        rowptr++;
        subDegArray++;
    }
}

/*Create deg array for given sub matrix */
void CreateSubDegArray(mat *subMatrix, ELEMENT *g_ListNodes, int *submatDegArray, int *deg_array){

    subMatrix -> deg_array = submatDegArray;

    while(g_ListNodes != NULL){
        *submatDegArray = deg_array[g_ListNodes->value];
        g_ListNodes = g_ListNodes -> next;
        submatDegArray++;
    }
}

/* Calculate NNZ of given sub matrix (use only her rowptr)*/
int CalculateNNZofSubMatrix(ELEMENT *g_ListNodes, int *rowptr){

    int nnz, i;
    nnz = 0, i=0;

    while (g_ListNodes != NULL) {
        if (i == g_ListNodes->value) {
            nnz += *(rowptr + 1) - *(rowptr);
            g_ListNodes = g_ListNodes->next;
        }
        rowptr++;
        i++;
    }
    return nnz;
}

/* Creates zero matrix representation */
void CreateZeroMatrix(mat *subMatrix) {

    free(subMatrix -> colind);
    free(subMatrix -> rowptr);
    subMatrix -> rowptr = NULL;     /*Representation of 0's matrix*/
    subMatrix -> colind = NULL;
}