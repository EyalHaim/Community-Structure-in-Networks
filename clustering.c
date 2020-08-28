#include "list.h"
#include "spmat.h"
#include <stdio.h>
#include <stdlib.h>
#include "clustering.h"
#include "errorHandler.h"
#include "divisionUtils.h"
#include "modularityMaximization.h"

GROUPS *FindOptimalDivsion(spmat *Amatrix, int *numOfClusters) {
    /*returns the number of clusters*/
    GROUPS *headP;
    GROUPS *headO;
    InitGroups(Amatrix->n, &headO, &headP);
    *numOfClusters = FindCluster(Amatrix, headP, &headO);
    return headO;
}


/*This function will initialize 2 Groups O and P.
  P will contain only 1 group with all the nodes of the graph.
  O will will empty.*/
void InitGroups(int size,GROUPS **headO, GROUPS **headP){
    int i;
    ELEMENT *tempHead;
    ELEMENT *temp, *node;
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

int FindCluster(spmat *Amatrix, GROUPS *headP, GROUPS **headO){
    /*returns the number of clusters*/
    int counter, isFirst;
    ELEMENT  *g1, *g2, *g;
    GROUPS  *node1, *node2, *toFree, *headOCurr;
    spmat *Bhat;

    counter = 0, isFirst = 1;
    while (headP != NULL){                                 /*REPEAT until P is empty. O pointer will contain the optimal Division at the end of the loop*/
        g = headP -> head;                                /*"remove" group G from P*/
        toFree = headP;
        headP = headP -> next;
        free(toFree);
        g1 = NULL, g2 = NULL;
        Bhat = BuildSubMatrix(Amatrix, g);                 /*Bhat will be represented with collection of elements in a struct?*/
        DivideG(Bhat,g, &g1, &g2);

        node1 = (GROUPS*)malloc(sizeof(GROUPS));
        EnsureMallocSucceeded(node1);
        if(g1 == NULL || g2 == NULL){
            node1->head = g; node1->next = NULL;
            if(isFirst){
                *headO = node1;
                headOCurr = *headO;
                isFirst = 0;}
            else{
                headOCurr -> next = node1;
                headOCurr = headOCurr->next;}
            counter++;}
        else {
            node2 = (GROUPS*)malloc(sizeof(GROUPS));
            EnsureMallocSucceeded(node2);
            node1->head = g1; node1 ->next = NULL;
            node2 -> head = g2; node2 -> next = NULL;
            if(g1->length == 1){
                if(isFirst){
                    *headO = node1;
                    headOCurr = *headO;
                    isFirst = 0;}
                else {
                    headOCurr->next = node1;
                    headOCurr = headOCurr->next;}
                counter ++;}
            else{
                node1->next = headP;
                headP = node1;}
            if(g2->length == 1){
                if(isFirst){
                    *headO = node2;
                    headOCurr = *headO;
                    isFirst = 0;}
                else {
                    headOCurr->next = node2;
                    headOCurr = headOCurr->next;}
                counter ++;}
            else{
                node2->next = headP;
                headP = node2;}
        }
        freeSparseMatrix(Bhat);
    }
    return counter;
}

void DivideG(spmat *Bhat, ELEMENT *g, ELEMENT **g1Ptr, ELEMENT **g2Ptr) {
    int subSizeGroup;
    double *eigenVec, *initialVec, *S_vec, *result, eigenVal, norma,Q0;

    subSizeGroup = Bhat -> n;
    eigenVec = (double *)malloc(sizeof(double ) *subSizeGroup);
    EnsureMallocSucceeded(eigenVec);
    S_vec = (double *)malloc(sizeof(double ) *subSizeGroup);
    EnsureMallocSucceeded( S_vec);
    initialVec = (double *) malloc(sizeof(double)*subSizeGroup);
    EnsureMallocSucceeded( initialVec);

    norma = CalculateNorma(Bhat);
    PowerIteration(Bhat, initialVec, eigenVec, &eigenVal, norma);

    if(eigenVal < 0){
        CalculateSvector(S_vec, eigenVec, subSizeGroup, 1);
    }
    else{
        result = (double *)malloc(sizeof(double ) *subSizeGroup);
        EnsureMallocSucceeded(result);

        CalculateSvector(S_vec, eigenVec, subSizeGroup, 0);
        BhatWithVectorMult(Bhat,S_vec,result, 0, 1, 0);
        Q0 = VecWithVecMult(S_vec, result, subSizeGroup);
        free(result);
        if(!IS_POSITIVE(Q0))
            CalculateSvector(S_vec, eigenVec, subSizeGroup, 1);
    }
    ModularityMaximization(Bhat, S_vec);

    UpdateListsByDivision(g, g1Ptr, g2Ptr, S_vec);


    free(eigenVec);
    free(S_vec);
    free(initialVec);

}

void UpdateListsByDivision(ELEMENT  *g, ELEMENT **g1,ELEMENT **g2, double *s){
    /*Divide g into g1,g2 by the division vector s*/

    int isFirst_g1, isFirst_g2;
    ELEMENT *g2CurrIndex, *g1CurrIndex, *g2Head, *g1Head, *gCurrIndex;
    double *sCurrIndex;

    gCurrIndex = g, g1Head = *g1, g2Head = *g2;
    sCurrIndex = s;
    isFirst_g1 = 1, isFirst_g2 =1;
    while(gCurrIndex != NULL){
        if( *sCurrIndex == 1)
           gCurrIndex = BuildListNode(&isFirst_g1, &g1Head, &g1CurrIndex, gCurrIndex);
        else
            gCurrIndex = BuildListNode(&isFirst_g2, &g2Head, &g2CurrIndex, gCurrIndex);

        sCurrIndex++;
    }
    *g1 = g1Head;
    *g2 = g2Head;

    if(isFirst_g1)
        *g1 = NULL;
    if (isFirst_g2)
        *g2 = NULL;

}

ELEMENT *BuildListNode(int *isFirst, ELEMENT **listHead, ELEMENT **listCurrent, ELEMENT *newNode){
   ELEMENT *nextIndex;
    if(*isFirst){
        *listHead = newNode;
        (*listHead) -> length = 1;
        *listCurrent = *listHead;
        nextIndex = (*listCurrent) -> next;
        (*listCurrent) -> next = NULL;
        *isFirst = 0;
    }

    else{
        (*listCurrent) -> next = newNode;
        (*listHead) -> length++;
        (*listCurrent) =  (*listCurrent)->next;
        nextIndex = (*listCurrent) -> next;
        (*listCurrent) -> next = NULL;
    }

    return nextIndex;
}


