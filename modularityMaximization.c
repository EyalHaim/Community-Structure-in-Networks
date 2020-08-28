#include "modularityMaximization.h"
#include "divisionUtils.h"
#include <stdlib.h>
#include <string.h>
#include "errorHandler.h"
#include <math.h>
#include <stdio.h>


void ModularityMaximization(spmat *Bhat, double *originalS){
    int *unMoved, *unMovedCurrIndex, *deg_array, *deg_arrayCurrIndex, *indices, *indicesCurrIndex;
    int subSizeGroup, outerIndex, innerIndex, maxScoreIndex, maxImproveIndex, M, i;
    double maxScore, currScore, deltaQ, maxImprove, currImprove, prevImprove, *sumsOfCols;
    double *Fig, *FigCurrIndex, *sumsOfColsCurrIndex, *originalSCurrIndex;

    subSizeGroup = Bhat -> n, M = Bhat -> M, deg_array = Bhat -> deg_array, Fig = Bhat -> Bg_rowSumsArray;
    InitVectorsForModularityMaximization(&unMoved, &sumsOfCols, &indices, subSizeGroup);

    do{
        memset(unMoved, 0, sizeof(int)*subSizeGroup);
        indicesCurrIndex = indices;
        maxImprove = -HUGE_VAL;
        for(outerIndex = 0; outerIndex < subSizeGroup; outerIndex++){
            BhatWithVectorMult(Bhat, originalS, sumsOfCols, 0,1,0);
            unMovedCurrIndex = unMoved, deg_arrayCurrIndex = deg_array, FigCurrIndex = Fig;
            originalSCurrIndex = originalS, sumsOfColsCurrIndex = sumsOfCols;
            maxScore = -HUGE_VAL;

            for(innerIndex=0; innerIndex <subSizeGroup; innerIndex++){
                currScore = 0.0;
                if(!(*unMovedCurrIndex)){                                 /*Check if current index is unmoved*/
                    currScore -= 4*(*originalSCurrIndex)*(*sumsOfColsCurrIndex);
                    currScore += 4*(0 - ((*deg_arrayCurrIndex)*(*deg_arrayCurrIndex)/M) - *FigCurrIndex);


                    if(currScore > maxScore) {
                        maxScore = currScore, maxScoreIndex = innerIndex;
                    }
                }
                 unMovedCurrIndex++, deg_arrayCurrIndex++, FigCurrIndex++;
                 originalSCurrIndex++, sumsOfColsCurrIndex++;
            }

            originalS[maxScoreIndex] *= -1;
            *indicesCurrIndex = maxScoreIndex;
            indicesCurrIndex++;

            currImprove = maxScore;
            if(outerIndex > 0)
                currImprove = prevImprove + maxScore;

            prevImprove = currImprove;
            if(currImprove > maxImprove) {
                maxImprove = currImprove, maxImproveIndex = outerIndex;
            }

            unMoved[maxScoreIndex] = 1;                                                   /*Moved relevant index*/
        }

        indicesCurrIndex = indices + (maxImproveIndex+1);
        for(i = maxImproveIndex + 1; i < subSizeGroup; i++){
            originalS[(*indicesCurrIndex)] *= -1;
            indicesCurrIndex++;
        }

        deltaQ = (maxImproveIndex == (subSizeGroup - 1) ? 0 : maxImprove);
    }while(IS_POSITIVE(deltaQ));

    FreeModularityMaximization(unMoved, sumsOfCols, indices);
}

void InitVectorsForModularityMaximization(int **unMoved, double **sumsOfCols,int **indices, int subSizeGroup){

    (*sumsOfCols) = (double *)malloc(sizeof(double)*subSizeGroup);
    EnsureMallocSucceeded(sumsOfCols);
    (*unMoved) = (int *)malloc(sizeof(int)*subSizeGroup);  /*0 if unmoved, 1 if moved*/
    EnsureMallocSucceeded(unMoved);
    (*indices) = (int *)malloc(sizeof(int)*subSizeGroup);
    EnsureMallocSucceeded(indices);
}

void FreeModularityMaximization(int *unMoved, double *sumsOfCols,int *indices){
    free(unMoved);
    free(sumsOfCols);
    free(indices);

}