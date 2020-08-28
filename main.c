#include "spmat.h"
#include "IO.h"
#include <stdlib.h>
#include <stdio.h>
#include "clustering.h"
#include <time.h> /*todo delete*/

int main(int argc, char* argv[]) {
    int numOfClusters;
    spmat *Amatrix;
    char *inputName, *outputName;
    GROUPS *clusteredGroupsHead;
    clock_t start,end; /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
    double timeall;

    start = clock();
    if(argc!=3) {
        printf("%s", "File Is Missing Error");
        exit(1);
    }

    inputName = argv[1];
    outputName = argv[2];
    ReadGraph(&Amatrix, inputName); /*Read Graph*/
    clusteredGroupsHead = FindOptimalDivsion(Amatrix, &numOfClusters);                          /*Find division of graph*/
    PrintDivision(clusteredGroupsHead, outputName, numOfClusters);
    Amatrix -> Bg_rowSumsArray = NULL; /*Non sub matrix */
    freeSparseMatrix(Amatrix);

    end = clock();
    timeall = ((double)(end-start) / CLOCKS_PER_SEC);
    printf("\n Our RunTime: %f", timeall);

    return 0;

}
