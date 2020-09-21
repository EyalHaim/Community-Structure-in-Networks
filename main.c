#include "graph.h"
#include "IO.h"
#include <stdlib.h>
#include <stdio.h>
#include "division.h"
#include <time.h>

/*Find Community Structure of given input graph*/
int main(int argc, char* argv[]) {

    int numOfClusters;
    mat *Amatrix;
    char *inputName, *outputName;
    GROUPS *clusteredGroupsHead;

    if(argc!=3) {
        printf("%s", "Error: Argument is missing! \n");
        exit(1);
    }

    inputName = argv[1];
    outputName = argv[2];

    srand(time(NULL));

    ReadGraph(&Amatrix, inputName);
    clusteredGroupsHead = FindOptimalDivision(Amatrix, &numOfClusters);       /*Find the optimal division of the graph*/
    PrintDivision(clusteredGroupsHead, outputName, numOfClusters);            /*Creates Output file */

    Amatrix -> Bg_rowSumsArray = NULL;                                        /*Amatrix is the original matrix and not a sub matrix */
    FreeMatrix(Amatrix, 1);

    return 0;                                                                 /*SUCCESS*/
}
