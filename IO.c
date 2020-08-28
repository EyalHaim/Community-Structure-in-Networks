#include <stdio.h>
#include "spmat.h"
#include "errorHandler.h"
#include <stdlib.h>
#include "IO.h"

void ReadGraph(spmat **Amatrix, char *inputName) {

    int read_assertion, node_index, n, nodeEdgesCount, M, size;
    int *neighbors_array, *k_index, *deg_array;
    FILE *inputFile;
    long numOfBytes;

    M = 0;
    inputFile = fopen(inputName, "r");                                                /*Open the input file*/
    if (inputFile == NULL)
        IOFailure();

    read_assertion = fread(&size, sizeof(int), 1, inputFile);                             /*Get the size of the graph*/

    /*TODO find out which error to show*/
    if (read_assertion == 0)                                                                  /*Case of empty Graph*/
        ZeroNodesError();

    EnsureReadingSucceeded(1, read_assertion);

    n = size;
    fseek(inputFile, 0,SEEK_END);                                                  /*TODO check stat.h*/
    numOfBytes = ftell(inputFile);
    M = ((numOfBytes - 4) - (n * 4)) / 4;
    deg_array = (int *) malloc(sizeof(int) *n);                                           /*deg_array will contain the degree of each node in the graph*/
    EnsureMallocSucceeded(deg_array);
    (*Amatrix) = AllocateSparseMatrix(n, M);
    k_index = deg_array;
    fseek(inputFile, 4,SEEK_SET);                                                       /*Move the file pointer to first node*/
    for (node_index = 0; node_index < n; node_index++) {                                        /*Create the sparse adjacency Matrix*/
        read_assertion = fread(&nodeEdgesCount, sizeof(int), 1, inputFile);
        EnsureReadingSucceeded(1, read_assertion);
        if (nodeEdgesCount > 0) {
            neighbors_array = (int *) malloc(sizeof(int) * nodeEdgesCount);
            EnsureMallocSucceeded(neighbors_array);
            read_assertion = fread(neighbors_array, sizeof(int), nodeEdgesCount, inputFile);
            EnsureReadingSucceeded(read_assertion, nodeEdgesCount);
            AddRow((*Amatrix), neighbors_array, nodeEdgesCount, node_index);
            free(neighbors_array);
            } else {
                AddRow((*Amatrix), NULL, 0, node_index);                           /*Case no neighbors*/
            }

            *k_index = nodeEdgesCount;
            k_index++;
        }

        if (M == 0)
            NoEdgesError();

        (*Amatrix)->M = M;
        (*Amatrix)->deg_array = deg_array;
        fclose(inputFile);


    /*     TESTPrintMatArrays(Amatrix,M);
        printf("\n***************************Done TestPrintMatArrays************************\n");
        TESTBgCreations(*Amatrix);
        printf("\n*****************************Done TESTBgCreations***************************\n");
        TESTMult(*Amatrix);
        printf("*****************************Done TESTMult**********************************\n");
        TESTBhatMult(*Amatrix);
        printf("*****************************Done TESTBHATMult**********************************\n");*//*TODO delete*/
}

void PrintDivision(GROUPS *headO, char *outputName, int numOfClusters){
    int n, length, value;
    FILE *outputFile;
    ELEMENT *list;
    GROUPS *toFree;

    outputFile = fopen(outputName, "w");
    if(outputFile == NULL)
        IOFailure();

    n = fwrite(&numOfClusters, sizeof(int),1,outputFile);           /*writes the header*/
    printf("\nnum of clusters is: %d\n", numOfClusters);
    EnsureWritingSucceeded(n, 1);
    while(headO != NULL){
        list = headO -> head;
        length = list -> length;

        n = fwrite(&length, sizeof(int),1,outputFile);           /*writes the length of group i*/
        printf("\n size of group: %d",length);
        EnsureWritingSucceeded(n, 1);
        printf("values are:");
        while(list != NULL){
            value = list->value;
            n = fwrite(&value, sizeof(int),1,outputFile);           /*writes the indices of group i*/
            printf(" %d,", value);
            EnsureWritingSucceeded(n, 1);
            list = list->next;
        }

        toFree = headO;
        headO = headO->next;
        freeGROUPSnode(toFree);

    }

    fclose(outputFile);
}
