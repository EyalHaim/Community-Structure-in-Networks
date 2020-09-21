#include <stdio.h>
#include "graph.h"
#include "errorHandler.h"
#include <stdlib.h>
#include "IO.h"

/*Reads the input graph file into sparse matrix representation */
void ReadGraph(mat **Amatrix, char *inputName) {

    int read_assertion, node_index, n, nodeEdgesCount, M, size;
    int *neighbors_array, *k_index, *deg_array;
    FILE *inputFile;
    long numOfBytes;

    M = 0;
    inputFile = fopen(inputName, "rb");                                                 /*Open the input file*/
    if (inputFile == NULL)
        IOFailure();

    read_assertion = fread(&size, sizeof(int), 1, inputFile);                          /*Get the size of the graph*/


    if (read_assertion == 0)                                                                  /*Case of empty Graph*/
        ZeroNodesError();

    EnsureReadingSucceeded(1, read_assertion);

    n = size;
    fseek(inputFile, 0,SEEK_END);
    numOfBytes = ftell(inputFile);
    M = ((numOfBytes - 4) - (n * 4)) / 4;                                                     /*Calculate the number of edges*/
    deg_array = (int *) malloc(sizeof(int) *n);                                         /*deg_array will contain the degree of each node in the graph*/
    EnsureMallocSucceeded(deg_array);
    (*Amatrix) = AllocateMatrix(n, M);
    k_index = deg_array;
    fseek(inputFile, 4,SEEK_SET);                                                      /*Move the file pointer to first node*/


    neighbors_array = (int *) malloc(sizeof(int) * n);
    EnsureMallocSucceeded(neighbors_array);

    for (node_index = 0; node_index < n; node_index++) {                                      /*Create the sparse adjacency Matrix*/
        read_assertion = fread(&nodeEdgesCount, sizeof(int), 1, inputFile);
        EnsureReadingSucceeded(1, read_assertion);

        if (nodeEdgesCount > 0) {
            read_assertion = fread(neighbors_array, sizeof(int), nodeEdgesCount, inputFile);
            EnsureReadingSucceeded(read_assertion, nodeEdgesCount);
            AddRowToMatrix((*Amatrix), neighbors_array, nodeEdgesCount, node_index);
        }
        else
            AddRowToMatrix((*Amatrix), NULL, 0, node_index);                  /*Case no neighbors*/

        *k_index = nodeEdgesCount;
        k_index++;
    }
    free(neighbors_array);
    if (M == 0)
        DivideByZeroError();

    (*Amatrix)->M = M;
    (*Amatrix)->deg_array = deg_array;
    fclose(inputFile);
}

/*Prints the Community Structure the program has found into output file*/
void PrintDivision(GROUPS *headO, char *outputName, int numOfClusters){

    int n, length, value;
    FILE *outputFile;
    ELEMENT *list;
    GROUPS *toFree;

    outputFile = fopen(outputName, "wb");
    if(outputFile == NULL)
        IOFailure();

    n = fwrite(&numOfClusters, sizeof(int),1,outputFile);           /*writes the number of clusters*/

    EnsureWritingSucceeded(n, 1);
    while(headO != NULL){
        list = headO -> head;
        length = list -> length;

        n = fwrite(&length, sizeof(int),1,outputFile);              /*writes the length of group i*/
        EnsureWritingSucceeded(n, 1);

        while(list != NULL){
            value = list->value;
            n = fwrite(&value, sizeof(int),1,outputFile);           /*writes the indices of group i*/
            EnsureWritingSucceeded(n, 1);
            list = list->next;
        }

        toFree = headO;
        headO = headO->next;
        freeGROUPSnode(toFree);                                            /*free resources*/
    }

    fclose(outputFile);
}
