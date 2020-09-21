/**
 * "IO" Summary:
 *
 * Module with the IO related functions
 *
 * This Module contains the following functions:
 *
 * ReadGraph                    - Reads the input graph file into sparse matrix representation
 * PrintDivision                - Prints the Community Structure the program has found into output file
 *
 */

#ifndef SPPROJECT_IO_H
#define SPPROJECT_IO_H

/**
 * On a success, the function translates the input graph into sparse matrix representation. If an error
 * occurs, a relevant error will be printed.
 *
 * @param Amatrix               - Double pointer to mat structure. Will contain the input graph representation.
 * @param inputName             - File name of the input file.
 */
void ReadGraph(mat **Amatrix, char *inputName);

/**
 * On a success, the function will create output file based on the Community Structure that the program has found. If an
 * error occurs, a relevant error will be printed.
 *
 * @param headO                - Pointer to the head of list that contains the final optimal division.
 * @param outputName           - File name of the output file.
 * @param numOfClusters        - Number of clusters in the final optimal division.
 */
void PrintDivision(GROUPS *headO, char *outputName, int numOfClusters);


#endif
