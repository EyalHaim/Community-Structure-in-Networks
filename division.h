/**
 * "division" Summary:
 *
 * Module with the main functions of the division algorithm, and their utils.
 *
 * This Module contains the following functions:
 *
 * FindOptimalDivision          - Creates the community structure of the graph
 * InitGroups                   - Initiate the algorithm with default division
 * FindDivision                 - Manages "Algorithm 3" - The main algorithm
 * DivideG                      - Manages "Algorithm 2" - Divide graph into 2 groups
 * ModularityMaximization       - Maximize the modularity of given division
 * UpdateSMultBArray            - Update SMultB array according to the moved vertex
 * CreateBgRowSumsArray         - Creates a vector with the sum of each column in B matrix(f_ig vector)
 * CreateSubDegArray            - Creates a deg array for sub-graph
 * CalculateNNZofSubMatrix      - Calculate the NNZ of given sub-graph
 * CreateZeroMatrix             - Creates a representation of the 0's matrix
 * BuildBhatSubMatrix           - Creates representation of sub-graph
 * PowerIteration               - Calculates the eigenvector and eigenvalue of given matrix.
 * CalculateEigenValue          - Calculates the eigenvalue of given matrix and given eigenvector
 * VecWithVecMult               - Calculates the multiplication of 2 vectors
 * RandomizedVec                - Creates a vector with randomized values
 * CheckEpsilon                 - Check if the "distance" between 2 vectors is less than epsilon
 * CalculateSvector             - Updates a vector to 1 or -1 according the positivity of his values
 * UpdateListsByDivision        - Creates a list of nodes which representing the division
 * BuildListNode                - Sub function used in UpdateListsByDivision to build node
 * AllocateFindDivisionVecs     - Allocate vectors that will be used in FindDivision function
 * FreeFindDivisionVecs         - Free the allocated dynamic memory used in FindDivision function
 */

#ifndef SPPROJECT_DIVISION_H
#define SPPROJECT_DIVISION_H
#include "list.h"
#include "graph.h"
#include "errorHandler.h"

/**
 * On success, the function returns the optimal division (Community Structure) of the graph.
 * If an error occurs, the function does nothing.
 *
 * @param Amatrix               - Representation of the input graph as Sparse Matrix
 * @param numOfClusters         - The function update this variable with the final number of clusters
 * @return
 * Return list of lists that represents the optimal division represented by headO
 */
GROUPS* FindOptimalDivision(mat *Amatrix, int *numOfClusters);

/**
 * On success, the function initialize O and P groups by updating headP to be list with 1 node that contains a list
 * with all the graph nodes, and headO to be Null pointer. If an error occurs, the function does nothing.
 *
 * @param size                  - Number of nodes in the graph
 * @param headO                 - Pointer to headO
 * @param headP                 - Pointer to headP
 */
void InitGroups(int size, GROUPS **headO,GROUPS **headP);

/**
 * On success, the function updates headO to represent the optimal division of the graph
 * and returns the number of clusters the algorithm has found. If an error occurs, a relevant error will be printed.
 *
 * @param Amatrix              - Representation of the input graph as sparse matrix
 * @param headP                - Pointer to headP
 * @param headO                - Pointer to pointer of headO
 * @return
 * Returns the number of groups in the division that the algorithm has found
 */
int FindDivision(mat *Amatrix, GROUPS *headP, GROUPS **headO);

/**On success, the function will update g1,g2 to be the optimal 2-group division of the sub-graph that declared
 * by the list of nodes - g. If an error occurs, a relevant error will be printed.
 *
 * @param Bhat                 - Pointer to mat variable that represents Bhat matrix of a subgraph (which built according to 'g').
 * @param g                    - List of nodes that represents sub group of vertices.
 * @param g1                   - Will contain the nodes of one of the groups after the division of g.
 * @param g2                   - Will contain the nodes of one of the groups after the division of g.
 * @param helpVec              - Pre-allocated vector that used for calculation (prevents unnecessary allocations).
 * @param SMultB               - Pre-allocated vector that used for modularity maximization algorithm.
 * @param unMoved              - Pre-allocated vector that used for modularity maximization algorithm.
 * @param indices              - Pre-allocated vector that used for modularity maximization algorithm.
 */
void DivideG(mat *Bhat, ELEMENT *g, ELEMENT **g1, ELEMENT **g2, double *helpVec, double *SMultB,
             int *unMoved, int *indices);

/**
 * On a success, this function will update originalS vector (which represent the division based on the leading
 * eigenvector of the sub-graph) according to the algorithm result.
 *
 * @param Bhat                  - Pointer to mat variable that represents Bhat matrix of a subgraph
 * @param originalS             - Represents the division based on the leading eigenvector of the sub-graph
 * @param unMoved               - Pre-allocated vec that declares if vertex was already moved (1) or unmoved (0).
 * @param SmultB                - Pre-allocated vec which his i's value contains the multiplication result of column i of mat B with vector s.
 * @param indices               - Pre-allocated vec which his i's values contains the vertex we moved in the i's loop.
 */
void ModularityMaximization(mat *Bhat, double *originalS, double *SMultB, int *Unmoved, int *indices);

/**
 * On a success, this function will update SMultB array according to the moved index.
 *
 * @param SMultBCurrIndex       - Current value of SMultB array
 * @param deg_arrayCurrIndex    - Current value of deg array
 * @param rowptr                - Array that used in a representation of a sparse matrix.
 * @param colind                - Array that used in a representation of a sparse matrix.
 * @param subGroupSize          - Number of vertices in the current sub-graph.
 * @param maxScoreIndex         - Index of the max score
 * @param Sj                    - contains 1/-1 according the vertex we moved
 * @param maxIndexDegDivideM    - The value of Kj / M where j is maxScoreIndex
 */
void UpdateSMultBArray(double *SMultBCurrIndex, int *deg_arrayCurrIndex, int *rowptr, int *colind, int subGroupSize, int maxScoreIndex, double Sj, double maxIndexDegDivideM);

/**
 * On a success, this function will updates Bg_rowSumsArray such that the i's value is the sum of row i in B[g]. If an
 * error occurs nothing will happen.
 *
 * @param subMatrix                - Pointer to mat variable that represents a subgraph of A.
 * @param M                        - Number of edges in the original input graph.
 * @param sumOfSubDegArray         - Number of edges in the sub graph.
 * @param subGroupSize             - Number of nodes in the sub graph.
 * @param Bg_rowSumsArray          - Update this array to contain values such that the i's value is the sum of row i in B[g]
 */
void CreateBgRowSumsArray(mat *subMatrix, int M, int sumOfSubDegArray, int subGroupSize, double *Bg_rowSumsArray);

/**
 * On a success, this function will updates submatDegArray to be the relevant degree's array according to the sub-graph.
 * If an error occurs nothing will happen.
 *
 * @param subMatrix                - Pointer to mat variable that represents a subgraph of A.
 * @param g_ListNodes              - List of nodes that represents sub group of vertices.
 * @param submatDegArray           - Vector that will contain the submat deg array.
 * @param deg_array                - Vertices degrees array of the original input graph.
 */
void CreateSubDegArray(mat *subMatrix, ELEMENT *g_ListNodes, int *submatDegArray, int *deg_array);

/**
 * On a success, this function will returns the NNZ of given sub-graph. If an error occurs nothing will happen.
 *
 * @param g_ListNodes              - List of nodes that represents sub group of vertices.
 * @param rowptr                   - Array that used in the representation of sparse matrix.
 * @return
 * Returns the calculated NNZ.
 */
int CalculateNNZofSubMatrix(ELEMENT *g_ListNodes, int *rowptr);

/**
 * On a success, this function will update subMatrix to represent a 0's matrix by initiating it's values to NULL.
 * If an error occurs nothing will happen.
 *
 * @param subMatrix             - Pointer to mat variable that represents a subgraph of A.
 */
void CreateZeroMatrix(mat *subMatrix);

/**
 * On a success, this function returns a pointer to mat struct that represents the B_g[hat] matrix according to current group g_ListNodes.
 * The relevant elements to represent the matrix are: Sub graph Ag(rowptr + colind), Vector Fig and the relevant deg array.
 * If an error occurs, nothing will happen.
 *
 * @param Amatrix                 - Representation of the input graph as Sparse Matrix
 * @param g                       - List of nodes that represents sub group of vertices.
 * @param submatDegArray          - Pre-allocated vector which will contain the deg_array of the built sub matrix.
 * @return
 * Returns pointer to the mat representation of the sub-graph.
 */
mat *BuildBhatSubMatrix(mat *Amatrix, ELEMENT *g, int *submatDegArray);

/**
 * On a success, this function will calculate eigen vector and eigen value of given B(hat) matrix.
 * If an error occurs nothing will happen.
 * @param Bhat                    - Pointer to mat that represents the subgraph.
 * @param initialVec              - Vector that used in the PowerIteration calculations
 * @param eigenVec                - Vector that will contain the values of the calculated eigen vector
 * @param eigenVal                - Will contain the eigen value of the shifted Bhat
 * @param norma                   - Contains the norma-1 of the first Bhat matrix
 * @param helpVec                 - Pre-allocated vector that used for calculations (prevents unnecessary allocations).
 */
void PowerIteration(mat *Bhat, double* initialVec, double *eigenVec, double *eigenVal, double norma, double *helpVec);

/**On a success, this function will return the eigenvalue of shifted Bhat, using eigenVec.
 * eigenValue = (eigenVec * [Bhat*eigenVec])/(eigenVec*eigenVec). If an error occurs nothing will happen.
 *
 * @param Bhat                   - Pointer to mat that represents the subgraph.
 * @param eigenVec               - Vector that will contain the values of the calculated eigen vector
 * @param norma                  - Contains the norma-1 of the first Bhat matrix
 * @param tempVec                - Pre-allocated vector which will contain the randomized values.
 * @return
 * Returns the eigenvalue of shifted Bhat.
 */
double CalculateEigenValue(mat *Bhat, double *eigenVec, double norma, double *tempVec);

/**
 * On a success, this function will return the multiplication result of vec1 and vec2. If an error occurs nothing will
 * happen.
 * @param vec1                   - Vector 1
 * @param vec2                   - Vector 2
 * @param size                   - Size of vector 1 and vector 2
 * @return
 * Returns the multiplication result
 */
double VecWithVecMult(double *vec1, double *vec2, int size);

/**
 * On a success, this function will update initialVec to contain randomized values. If an error occurs nothing will happen.
 *
 * @param initialVec               - Pre-allocated vector which will contain the randomized values.
 * @param size                     - Size of initialVec
 */
void RandomizedVec(double *initialVec,int size);

/**
 * On a success, this function will return 0 if |vecbOld[i] - vecbNew[i]| smaller than epsilon for each i. Otherwise, it will return 1.
 *  If an error occurs nothing will happen.
 *
 * @param vecbOld                  - Vector 1
 * @param vecbNew                  - Vector 2
 * @param size                     - Size of vector 1 and vector 2.
 * @return
 * 0 - if for each i - |vecbOld[i] - vecbNew[i]| < epsilon
 * 1 - Otherwise
 */
int CheckEpsilon(double *vecbOld, double *vecbNew, int size);

/**
 * On a success, this function will update s[i] to be 1 if s > 0, -1 otherwise, for each i. S will represent division.
 * If an error occurs nothing will happen.
 * @param s                     - Vector that will represent a division of a graph(-1/1)
 * @param size                  - Size of vector
 * @param isOnesVector          - 1 to initiate all vector with 1's, 0 to initiate vector by eigenVec (-1/1)
 */
void CalculateSvector(double *s, int size, int isOnesVector);

/**
 * On success, g1 and g2 will contain the vertices according to the division that the vector s represents.
 * If an error occurs, the function does nothing.
 *
 * @param g                -List of nodes that represents sub group of vertices.
 * @param g1               -Will contain the nodes of one of the groups after the division of g according to s vector.
 * @param g2               -Will contain the nodes of one of the groups after the division of g according to s vector.
 * @param s                -Vector of 1's and -1's that represents a division.
 */
void UpdateListsByDivision(ELEMENT  *g, ELEMENT **g1,ELEMENT **g2, double *s);

/**
 * A sub function that used in UpdateListsByDivision in order to add a node to g1/g2.
 * If an error occurs, the function does nothing
 * @param isFirst          -1 if newNode is the first node added to g1/g2, else 0
 * @param listHead         -Pointer to the head of g1/g1.
 * @param listCurrentNode  -Pointer to the last added Node in g1/g2
 * @param newNode          -New node to be added
 * @return
 * Returns a pointer to the recently added node
 */
ELEMENT *BuildListNode(int *isFirst, ELEMENT **listHead, ELEMENT **listCurrentNode, ELEMENT *newNode);

/**
 * On success, all the vectors used in  FindDivision will be allocated. If an error occurs,
 * relevant error will be printed.
 * @param submatDegArray   -Pre-allocated vector which will contain the deg_array of the built sub matrix.
 * @param helpVec          -Pre-allocated vector that used for calculation (prevents unnecessary allocations).
 * @param SMultB           -Pre-allocated vector that used for modularity maximization algorithm.
 * @param unMoved          -Pre-allocated vector that used for modularity maximization algorithm.
 * @param indices          -Pre-allocated vector that used for modularity maximization algorithm.
 * @param n                -Number of nodes in Amatrix.
 */
void AllocateFindDivisionVecs(int **submatDegArray, double **helpVec, double **SMultB, int **indices, int **unMoved, int n);

/**
 * On success, all the allocated vectors used in FindDivision will be freed.
 * @param submatDegArray   -Pre-allocated vector which will contain the deg_array of the built sub matrix.
 * @param helpVec          -Pre-allocated vector that used for calculation (prevents unnecessary allocations).
 * @param SMultB           -Pre-allocated vector that used for modularity maximization algorithm.
 * @param unMoved          -Pre-allocated vector that used for modularity maximization algorithm.
 * @param indices          -Pre-allocated vector that used for modularity maximization algorithm.
 */
void FreeFindDivisionVecs(int **submatDegArray, double **helpVec, double **SMultB, int **indices, int **unMoved);


#endif
