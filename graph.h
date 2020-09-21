/**
 * "mat" Summary:
 *
 * Module with functions and data structure related to A, B and Bhat matrices of a specific graph as defined in the algorithm instructions.
 *
 * A - Adjacency sparse matrix of the input graph.
 *   A[g] - Sub matrix of A based of subgroup of nodes - g.
 *
 * B - Modularity Matrix of the input graph.
 *   B[g] - Sub matrix of B based of subgroup of nodes - g.
 *
 * Bhat - Modularity Matrix of input graph where the sum of row Bi is subtracted from the [i,i] value.
 *   Bhat[g] - Sub matrix of Bhat base of subgroup of nodes - g.
 *
 * This Module contains the following functions:
 *
 * AddRowToSparseMatrix                     - Add an array as a row for A matrix
 * BWithVectorMult                          - Multiplication of B/Bhat/ShiftedBhat with a vector.
 * CalculateSumOfKiVi                       - Calculates the multiplication of deg array with given vector V.
 * SparseWithVectorMult                     - Multiplication of sparse matrix with a vector
 * CalculateNorma                           - Calculates the norma-1 of given matrix.
 * CalculateSumOfAbsValuesOfBhatCol         - Calculate the absolute sum of specific column in Bhat
 * AllocateMatrix                           - Allocate the resources used in mat structure
 * FreeMatrix                               - Free the resources used in mat structure.
 */

#ifndef _GRAPH_H
#define _GRAPH_H
#include "list.h"

/**
 * Structure that used to represent A, B and Bhat matrices of a specific graph as defined in the project's instructions (and explained above).
 * B_hat matrix - colind + rowptr (sparse matrix) + deg array + Fig array (a.k.a Bg_rowSumsArray)
 * A matrix     - sparse matrix. the value of Bg_rowSumsArray will be NULL
 * 0's matrix   - the value of rowptr and colind will be NULL
 */
typedef struct _mat {
    int      *colind;                           /*Array that used in the representation of a sparse matrix*/
    int      *rowptr;                           /*Array that used in the representation of a sparse matrix (contains the boundaries of indexes in colind array for each row)*/
    int      lastInsertedIndex;                 /*Contains the last used index in colind array while building a sparse matrix*/
    int	     n;                                 /*Number of nodes in the current graph*/
    int      M;                                 /*Number of edges in the original input graph*/
    int      *deg_array;                        /*Array that contains the degree of each node in the current graph */
    double   *Bg_rowSumsArray;                  /*Array such that the i's value is the sum of row i in B[g] */
    double   norma;                             /*norma-1 of the first created Bhat */

} mat;

/**
 * On a success, this function will add a row (that represented by neighbors array) to A matrix.
 * If an error occurs, nothing will happen.
 *
 * @param A                     - Pointer to sparse matrix representation of a graph.
 * @param neighbors             - Array of edges we will add to A.
 * @param rowSize               - The size of neighbors array.
 * @param rowIndex              - Row number.
 */
void AddRowToMatrix(mat *A, const int *neighbors, int rowSize, int rowIndex);

/**
 * On a success, this function will update result vector to contain the multiplication result of
 * B/Bhat matrix with specific vector. If and error occurs, nothing will happen.
 *
 * @param Bhat                   - Pointer to matrix representation of a graph.
 * @param vector                 - Vector.
 * @param result                 - Result vector.
 * @param isNormalized           - 1 if the result vec is divided by his normalized value, 0 - otherwise.
 * @param isShiftedMat           - 1 if Shifted mat is multiplied, 0 if regular mat is multiplied.
 * @param isBhat                 - 1 if Bhat is multiplied, 0 if B is multiplied.
 * @param norma                  - norma-1 of Bhat.
 */
void BWithVectorMult(mat *Bhat, double *vector, double *result,int isNormalized, int isShiftedMat, int isBhat, double norma);

/**
 * On a success, this function will return the multiplication result of deg array with given vector.
 * If an error occurs nothing will happen. The result will be used fot efficient B(hat) multiplication.
 *
 * @param vector                  - Given vector
 * @param deg_array               - Vertices degrees array of the original input graph.
 * @param M                       - Number of edges in the original input graph.
 * @param subGroupSize            - Number of nodes in the sub graph.
 * @return
 * Returns the multiplication's value.
 */
double CalculateSumOfKiVi(double *vector, int *deg_array, int M, int subGroupSize);

/**
 * On a success, this function will update result vector to contain the multiplication result
 * of Sparse matrix Ag (which is part of B[hat] representation) with vector v.
 * If an error occurs, nothing will happen. This function is used
 * as part of BWithVectorMult function
 * @param A                      - Pointer to sparse matrix representation of a graph.
 * @param v                      - Vector.
 * @param result                 - Result vector.
 */
void SparseWithVectorMult(mat *A, double *v, double *result);

/**
 * On a success, this function will return the norma-1 of Bhat matrix. If an error occurs nothing will happen.
 * @param Bhat                  - Pointer to mat that represents the subgraph
 * @param helpVec               - Pre-allocated vector that used for calculations (prevents unnecessary allocations).
 * @return
 * Returns the calculated norma
 */
double CalculateNorma(mat *Bhat, double *helpVec);

/**
 * Used as part of calculateNorma. On a success, this function will calculate the sum of Bhat specific col where we take
 * each value in his absolute value. If an error occurs, nothing will happen.
 *
 * @param Bhat                   - Pointer to matrix representation of a graph.
 * @param vector                 - Vector.
 * @param result                 - Result vector.
 * @param col                    - Column number.
 * @return
 *
 * Returns the sum of the absolute values of specific column in Bhat.
 */
double CalculateSumOfAbsValuesOfBhatCol(mat *Bhat, double *vector, double *result,int col);

/**
 * On a success, this function will allocate dynamic memory that used in mat structure.
 * If an error occurs, relevant error will be printed.
 *
 * @param size                  - Number of nodes
 * @param nnz                   - NNZ
 * @return
 * Returns pointer to the allocated structure.
 */
mat* AllocateMatrix(int size, int nnz);

/**
 * On a success, this function will free all the resources of given mat structure Matrix.
 * If an error occurs, nothing will happen.
 *
 * @param A                      - Pointer to matrix representation of a graph.
 * @param isOriginal             - 1 if we send original A matrix, 0 otherwise.
 */
void FreeMatrix(mat *A, int isOriginal);


#endif
