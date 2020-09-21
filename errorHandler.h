/**
 * "errorHandler" Summary:
 *
 * Module with the main macros of error handlers.
 *
 * This Module contains the following macros:
 *
 * IOFailure                    - Prints IO Failure
 * DivideByZeroError            - Prints divide by zero error (M=0)
 * ZeroNodesError               - Prints error if graph has no nodes
 * EnsureReadingSucceeded       - Prints error if reading from file has failed
 * EnsureWritingSucceeded       - Prints error if writing to file has failed
 * EnsureMallocSucceeded        - Prints error if memory allocation has failed
 * InfinityLoopError            - Prints error if infinite loop detected
 */

#ifndef SPPROJECT_ERRORHANDLER_H
#define SPPROJECT_ERRORHANDLER_H
#define EPSILON 0.00001
#define IS_POSITIVE(X) ((X) > 0.00001)


/**
 * Prints ERROR if opening the input file or output file has failed.
 */
#define IOFailure() do {printf("%s", "Error: IO Failure!\n"); exit(1);} while(0)

/**
 * Prints ERROR if the number of edges in the graph is 0, which will cause to division by zero.
 */
#define DivideByZeroError() do{printf("%s", "Error: Divide by zero!\n"); exit(1);} while(0)

/**
 * Prints ERROR if the graph contains no nodes.
 */
#define ZeroNodesError() do{printf("%s", "Error: Graph has no nodes!\n");exit(1);} while(0)

/**
 * Prints ERROR if reading bytes from input file has failed.
 */
#define EnsureReadingSucceeded(i, expected) do{if ((i)!= (expected)){\
            printf("%s", "Error: Read From File Failure!\n"); \
            exit(1);\
            }}while(0)

/**
 * Prints ERROR if writing bytes to output file has failed.
 */
#define EnsureWritingSucceeded(i, expected) do{if ((i)!= (expected)){\
            printf("%s", "Error: Write To File Failure!\n"); \
            exit(1);\
            }}while(0)

/**
 * Prints ERROR if dynamic memory allocation has failed.
 */
#define EnsureMallocSucceeded(ptr) do{if((ptr)==NULL){ \
            printf("%s", "Error: Failed Malloc!\n"); \
            exit(1); \
            }}while(0)

/**
 * Prints ERROR if infinite loops is detected
 */
#define InfinityLoopError() do{printf("%s", "Error: Infinite loop detected! \n"); exit(1);} while(0)

#endif
