#ifndef SPPROJECT_ERRORHANDLER_H
#define SPPROJECT_ERRORHANDLER_H

#define IOFailure() do {printf("%s", "IO Failure Error"); exit(1);} while(0)
#define NoEdgesError() do{printf("%s", "Divide by zero Error (M=0)"); exit(1);} while(0)
#define ZeroNodesError() do{printf("%s", "Graph has no Nodes Error");exit(1);} while(0)
#define EnsureReadingSucceeded(i, expected) do{if ((i)!= (expected)){\
            printf("%s", "Read From File Failure Error"); \
            exit(1);\
            }}while(0)
#define EnsureWritingSucceeded(i, expected) do{if ((i)!= (expected)){\
            printf("%s", "Write to File Failure Error"); \
            exit(1);\
            }}while(0)
#define EnsureMallocSucceeded(ptr) do{if((ptr)==NULL){ \
            printf("%s", "Failed Malloc Error"); \
            exit(1); \
            }}while(0)


#endif
