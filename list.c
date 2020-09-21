#include "list.h"
#include <stdlib.h>

/*Free all the resources of ELEMENT list structure*/
void freeELEMENTList(ELEMENT *head){

    ELEMENT *current, *next;

    current = head;
    while(current){
        next = current -> next;
        free(current);
        current = next;
    }
}

/*Free all the resources of specific node in GROUPS structure*/
void freeGROUPSnode(GROUPS *node){
    freeELEMENTList(node->head);
    free(node);
}
