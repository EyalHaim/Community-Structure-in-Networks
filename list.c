#include "list.h"
#include <stdlib.h>

void freeELEMENTList(ELEMENT *head){
    /*free the node of the list*/
    ELEMENT *current, *next;
    current = head;
    while(current){
        next = current -> next;
        free(current);
        current = next;
    }
}

void freeGROUPSnode(GROUPS *node){
    freeELEMENTList(node->head);
    free(node);
}
