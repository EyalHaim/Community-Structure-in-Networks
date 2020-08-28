#ifndef SPPROJECT_LIST_H
#define SPPROJECT_LIST_H

typedef struct linked_list{
    int value;
    int length;
    struct linked_list *next;
} ELEMENT;


/*Will be used to represent the set of Groups in P and the set of Groups in O.*/
/*Each Group is a linked list represented by ELEMENT nodes.*/
typedef struct list_of_lists{
    ELEMENT * head;
    struct list_of_lists *next;
} GROUPS;


void freeELEMENTList(ELEMENT *);
void freeGROUPSnode(GROUPS *node);

#endif /*SPPROJECT_LIST_H*/
