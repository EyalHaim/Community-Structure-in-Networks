/**
 * "list" Summary:
 *
 * Module with list related data structures and their related functions.
 *
 * This Module contains the following functions:
 *
 * freeELEMENTList              - Free all the resources of ELEMENT list structure.
 * freeGROUPSnode               - Free all the resources of specific node in GROUPS structure.
 *
 */

#ifndef SPPROJECT_LIST_H
#define SPPROJECT_LIST_H

/**
 * Structure used to represent list of integers.
 */
typedef struct linked_list{
    int value;                /*value of a node - represents vertex of the graph*/
    int length;               /*the first node in each list will contain the length of the list*/
    struct linked_list *next;
} ELEMENT;


/**
 * List of lists structure. Will be used to represent the division of the graph into clusters.
 * each cluster (ELEMENT list) will be node in this structure.
 */
typedef struct list_of_lists{
    ELEMENT * head;
    struct list_of_lists *next;
} GROUPS;

/**
 * On a success, this function will free all the nodes of linked list which her head is headNode. If an error occurs,
 * nothing will happen.
 *
 * @param headNode            - Pointer to the head node of the list
 */
void freeELEMENTList(ELEMENT *headNode);

/**
 * On a success, this function will free all the resources of given node of GROUPS. If an error occurs, nothing will happen.
 *
 * @param node                - Pointer to the node.
 */
void freeGROUPSnode(GROUPS *node);

#endif
