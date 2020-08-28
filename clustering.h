#ifndef SPPROJECT_CLUSTERING_H
#define SPPROJECT_CLUSTERING_H
#include "list.h"
#include "spmat.h"

 GROUPS* FindOptimalDivsion(spmat *Amatrix, int *numOfClusters);
 void InitGroups(int size, GROUPS **headO,GROUPS **headP);
 int FindCluster(spmat *Amatrix, GROUPS *headP, GROUPS **headO);
 void DivideG(spmat *Bhat, ELEMENT *g, ELEMENT **g1, ELEMENT **g2);
 ELEMENT *BuildListNode(int *isFirst, ELEMENT **listHead, ELEMENT **listCurrent, ELEMENT *newNode);
 void UpdateListsByDivision(ELEMENT  *g, ELEMENT **g1,ELEMENT **g2, double *s);

#endif /*SPPROJECT_CLUSTERING_H*/
