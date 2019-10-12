/*  L.c
 
    List data structure and utility functions.  
    Elements can only be added or deleted from end
  
    To create a new linked list:
       L = NEW_L(N);
       where N is the number of nodes initially allocated

    Functions:
    NEW_L - make a new linked list
    L_GET_ID - get the id of the i-th node
    L_GET_VAL - get the value of the i-th node
    L_FIND - get the value in a node with specified id
    L_FIND_ID - get the id of a node with specified value
    L_ADD - add a node to the end of the list
    L_POP - remove last node from list (get id + val)
    L_TO_STRING - print out a list
    L_CLEAR - clear a linked list data structure from memory

    Internal Functions:
    double_size - double the size of the allocated array

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */

#include "L.h"
#include <stdio.h>
#include <stdlib.h>


void double_size(struct L * L)
{
    len_t i; //counter
    struct LNODE * nodes; //pointer to new larger node array
    nodes = (struct LNODE *) malloc(2*L->na*sizeof(struct LNODE)); //allocate node array with 2x size
    for (i=0; i<L->n; i++) { //copy old list into new
        nodes[i].id = L->nodes[i].id;
        nodes[i].v = L->nodes[i].v;
    }
    for (i=L->n; i<2*L->n; i++) { //initialize new empty nodes
        nodes[i].id = NULL_ID;
        nodes[i].v = NULL_VAL;
    }
    L->na = 2*L->na; //new number of allocated spots
    free(L->nodes); //free old smaller array
    L->nodes = nodes; //point to new larger arrray
}


struct L * NEW_L(len_t N)
{
    len_t i; //counter
    struct L * L; //pointer to list header
    struct LNODE * nodes; //pointer to node array
    if (N<1) { N = 1; } //give em something
    L = (struct L *) malloc(sizeof(struct L)); //allocate new linked list head
    nodes = (struct LNODE *) malloc(N*sizeof(struct LNODE)); //allocate node array
    for (i=0; i<N; i++){ //initialize node elements
        nodes[i].id = NULL_ID;
        nodes[i].v = NULL_VAL;
    }
    L->nodes = nodes; //assign that array to that list head
    L->n = 0; //store number of nodes
    L->na = N; //store number of allocated nodes
    return L; //return the new list
}


// Return the id of the i-th node
nid_t L_GET_ID(struct L * L, len_t ind)
{
    if (ind > L->n-1) { 
        return NULL_VAL; //ind is past end of list
    } else {
        return L->nodes[ind].id; //return the ind-th id
    }
}


// Return the value stored by the i-th node
val_t L_GET_VAL(struct L * L, len_t ind)
{
    if (ind > L->n-1) { 
        return NULL_VAL; //ind is past end of list
    } else {
        return L->nodes[ind].v; //return the ind-th value
    }
}


// Set the id of the i-th node
void L_SET_ID(struct L * L, len_t ind, nid_t id)
{
    if (ind < L->n-1) { 
        L->nodes[ind].id = id; //set the ind-th id
    }
}


// Set the value stored by the i-th node
void L_SET_VAL(struct L * L, len_t ind, val_t val)
{
    if (ind < L->n-1) { 
        L->nodes[ind].v = val; //set the ind-th value
    }
}


// Return the value stored by node with specified id
val_t L_FIND(struct L * L, nid_t id)
{
    int i; //counter
    for (i=0; i<L->n; i++) {
        if (L->nodes[i].id == id){ return L->nodes[i].v; } //return value if id is found
    }
    return NULL_VAL; //value was not found
}


// Return the index of a node with specified value
nid_t L_FIND_ID(struct L * L, val_t v)
{
    int i; //counter
    for (i=0; i<L->n; i++) {
        if (L->nodes[i].v == v){ return L->nodes[i].id; } //return id if val is found
    }
    return NULL_ID; //value was not found
}


// Add a node to the end of the list
struct LNODE * L_ADD(struct L * L, nid_t id, val_t v)
{
    if (L->n+1 > L->na) { double_size(L); } //double size if we need more space
    L->nodes[L->n].id = id; //assign id
    L->nodes[L->n].v = v; //assign value
    L->n++; //add to length of list
    return &L->nodes[L->n]; //return pointer to new node
}


// Add a value to the existing value stored by node with specified id (adds new node if not found)
int L_ADD_VAL_BY_ID(struct L * L, nid_t id, val_t v)
{
    int i; //counter and found flag
    for (i=0; i<L->n; i++) { //search thru list
        if (L->nodes[i].id == id){ //if node is found
            L->nodes[i].v += v; //add to that value
            return 0; //didn't have to add an extra node
        }
    } //node was not found, so
    L_ADD(L, id, v); //add new node with that value
    return 1; //added an extra node
}


// Get and remove a node from the end of the list
int L_POP(struct L * L, nid_t * id, val_t * v)
{
    if (!L){ return 0; } //can't pop a null pointer
    if (L->n < 1){ return 0; } //can't pop an empty list
    *id = L->nodes[L->n-1].id; //return this node's id
    *v = L->nodes[L->n-1].v; //return this node's value
    L->n--; //decrease list length counter
    return 1; //sucessfully popped
}


// Get the length of a list
len_t L_LENGTH(struct L * L)
{
    if (L) { //if it's not a null pointer
        return L->n;
    } else {
        return 0;
    }
}


// Print out a list
void L_TO_STRING(struct L * L)
{
    len_t i; //counter
    if (L->n < 1){
        printf("Empty List");
    } else { 
        for (i=0; i<L->n; i++) {
            printf("%d\t%f\n", L->nodes[i].id, L->nodes[i].v); //print value and id
        }
    } 
}


// Empty a list without de-allocating
void L_EMPTY(struct L * L)
{
    if (L) { //if L is a valid pointer
        L->n = 0; //set list to empty
    }
}


// Clear a list data structure from memory
void L_CLEAR(struct L * L)
{
    free(L->nodes); //free mem for array
    free(L); //free the list head
}


// Add a node to the end of the list
struct LNODE * L_PUSH_ID(struct L * L, nid_t id, nid_t id2)
{
    if (L->n+1 > L->na) { double_size(L); } //double size if we need more space
    L->nodes[L->n].id = id; //assign id
    L->nodes[L->n].v = (val_t) id2; //assign value
    L->n++; //add to length of list
    return &L->nodes[L->n]; //return pointer to new node
}


// Get and remove a node from the end of the list
int L_POP_ID(struct L * L, nid_t * id, nid_t * id2)
{
    if (!L){ return 0; } //can't pop a null pointer
    if (L->n < 1){ return 0; } //can't pop an empty list
    *id = L->nodes[L->n-1].id; //return this node's id
    *id2 = (nid_t) L->nodes[L->n-1].v; //return this node's value
    L->n--; //decrease list length counter
    return 1; //sucessfully popped
}


// sum the values of a list
val_t L_SUM_VALS(struct L * L)
{
    len_t i; //counter
    val_t s = NULL_VAL; //sum
    if (L && L->n > 0){
        s=0;
        for (i=0; i<L->n; i++) {
            s += L->nodes[i].v; //sum the values
        }
    } 
    return s; //return sum of list's values
}

