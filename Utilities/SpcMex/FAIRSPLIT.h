/*  FAIRSPLIT.h

    Fair-split decomposition and application functions

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */

// Define what graph data structure to use
#ifndef DEF_EDGELIST_LEN
    #include "ALGRAPH.h"
    typedef struct ALGRAPH Graph; //define the graph datastruct to use
#endif

typedef nid_t fsp_t; //type for "pointers"
typedef edge_t val_t; //type for data values
#define DEF_LIST_LEN 2


// Struct for list/queue
typedef struct
{
    fsp_t * list; //pointer to list of pairs
    fsp_t n; //number of pairs
    fsp_t na; //number of spots allocated
    fsp_t qb; //beginning index of the queue (always 0 for lists)
    fsp_t qe; //end index of the queue (will add to spot qe)
} FS_LIST;

// Struct for node in FS algorithm's doubly-linked lists
typedef struct 
{
    fsp_t data; //pointer to data for this point
    char half;   //what half of list this point is being divided into
    fsp_t adj_n;  //pointer to this point in next dimension's list
    fsp_t adj_p;  //pointer to this point in previous dimension's list
    fsp_t prev; //pointer to previous point in this dimension's list
    fsp_t next; //pointer to next point in this dimension's list
} FS_LNODE;


// Struct for node in FS tree
typedef struct 
{
    fsp_t p; //pointer to parent node
    fsp_t c1; //pointer to first child (or if leaf, pointer to this point's data)
    fsp_t c2; //pointer to second child (or if leaf, negative value)
    val_t * Ra; //pointer to array with lower bound for each dimension
    val_t * Rb; //pointer to array with upper bound for each dimension
} FS_TNODE;


// Struct for FS tree
typedef struct
{
    FS_TNODE * tree; //pointer to the tree array
    fsp_t root; //ind of the root of the tree in T (also = num elements in tree)
    fsp_t n; //number of nodes in the tree
    int d; //number of dimensions
    val_t * V; //pointer to values
} FS_TREE;


// Functions
FS_TREE * FS_DECOMP(val_t *, fsp_t, int);
FS_TREE * FAIR_SPLIT(val_t *, fsp_t, int);
void CLEAR_FS_TREE(FS_TREE *);
void CLEAR_FS_LIST(FS_LIST *);
void FS_TREE_TO_STRING(FS_TREE *);
void FS_LIST_POINTS_TO_STRING(FS_TREE *, FS_LIST *);
void to_string_tree_as_list(FS_TREE *);

// Spanning graph
void FS_AEMSG(FS_TREE *, Graph *);

// AKNN
//void FS_AKNN(FS_TREE *, fsp_t, fsp_t *, val_t *, int);
void FS_AMKNN(FS_TREE *, Graph *, fsp_t, int);

/*
// List functions
FS_LIST * new_wsp_list(fsp_t);
void fs_list_add(FS_LIST *, fsp_t, fsp_t);
fsp_t fs_list_numpairs(FS_LIST *);
void fs_list_get(FS_LIST *, fsp_t, fsp_t *, fsp_t *);
void fs_list_dequeue(FS_LIST *, fsp_t *, fsp_t *);

// single element list functions
void fs_list_add_one(FS_LIST *, fsp_t);
fsp_t fs_list_length(FS_LIST *);
void fs_list_get_one(FS_LIST *, fsp_t, fsp_t *);
void fs_list_dequeue_one(FS_LIST *, fsp_t *);
void FS_LIST_TO_STRING_SINGLES(FS_LIST *);
*/

