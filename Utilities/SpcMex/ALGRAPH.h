/*  ALGRAPH.h
 
    Adjacency-list representation of graph, and various utility functions
  
    To create a new graph:
       G = NEW_ALGRAPH(N, d);
       where N is the number of nodes
       d is directed (1) or undirected (0)

    Basic Functions:
    NEW_ALGRAPH  - create a new graph data structure
    ALG_ADD_NODE - add a node with a specified value
    ALG_ADD_EDGE - add an edge with a specified weight
    ALG_GET_NUM_NODES - get the number of nodes in the graph
    ALG_GET_NUM_EDGES - get the number of edges in the graph
    ALG_WEIGHT - find the weight of the edge between 2 points
    ALG_ADD_VAL_TO_EDGE - add value to a preexisting edge (or make new one if empty)
    ALG_EDGE_WEIGHTS - return vector with all the edge weights
    ALG_CLEAR - clear and free memory for an ALGRAPH data structure
    ALG_TO_STRING - print out all the edges as: FROM \t TO \t WEIGHT \n

    Utility functions:
    ALG_CYCLES - determine if there are cycles in the graph
    ALG_MST - find the minimal spanning tree of a graph

    TODO:
    ALG_KNN - find K Nearest Neighbors of each node
    ALG_KNNS - find K Nearest Neighbors of a single node
    ALG_KNNM - find the graph connecting Mutual K Nearest Neighbors
    ALG_KNNA - find KNN using approximation algorithm
    ALG_KNNS - find KNN of a single node using approx. algorithm
    ALG_KNNM - find the graph connecting Mutual KNN using approx. algorithm
    ALG_SUB - find connected subgraphs of a graph
    ALG_SUBAA - connect all nodes of each subgraph to each other
    ALG_AND - AND together 2 graphs with the same number of nodes
    ALG_OR - OR together 2 graphs with the same # of nodes ("superimpose")
    ALG_ADD - add the edge weights of 2 graphs with same number of nodes
    ALG_THRESH - remove all edges whose weights are less than some threshold
    ALG_DFS - perform depth-first search for vertex values
    ALG_SHORTEST_PATH - find shortest path between 2 points
    ALG_DELTR - make graph representing deluanay triangulation of points

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */


// What list datastructure are we using


// Rename the functions so we can easily switch out List implementation
// Binary tree
/*
#include "Trash.h" //List data structure
typedef TRASH_TREE List; //define what kind of list we're using
typedef trval_t edge_t; //type of data as edge weights, nid_t is still nid_t
typedef trid_t nid_t; //type of data as edge weights, nid_t is still nid_t
typedef trid_t len_t; //type of data as edge weights, nid_t is still nid_t
#define NEW_LIST(n) TRASH_NEW(n)
#define EMPTY_LIST(L) TRASH_EMPTY(L)
#define CLEAR_LIST(L) TRASH_CLEAR(L)
#define LIST_LENGTH(L) TRASH_LENGTH(L)
#define LIST_ADD(L, i, v) TRASH_ADD(L, i, v)
#define LIST_ADD_VAL(L, i, v) TRASH_ADD_VAL(L, i, v)
#define LIST_GET_VAL(L, i) TRASH_GET_VAL(L, i)
#define LIST_GET_ID(L, i) TRASH_GET_ID(L, i)
#define LIST_SET_VAL(L, i, v) TRASH_SET_VAL(L, i, v)
#define LIST_SET_ID(L, i) TRASH_SET_ID(L, i)
#define LIST_FIND(L, ind) TRASH_FIND(L, ind)
#define LIST_POP(L, ip, vp) TRASH_POP(L, ip, vp)
#define LIST_PUSH(L, ip, vp) TRASH_PUSH(L, ip, vp)
#define LIST_POP_ID(L, ip, vp) TRASH_POP_ID(L, ip, vp)
#define LIST_PUSH_ID(L, ip, vp) TRASH_PUSH_ID(L, ip, vp)
#define LIST_SUM_VALS(L) TRASH_SUM_VALS(L)
#define LIST_TO_STRING(L) TRASH_TO_STRING(L)
*/

// Normal list
#include "L.h" //List data structure
typedef struct L List; //define what kind of list we're using
typedef val_t edge_t; //type of data as edge weights, nid_t is still nid_t
#define NEW_LIST(n) NEW_L(n)
#define LIST_GET_VAL(L, i) L_GET_VAL(L, i)
#define LIST_GET_ID(L, i) L_GET_ID(L, i)
#define LIST_SET_VAL(L, i, v) L_SET_VAL(L, i, v)
#define LIST_SET_ID(L, i) L_SET_ID(L, i)
#define LIST_FIND(L, ind) L_FIND(L, ind)
#define LIST_ADD(L, i, v) L_ADD(L, i, v)
#define LIST_ADD_VAL(L, i, v) L_ADD_VAL_BY_ID(L, i, v)
#define LIST_POP(L, ip, vp) L_POP(L, ip, vp)
#define LIST_LENGTH(L) L_LENGTH(L)
#define LIST_TO_STRING(L) L_TO_STRING(L)
#define EMPTY_LIST(L) L_EMPTY(L)
#define CLEAR_LIST(L) L_CLEAR(L)
#define LIST_PUSH(L, ip, vp) L_ADD(L, ip, vp)
#define LIST_PUSH_ID(L, ip, vp) L_PUSH_ID(L, ip, vp)
#define LIST_POP_ID(L, ip, vp) L_POP_ID(L, ip, vp)
#define LIST_SUM_VALS(L) L_SUM_VALS(L)

// Defaults
#define DEF_EDGELIST_LEN 14
#define NO_EDGE_VAL NULL_VAL //value of an edge which is nonexistant
#define TRUE_EDGE_VAL 1.0 //value of a true edge

struct ALGRAPH
{
    List ** nodes; //pointer to an array of linked lists containing edges
    nid_t n;  //number of nodes
    nid_t ne;  //number of edges
    int dir; //directed (1) or undirected (0)
};

struct ALG_EDGE
{
    nid_t a; //startpoint vertex id
    nid_t b; //endpoing vertex id
    edge_t w; //weight of edge
};

// Basic functions
struct ALGRAPH * NEW_ALGRAPH(nid_t, int);
nid_t ALG_GET_NUM_NODES(struct ALGRAPH *);
nid_t ALG_GET_NUM_EDGES(struct ALGRAPH *);
void ALG_ADD_EDGE(struct ALGRAPH *, nid_t, nid_t, edge_t);
void ALG_OR_EDGE(struct ALGRAPH *, nid_t, nid_t, edge_t);
void ALG_ADD_VAL_TO_EDGE(struct ALGRAPH *, nid_t, nid_t, edge_t);
edge_t ALG_WEIGHT(struct ALGRAPH *, nid_t, nid_t);
int ALG_CONNECTED(struct ALGRAPH *, nid_t, nid_t);
edge_t * ALG_EDGE_WEIGHTS(struct ALGRAPH *);
edge_t ALG_AVG_EDGE_WEIGHT(struct ALGRAPH *);
edge_t ALG_AVG_SQ_EDGE_WEIGHT(struct ALGRAPH *);
void ALG_EMPTY(struct ALGRAPH *);
void ALG_CLEAR(struct ALGRAPH *);
void ALG_TO_STRING(struct ALGRAPH *);

// Cycles
int cycles_undir(struct ALGRAPH *, nid_t, List *, nid_t *);
int cycles_dir(struct ALGRAPH *, nid_t, List *, nid_t *);
int cycles_unalloc(struct ALGRAPH *, nid_t, List *, nid_t *);
int ALG_CYCLES(struct ALGRAPH *);

// Min spanning tree
int delete_mra_edge(struct ALGRAPH *, nid_t);
edge_t es_get_edge_val(struct ALG_EDGE *);
void es_set_edges_equal(struct ALG_EDGE *, struct ALG_EDGE *);
struct ALG_EDGE * es_get_es(struct ALGRAPH *) ;
struct ALGRAPH * kruskals(struct ALGRAPH *);
struct ALGRAPH * ALG_MST(struct ALGRAPH *);
void ALG_MST_pa(struct ALGRAPH *, struct ALGRAPH *);

// Subgraphs
void THRESH_SUBGRAPHS_SIZES(struct ALGRAPH *, nid_t *, edge_t, int, List *, nid_t *, nid_t *, nid_t *);
void THRESH_SUBGRAPHS(struct ALGRAPH *, edge_t, nid_t *, List *);
int label_subgraph(struct ALGRAPH *, List *, List *, nid_t *, nid_t, nid_t);
nid_t subgraphs_n(struct ALGRAPH *, List *, nid_t *, nid_t *);
nid_t subaa_additive(struct ALGRAPH *, struct ALGRAPH *, List *, List *, nid_t *, nid_t *);
nid_t ALG_SUB(struct ALGRAPH *, nid_t *);
struct ALGRAPH * ALG_SUBAA(struct ALGRAPH *);
struct ALGRAPH * NAIVE_EMST(edge_t *, nid_t, int);

// KNN
void brute_NN(edge_t *, struct ALG_EDGE *, nid_t, int, nid_t, nid_t);
void ALG_KNNS(edge_t *, nid_t *, edge_t *, nid_t, int, nid_t, nid_t);
void ALG_KNN(edge_t *, nid_t *, edge_t *, nid_t, int, nid_t);
struct ALGRAPH * ALG_KNNM(edge_t *, nid_t, int, nid_t);

// AND / OR / ADD / THRESH
struct ALGRAPH * ALG_AND(struct ALGRAPH *, struct ALGRAPH *);
struct ALGRAPH * ALG_OR(struct ALGRAPH *, struct ALGRAPH *);
void ALG_ORIP(struct ALGRAPH *, struct ALGRAPH *);
void ALG_ADD(struct ALGRAPH *, struct ALGRAPH *);
void ALG_THRESH(struct ALGRAPH *, struct ALGRAPH *,  edge_t);




