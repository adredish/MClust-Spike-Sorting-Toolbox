/*  ALGRAPH.c
 
    Adjacency-list representation of graph, and various utility functions
    NOTE: in this case it might not be best to use adjacency binary heap
    (instead of list) because the list len is max of k (so for us, 11)
  
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
    ALG_KNNS - find K Nearest Neighbors of a single node
    ALG_KNN - find K Nearest Neighbors of each node
    ALG_KNNM - find the graph connecting Mutual K Nearest Neighbors
    ALG_SUB - find connected subgraphs of a graph
    ALG_SUBAA - connect all nodes of each subgraph to each other
    ALG_AND - AND together 2 graphs with the same number of nodes
    ALG_OR - OR together 2 graphs with the same # of nodes ("superimpose")
    ALG_ADD - add the edge weights of 2 graphs with same number of nodes
    ALG_THRESH - remove all edges whose weights are less than some threshold

    TODO:
    ALG_DELTR - make graph representing deluanay triangulation of points
        NOTE: use this instead of KNNM for finding neighbors in low dimensions
    ALG_KNNA - find KNN using approximation algorithm
    ALG_KNNAS - find KNN of a single node using approx. algorithm
    ALG_KNNAM - find the graph connecting Mutual KNN using approx. algorithm
    ALG_DFS - perform depth-first search for vertex values
    ALG_SHORTEST_PATH - find shortest path between 2 points

    Notes:
    Have to store edges in both directions or the subgraphs functions won't work

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ALGRAPH.h"
#include "MergeSort.h"

//Allocate and initialize a new empty graph
struct ALGRAPH * NEW_ALGRAPH(nid_t N, int d)
{
    nid_t i; //counter
    struct ALGRAPH * G; //pointer to graph header
    List ** nodes; //pointer to node array
    G = (struct ALGRAPH *) malloc(sizeof(struct ALGRAPH)); //allocate graph header
    nodes = malloc(N*sizeof(List *)); //allocate array of pointers to list headers
    for (i=0; i<N; i++) {
        nodes[i] = NULL; //initially empty (null pointers)
    }
    G->nodes = nodes; //assign that array to the graph head
    G->n = N; //number of nodes
    G->ne = 0; //no edges yet
    G->dir = d; //directed (1) or undirected (0) graph
    return G; //return pointer to the graph head
}


// Get the number of nodes
nid_t ALG_GET_NUM_NODES(struct ALGRAPH * G)
{
    return G->n;
}


// Get the number of edges
nid_t ALG_GET_NUM_EDGES(struct ALGRAPH * G)
{
    return G->ne;
}


// Add an edge between nodes i and j of weight w
void ALG_ADD_EDGE(struct ALGRAPH * G, nid_t i, nid_t j, edge_t w)
{
    if (!G->nodes[i]) { //if the list is empty
        G->nodes[i] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
    }
    LIST_ADD(G->nodes[i], j, w); //add to the list
    G->ne++; //add to number of edges
    if (!G->dir) { //if the graph is undirected, have to add edge in opposite dir too (for cycles search)
        if (!G->nodes[j]) { //if j's list is empty
            G->nodes[j] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
        }
        LIST_ADD(G->nodes[j], i, w); //add to the list
    }
}


// Add an edge between nodes i and j of weight w
void ALG_OR_EDGE(struct ALGRAPH * G, nid_t i, nid_t j, edge_t w)
{
    if (!G->nodes[i]) { //if the list is empty
        G->nodes[i] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
        LIST_ADD(G->nodes[i], j, w); //add edge to list
        G->ne++; //add to number of edges
    } else {
        if ( LIST_FIND(G->nodes[i], j) == NO_EDGE_VAL ) { //if the edge doesn't already exist
            LIST_ADD(G->nodes[i], j, w); //add to the list
            G->ne++; //add to number of edges
        }
    }
    if (!G->dir) { //if undirected graph, have to add edge in other dir
        if (!G->nodes[j]) { //if the list is empty
            G->nodes[j] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
            LIST_ADD(G->nodes[j], i, w); //add edge to list
        } else {
            if ( LIST_FIND(G->nodes[j], i) == NO_EDGE_VAL ) { //NOTE: this could probably be optimized; if j was not found in i's list then we don't have to check j's list for i, since it will not be there
                LIST_ADD(G->nodes[j], i, w); //add to the list
            }
        }
    }
}


// Add a value to the weight between nodes i and j of weight w
void ALG_ADD_VAL_TO_EDGE(struct ALGRAPH * G, nid_t i, nid_t j, edge_t w)
{
    if (!G->nodes[i]) { //if the list is empty
        G->nodes[i] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
        LIST_ADD(G->nodes[i], j, w); //add to the list
        G->ne++; //add to number of edges
    } else {
        if ( LIST_ADD_VAL(G->nodes[i], j, w) ) { //add value to list (or add new element)
            G->ne++;  //but if it had to add an edge, add 1 to num edges
        } 
    }
    if (!G->dir) { //if undirected graph, have to add edge in other dir
        if (!G->nodes[j]) { //if the list is empty
            G->nodes[j] = NEW_LIST(DEF_EDGELIST_LEN); //allocate new list
            LIST_ADD(G->nodes[j], i, w); //add to the list
        } else {
            LIST_ADD_VAL(G->nodes[j], i, w); //add value to list (or add new element)
        }
    }
}


// Get the edge weight between nodes i and j
edge_t ALG_WEIGHT(struct ALGRAPH * G, nid_t i, nid_t j)
{
    if (!G->nodes[i]) { //if the list is empty
        return NO_EDGE_VAL; //return no edge weight
    } else {
        return LIST_FIND(G->nodes[i], j); //return the edge weight (value of id=j in i's list)
    }
}


// Find if two nodes are connected
int ALG_CONNECTED(struct ALGRAPH * G, nid_t i, nid_t j)
{
    if ( ALG_WEIGHT(G, i, j) == NO_EDGE_VAL ) { //if there is no edge or it has null weight,
        return 0;
    } else {
        return 1;
    }
}


// Get an array containing all the edge weights
edge_t * ALG_EDGE_WEIGHTS(struct ALGRAPH * G)
{
    nid_t i, k=0; //counter
    len_t j;
    edge_t * W; //pointer to array of weights
    W = (edge_t *) malloc(G->ne*sizeof(edge_t)); //allocate the array
    for (i=0; i<G->n; i++) { //for each node's list,
        for (j=0; j<LIST_LENGTH(G->nodes[i]); j++) { //for each element of that list
            if ( i < LIST_GET_ID(G->nodes[i], j) ) { //only keep the edge originating from the lesser id
                W[k] = LIST_GET_VAL(G->nodes[i], j); //save that edge weight
                k++; //add to counter
            }
        }
    }
    return W; //return pointer to the array of edge weights
}


// Get the sum of all the edge weights, recursively
edge_t alg_sum_recurse(struct ALGRAPH * G, nid_t a, nid_t b)
{
    edge_t sv1, sv2;
    nid_t m;
    if (b-a > 0) { //several lists, recursively split to avoid mag. error
        m = (a+b)/2; //midpoint
        sv1 = alg_sum_recurse(G, a, m);
        sv2 = alg_sum_recurse(G, m+1, b);
        return sv1+sv2; //sum of both halves
    } else { //one list
        return LIST_SUM_VALS(G->nodes[b]);
    }
}


// Get the sum of all the edge weights
edge_t ALG_SUM_EDGE_WEIGHTS(struct ALGRAPH * G)
{
    if (G->dir) {
        return alg_sum_recurse(G, 0, G->n-1);
    } else {
        return 0.5*alg_sum_recurse(G, 0, G->n-1); //halve b/c duplicate edges in dir graph
    }
}


// Get the average of all the edge weights
edge_t ALG_AVG_EDGE_WEIGHT(struct ALGRAPH * G)
{
    return ALG_SUM_EDGE_WEIGHTS(G)/((edge_t) G->ne);
}


// Get the sum of all the squared edge weights, recursively
edge_t alg_avg_sq_recurse(struct ALGRAPH * G, nid_t a, nid_t b, edge_t mult)
{
    edge_t sv1, sv2;
    int m;
    nid_t i;
    if (b-a > 0) { //several lists, recursively split to avoid mag. error
        m = (a+b)/2; //midpoint
        sv1 = alg_avg_sq_recurse(G, a, m, mult);
        sv2 = alg_avg_sq_recurse(G, m+1, b, mult);
        return sv1+sv2; //sum of both halves
    } else { //one list, find the sum of squares
        sv1 = 0.0;
        for (i=0; i<LIST_LENGTH(G->nodes[a]); i++) {
            sv2 = LIST_GET_VAL(G->nodes[a], i);
            sv1 += sv2*sv2*mult; //squared average
        }
        return sv1; //return the sum of squares (averaged)
    }
}


// Get the sum of all the squared edge weights
edge_t ALG_AVG_SQ_EDGE_WEIGHT(struct ALGRAPH * G)
{
    if (G->dir) {
        return alg_avg_sq_recurse(G, 0, G->n-1, 1.0/((edge_t) G->ne));
    } else {
        return alg_avg_sq_recurse(G, 0, G->n-1, 0.5/((edge_t) G->ne));
    }
}


// Empties a graph datastructure without deallocating it from memory
void ALG_EMPTY(struct ALGRAPH * G)
{
    nid_t i; //counter
    for (i=0; i<G->n; i++) { //for each node's list,
        if (G->nodes[i]) { //if it's not a null pointer
            EMPTY_LIST(G->nodes[i]); //empty but don't deallocate
        }
    }
    G->ne = 0; //no edges (but same number of nodes)
}


// Clears a graph data structure from memory
void ALG_CLEAR(struct ALGRAPH * G)
{
    nid_t i; //counter
    for (i=0; i<G->n; i++) { //for each node's list,
        if ( G->nodes[i] ) { //if that list isn't empty, (not null pointer)
            CLEAR_LIST(G->nodes[i]);
        }
    }
    free(G->nodes); //free the list of pointers to nodes
    free(G); //free the graph header
}


// Prints the edges of a graph
void ALG_TO_STRING(struct ALGRAPH * G)
{
    nid_t i; //counters
    len_t j;
    if ( G->ne < 1 ) { //if no edges
        printf("Empty Graph\n"); 
        return;
    }
    printf("FROM\tTO\tWEIGHT\n");
    for (i=0; i<G->n; i++) { //for each node's list,
        for (j=0; j<LIST_LENGTH(G->nodes[i]); j++) { //for each element of that list
            if ( i < LIST_GET_ID(G->nodes[i], j) ) { //only print the edge originating from the lesser id
                printf("%d\t%d\t%f\n", i, LIST_GET_ID(G->nodes[i],j), LIST_GET_VAL(G->nodes[i],j));
            }
        }
    }
}





/*****************************************************************************
                               CYCLES AND DFS
*****************************************************************************/

// Determine if there are cycles in a graph, using depth-first search, 
// starting at node with ind vi, S is stack (preferably empty but doesn't
// have to be), and L is list of already-visited nodes
// NOTE: this works only with undirected graphs
int cycles_undir(struct ALGRAPH * G, nid_t vi, List * S, nid_t * L)
{
    len_t i; //counter
    nid_t pid, cid, pval; //popped id and child's id, and parents id
    LIST_PUSH_ID(S, vi, 0); //push first node to stack
    while ( LIST_LENGTH(S) ) { //while the stack is not empty
        LIST_POP_ID(S, &pid, &pval);  //pop top element on stack
        if ( L[pid] ) { return 1; } //if we've already visited this node, a cycle was found
        L[pid] = 1; //mark it as visited
        for (i=0; i<LIST_LENGTH(G->nodes[pid]); i++) { //for each child of this node
            cid = LIST_GET_ID(G->nodes[pid], i) ; //get this child's id
            if ( pval != cid ) { //if this node did not refer it
                LIST_PUSH_ID(S, cid, pid);  //add child to stack w/ parent id as value
            }
        }
    }
    return 0;  //we got to the end, no cycles apparently!
}


// Detects cycles in directed graphs
int cycles_dir(struct ALGRAPH * G, nid_t vi, List * S, nid_t * L)
{
    len_t i; //counter
    nid_t pid, cid, pval; //popped id and child's id, and parents id
    LIST_PUSH_ID(S, vi, 0); //push first node to stack
    while ( LIST_LENGTH(S) ) { //while the stack is not empty
        LIST_POP_ID(S, &pid, &pval);  //pop top element on stack
        if ( L[pid] ) { return 1; } //if we've already visited this node, a cycle was found
        L[pid] = 1; //mark it as visited
        for (i=0; i<LIST_LENGTH(G->nodes[pid]); i++) { //for each child of this node
            cid = LIST_GET_ID(G->nodes[pid], i) ; //get this child's id
            LIST_PUSH_ID(S, cid, 0);  //add child to stack
        }
    }
    return 0;  //we got to the end, no cycles apparently!
}

// Searches for cycles in either dir or undir, without allocating new stack + list
// however it re-initializes stack and list to empty and all false
int cycles_unalloc(struct ALGRAPH * G, nid_t vi, List * S, nid_t * L)
{
    nid_t i; //counter
    nid_t pid, cid; //popped id and child's id
    edge_t pval; //popped value
    for (i=0; i<G->n; i++) { L[i] = 0; } //initialize the labels to false
    EMPTY_LIST(S); //initialize stack to empty
    if (G->dir) { //if directed graph,
        return cycles_dir(G, vi, S, L); //do a directed cycles search from that node
    } else {
        return cycles_undir(G, vi, S, L); //do an undirected cycles search
    }
}


// Works on either directed or undirected, for user calling
int ALG_CYCLES(struct ALGRAPH * G)
{
    nid_t i,j; //counter
    int c=0; //cycles flag
    List * S; //stack
    nid_t * L, * V;
    S = NEW_LIST(10); //new stack
    L = malloc(G->n*sizeof(nid_t)); //label as to whether we've visited each node
    V = malloc(G->n*sizeof(nid_t)); //and one for the cycles func
    for (i=0; i<G->n; i++) { V[i]=0; L[i]=0; } //initialize the labels to false
    if ( G->dir ) { //if this is a directed graph,
        for (i=0; i<G->n; i++) { //for each node
            if ( !V[i] ) { //if we haven't cycle searched from this node yet,
                c = cycles_unalloc(G, i, S, L); //then we need to have a fresh L each time
                if (c) { break; } //found a cycle!!!
                for (j=0; j<G->n; j++) { V[j] += L[j]; } //add to who has been visited (since we only need to do a new DFS from nodes who haven't been visited yet)
            }
        }
    } else { //undirected, can use same L array
        for (i=0; i<G->n; i++) { //for each node
            if ( !L[i] ) { //if this hasn't been hit by a dfs,
                c = cycles_undir(G, i, S, L); //then we don't need to have a fresh L each time
                if (c) { break; } //found a cycle!!!
            }
        }
    }
    CLEAR_LIST(S); //free the stack
    free(L); //and the lists
    free(V);
    return c;
}




/*****************************************************************************
                           MINIMUM SPANNING TREE
*****************************************************************************/

// Delete the most recently added edge from a given node
int delete_mra_edge(struct ALGRAPH * G, nid_t ind)
{
    edge_t v;  nid_t id;  nid_t id2; //for pop return
    if (LIST_LENGTH(G->nodes[ind]) == 0) { //if the list is empty
        return 0; //no edge was deleted b/c there was no edge to delete
    } else { 
        LIST_POP(G->nodes[ind], &id, &v); //pop the end of that list off
        G->ne--; //subtract from edge counter
        if (!G->dir) { //if this is an undirected graph,
            LIST_POP(G->nodes[id], &id2, &v); //pop the end of other node's list too
        }
        return 1; //success!
    }
}


// Get edge value from an edge struct
edge_t es_get_edge_val(struct ALG_EDGE * E)
{
    return E->w;
}


// Set an edge struct equal to another edge struct
void es_set_edges_equal(struct ALG_EDGE * E1, struct ALG_EDGE * E2)
{
    E1->a = E2->a; //set E1=E2
    E1->b = E2->b;
    E1->w = E2->w;
}


// Get array of Edge Struct
struct ALG_EDGE * es_get_es(struct ALGRAPH * G) 
{
    nid_t i, tid, k; //counter
    len_t j;
    struct ALG_EDGE * W; //pointer to array of weights
    k=0;
    W = malloc(G->ne*sizeof(struct ALG_EDGE)); //allocate the array
    for (i=0; i<G->n; i++) { //for each node's list,
        for (j=0; j<LIST_LENGTH(G->nodes[i]); j++) { //for each element of that list
            tid = LIST_GET_ID(G->nodes[i], j); //save that edge endpoint's id
            if (i < tid) {
                W[k].a = (nid_t) i;
                W[k].b = tid; //save that edge endpoint's id
                W[k].w = LIST_GET_VAL(G->nodes[i], j); //save that edge weight
                k++; //add to counter
            }
        }
    }
    return W; //return pointer to the array of edge structs
}

// Find the minimum spanning tree of a graph, using Kruskal's algorithm
void ALG_MST_pa(struct ALGRAPH * G, struct ALGRAPH * M)
{
    // Declare variables
    struct ALG_EDGE * E;
    nid_t ec, ei, i; //num. edges in mst so far, and the edge we're on
    nid_t * L;
    List * S; //stack for cycles
    nid_t * CL;

    // Make a sorted list of the edges
    E = es_get_es(G); //get the array of edge structs for this graph
    MergeSort(E, G->ne); //sort the structs (2* b/c it's an undirected graph)

    // Make the minimum spanning tree
    ec=0; ei=0; //num. edges in mst so far, and the edge we're on
    L = malloc(G->n*sizeof(nid_t)); //label vector for whether we've visited a node already
    for (i=0; i<G->n; i++) { L[i] = 0; } //initialize all to false
    S = NEW_LIST(10); //new stack
    CL = malloc(G->n*sizeof(nid_t)); //label vector for whether we've visited a node already in DFS
    while (ec < G->n-1) { //while there aren't enough edges to connect all vertexes,
        if (E[ei].a < E[ei].b) { //only add smaller id->larger id (since one add adds both ways)
            ALG_ADD_EDGE(M, E[ei].a, E[ei].b, E[ei].w); //add next edge from sorted list
            if ( L[E[ei].a] && L[E[ei].b] ) { //but if both those vertexes already have edges coming from them, there's a risk of causing a cycle
                if ( cycles_unalloc(M, E[ei].a, S, CL) ) { //so check if adding that edge causes there to be cycles, and if it does,
                    delete_mra_edge(M, E[ei].a); //delete that edge you just added
                    ec--; //update counter
                }
            } else { //at least one previously had no edges coming out of it
                L[E[ei].a] = 1; L[E[ei].b] = 1; //mark these nodes as visited
            }
            ec++; //add to counter for how many edges we've added to MST
        }
        ei++; //add to counter for where we are in sorted list
    }
    CLEAR_LIST(S); //free the list from mem
    free(E);
    free(L);
    free(CL);
}



// Find the minimum spanning tree of a graph, using Kruskal's algorithm
struct ALGRAPH * kruskals(struct ALGRAPH * G)
{
    // Declare variables
    struct ALG_EDGE * E;
    nid_t ec, ei, i; 
    nid_t * L;
    struct ALGRAPH * M;
    List * S;
    nid_t * CL;

    // Make a sorted list of the edges
    E = es_get_es(G); //get the array of edge structs for this graph
    MergeSort(E, G->ne); //sort the structs (2* b/c it's an undirected graph)

    // Make the minimum spanning tree
    ec=0; ei=0; //num. edges in mst so far, and the edge we're on
    L = malloc(G->n*sizeof(nid_t)); //label vector for whether we've visited a node already
    for (i=0; i<G->n; i++) { L[i] = 0; } //initialize all to false
    M = NEW_ALGRAPH(G->n, 0); //allocate new graph for MST w/ same number nodes + undirected
    S = NEW_LIST(10); //new stack
    CL = malloc(G->n*sizeof(nid_t)); //label vector for whether we've visited a node already in DFS
    while (ec < G->n-1) { //while there aren't enough edges to connect all vertexes,
        if (E[ei].a < E[ei].b) { //only add smaller id->larger id (since one add adds both ways)
            ALG_ADD_EDGE(M, E[ei].a, E[ei].b, E[ei].w); //add next edge from sorted list
            if ( L[E[ei].a] && L[E[ei].b] ) { //but if both those vertexes already have edges coming from them, there's a risk of causing a cycle
                if ( cycles_unalloc(M, E[ei].a, S, CL) ) { //so check if adding that edge causes there to be cycles, and if it does,
                    delete_mra_edge(M, E[ei].a); //delete that edge you just added
                    ec--; //update counter
                }
            } else { //at least one previously had no edges coming out of it
                L[E[ei].a] = 1; L[E[ei].b] = 1; //mark these nodes as visited
            }
            ec++; //add to counter for how many edges we've added to MST
        }
        ei++; //add to counter for where we are in sorted list
    }
    CLEAR_LIST(S); //free the list from mem
    free(E);
    free(L);
    free(CL);
    return M; //return the minimum spanning tree
}


// Find the minimum spanning tree of a graph
struct ALGRAPH * ALG_MST(struct ALGRAPH * G)
{
    return kruskals(G);  //just use kruskal's for now, but would like to implement Chazelle's soft heap algorithm or the well-separated pairs decomp
}


/*  Find the minimum spanning tree given points
    NOTES: we can't find the (NxN-N)/2 upper triangular matrix of distances and
    then use the previous MST function, b/c that would be 1TB of space for 500k points...
    That is, it would be N^2 space AND N^2 time.  :(
    which means we can't use Kruskal's algorithm (for which we need to SORT said N^2 data)
    also can't use prims which is O(E) and here E=N^2
    
    The minimum spanning tree of a set of points is contained in the Delaunay triangulation
    graph, so we could just 1) find the Delaunay triangulation, and then 2) run our existing
    ALG_MST function on it to get the MST of the points.
    It looks like DWYER's divide and conquer delaunay triangluation algorithm is the fastest
    for higher dimensions
    See 
        Dwyer 1987 "A Faster Divide-and-Conquer Algorithm for Constructing Delaunay Triangulations"
        Leach 1992 "Improving Worst-Case Optimal Delaunay Triangulation Algorithms"
        Cignoni et al 1998 "DeWall: A fast divide and conquer Delaunay triangulation algorithm in Ed"

    Alternatively, we can use an approximation algorithm, since the purpose of finding
    the MST for SPC is not to find the exact mst, but rather to ensure that all points 
    are connected.  This can be done in O(n log n) time, supposedly
    See
        Smid 2005 (unpublished) "The well-separated pair decomp. and its applications", section 5.3
        Callahan and Kosaraju (Ch 32 in a book?) "Faster Algorithms for Some geometric graph problems in higher dimensions"
    Should probably use that last approx one w/ well-sep pair decomp.
*/
// Find a euclidian minimum spanning tree for a set of points, using O(n^2) naive algorithm
struct ALGRAPH * NAIVE_EMST(edge_t * V, nid_t n, int d)
{
    int di;
    edge_t ddiff, dsum;
    nid_t i, j;
    struct ALGRAPH * CEG; //the complete euclidian graph
    CEG = NEW_ALGRAPH(n, 0); //allocate new graph for MST w/ same number nodes + undirected

    // Connect edges all-to-all
    for (i=0; i<n; i++) { //connect all points to all other points
        for (j=i+1; j<n; j++) {
            dsum = 0;
            for (di=0; di<d; di++) { //compute distance between two points
                ddiff = *(V+d*i+di)-*(V+d*j+di);
                dsum += ddiff*ddiff;
            }
            dsum = sqrt(dsum); //euclidian dist between points
            ALG_ADD_EDGE(CEG, i, j, dsum); //add edge between points i and j
        }
    }

    // Find the MST of that graph
    return ALG_MST(CEG);
}





/*****************************************************************************
                                  SUBGRAPHS
*****************************************************************************/


/*  Finds the subgraphs of a graph with the edges thresholded at
    some value, and assigns subgraph ids.  Given the number of largest
    clusters to keep track of, also returns the sizes and ids of the
    largest clusters

    INPUTS
        Tb - the graph
        maxCor - list of maximum-correlation points
        t - the threshold (value above which the edge will be deleted)
        Nmc - number of largest clusters to keep track of
        S - a List datastructure which serves as the breadth-first search stack
    OUTUPTS
        C - subgraph ids
        MCS - sizes of the Nms-largest clusters
        MCI - ids of the Nms-largest clusters
*/
void THRESH_SUBGRAPHS_SIZES(struct ALGRAPH * Tb, nid_t * maxCor, edge_t t, int Nmc, List * S, nid_t * C, nid_t * MCS, nid_t * MCI)
{
    // Declare variables
    nid_t i, j, pid, pval, cid, tl=1; 
    edge_t cval;
    nid_t * clusSizes = calloc(tl,sizeof(nid_t)); //the size of each cluster id
    nid_t mid, mval;
    
    // Do the subgraphs
    for (i=0; i<Tb->n; i++) { C[i]=0; } //initialize to false
    for (i=0; i<Tb->n; i++) { //for each node
        if ( !C[i] ) { //if we haven't visited this node before in dfs,
            //label the subgraph starting at that point
            LIST_PUSH_ID(S, i, -1); //push this node to dfs stack
            while ( LIST_LENGTH(S) ) { //while the stack is not empty,
                LIST_POP_ID(S, &pid, &pval); //pop top element of stack
                if ( !C[pid] ) { //as long as we haven't already visited this node,
                    C[pid] = tl; //mark it with current label
                    for (j=0; j<LIST_LENGTH(Tb->nodes[pid]); j++) { //for each child of this ndoe,
                        cid = LIST_GET_ID(Tb->nodes[pid], j); //get this child's id
                        cval = LIST_GET_VAL(Tb->nodes[pid], j); //get this child's value
                        // old version (which worked)
                        /*
                        if ( pval != cid && cval > t ) { //if this node did not refer it and the edge value is less than the threshold
                            LIST_PUSH_ID(S, cid, pid); //add this child to dfs stack
                        }
                        */
                        
                        // new version (which always connects a point along it's max-cor edge)
                        if ( (pval != cid && cval > t ) || (maxCor[cid] == pid) ) { //if this node did not refer it and the edge value is less than the threshold  OR if I am the max-cor edge of that point,
                            LIST_PUSH_ID(S, cid, pid); //add this child to dfs stack
                        }

                    }
                }
            }
            tl++; //increment the subgraph label id
        }
    }

    // Now find the sizes of the largest clusters
    clusSizes = calloc(tl,sizeof(nid_t)); //the size of each cluster id
    for (i=0; i<Tb->n; i++) { //for each point
        clusSizes[C[i]]++; //keep track of how many points are in each cluster
    }
    for (i=0; i<Nmc; i++) { //for however many max clusters were requested
        mval=0; //reset to 0
        for (j=0; j<tl; j++) { //for each cluster,
            if (clusSizes[j] > mval) { //if this is largest so far,
                mval = clusSizes[j]; //store the size
                mid = j; //and the id
            }
        }
        MCS[i] = mval; //save the i-th largest cluster size
        MCI[i] = mid; //and its corresponding id
        clusSizes[mid] = 0; //ensure we don't pick it out again
    }
    free(clusSizes);
}


/*  Finds the subgraphs of a graph with the edges thresholded at
    some value, and assigns subgraph ids

    INPUTS
        Tb - the graph
        t - the threshold (value above which the edge will be deleted)
    OUTUPTS
        C - subgraph ids
*/
void THRESH_SUBGRAPHS(struct ALGRAPH * Tb, edge_t t, nid_t * C, List * S)
{
    nid_t i, j, pid, pval, cid, tl=1; 
    edge_t cval;
    for (i=0; i<Tb->n; i++) { C[i]=0; } //initialize to false
    for (i=0; i<Tb->n; i++) { //for each node
        if ( !C[i] ) { //if we haven't visited this node before in dfs,
            //label the subgraph starting at that point
            LIST_PUSH_ID(S, i, -1); //push this node to dfs stack
            while ( LIST_LENGTH(S) ) { //while the stack is not empty,
                LIST_POP_ID(S, &pid, &pval); //pop top element of stack
                if ( !C[pid] ) { //as long as we haven't already visited this node,
                    C[pid] = tl; //mark it with current label
                    for (j=0; j<LIST_LENGTH(Tb->nodes[pid]); j++) { //for each child of this ndoe,
                        cid = LIST_GET_ID(Tb->nodes[pid], j); //get this child's id
                        cval = LIST_GET_VAL(Tb->nodes[pid], j); //get this child's value
                        if ( pval != cid && cval > t ) { //if this node did not refer it and the edge value is less than the threshold
                            LIST_PUSH_ID(S, cid, pid); //add this child to dfs stack
                        }
                    }
                }
            }
            tl++; //increment the subgraph label id
        }
    }
}




// Sets all labels (in L) for nodes in this subgraph to tl, returns number of elements
// in this subgraph.  If EL is a pointer to a list (and not null), will return a list
// of the ids for all nodes in this subgraph
int label_subgraph(struct ALGRAPH * G, List * EL, List * S, nid_t * L, nid_t vi, nid_t tl)
{
    len_t i; //counter thru list
    nid_t nl=0; //number of nodes labeled
    nid_t pid, cid, pval; //popped id and child's id, and parent id
    LIST_PUSH_ID(S, vi, -1); //push first node to stack
    while ( LIST_LENGTH(S) ) { //while the stack is not empty
        LIST_POP_ID(S, &pid, &pval);  //pop top element on stack
        if ( !L[pid] ) { //as long as we haven't visited this node already,
            L[pid] = tl; //mark it as visited w/ requested id
            nl++; //add to count of how many we labeled
            if (EL) { LIST_PUSH_ID(EL, pid, pval); } //if we want list of ids in this subgraph, add to it
            for (i=0; i<LIST_LENGTH(G->nodes[pid]); i++) { //for each child of this node
                cid = LIST_GET_ID(G->nodes[pid], i); //get this child's id
                if ( pval != cid ) { //if this node did not refer it
                    LIST_PUSH_ID(S, cid, pid);  //add child to stack w/ parent id as value
                }
            }
        }
    }
    return nl;  //return how many nodes were in this subgraph
}


// Label subgraphs, find max subgraph size, and 
nid_t subgraphs_n(struct ALGRAPH * G, List * S, nid_t * L, nid_t * ms)
{
    nid_t i, si=1, ss; //counters for nodes, subgraphs, and max subgraph size
    *ms = 0; //default is size-1 max subgraph size
    EMPTY_LIST(S); //ensure stack is empty
    for (i=0; i<G->n; i++) { L[i]=0; } //initialize to false
    for (i=0; i<G->n; i++) {
        if ( !L[i] ) { //if we haven't visited this node before in dfs,
            ss = label_subgraph(G, NULL, S, L, i, si); //label the subgraph at that point with unique id
            si++; //increment id
            if (ss>*ms) { *ms=ss; } //if this is larger than previous largest, remember this one
        }
    }
    return si-1; //return number of subgraphs that were labelled

}


/* Connects all nodes in each subgraphs of an UNDIRECTED graph to all other nodes in that subgraph
    INPUTS:
        G - graph of interest
        AAG - add all-to-all connections between subgraphs in G to AAG
        Sl - list to store nodes in this subgraph
        S - stack for DFS
        L - array of cluster label ids for each node
        ms - will set this to the size of the largest cluster
    Returns the number of clusters
*/
nid_t subaa_additive(struct ALGRAPH * G, struct ALGRAPH * AAG, List * Sl, List * S, nid_t * L, nid_t * ms)
{
    nid_t i, j, k, ss, si=1; //counters for nodes, subgraphs, and max subgraph size
    nid_t pid; //popped id
    edge_t pval; //popped value
    EMPTY_LIST(S); //ensure stack is empty
    EMPTY_LIST(Sl); //ensure stack is empty
    for (i=0; i<G->n; i++) { L[i]=0; } //initialize to false
    for (i=0; i<G->n; i++) { //NOTE this loop cannot be paralellized, but some of the nested loops can be
        if ( !L[i] ) { //if we haven't visited this node before in dfs,
            ss = label_subgraph(G, Sl, S, L, i, si); //label the subgraph at that point with unique (nonzero) id - the subgraph search cannot be paralellized
            si++; //increment id
            if (ss>*ms) { *ms=ss; } //if this is larger than previous largest, remember this one
            for (j=0; j<ss; j++) { //connect all elements in this subgraph - NOTE: this loop can be paralellized!!!
                for (k=j+1; k<ss; k++) { //to all other elements
                    if (j!=k) { //but not to themselves
                        ALG_ADD_VAL_TO_EDGE(AAG, LIST_GET_ID(Sl, j), LIST_GET_ID(Sl, k), TRUE_EDGE_VAL);

                    }
                }
            }
            EMPTY_LIST(S); //initialize stack to empty (without deallocating)
            EMPTY_LIST(Sl); //initialize stack to empty (")
        }
    }
    return si-1; //return number of clusters labelled
}


/*  Finds all subgraphs in an UNDIRECTED graph
    G is graph, L is N-length vector which will store cluster ids for each node
    Returns the number of subgraphs
*/
nid_t ALG_SUB(struct ALGRAPH * G, nid_t * L)
{
    nid_t i, si, ms; //counters for nodes, subgraphs, and max subgraph size
    List * S; //stack for dfs
    S = NEW_LIST(10); //new stack
    si = subgraphs_n(G, S, L, &ms);
    CLEAR_LIST(S); //free the list from memory
    return si; //return number of subgraphs that were labelled
}


// Clean version for user calling
struct ALGRAPH * ALG_SUBAA(struct ALGRAPH * G)
{
    nid_t * L;
    nid_t ms; //max subgraph size
    List * S;
    List * Sl;
    struct ALGRAPH * AAG; 
    L = malloc(G->n*sizeof(nid_t)); //label vector
    S = NEW_LIST(10); //stack for DFS
    Sl = NEW_LIST(10); //list of ids in this subgraph
    AAG = NEW_ALGRAPH(G->n, 0); //graph of all-to-all connections
    subaa_additive(G, AAG, Sl, S, L, &ms);
    CLEAR_LIST(S); //free the list from mem
    CLEAR_LIST(Sl); //free the list from mem
    free(L);
    return AAG;
}





/*****************************************************************************
                              K NEAREST NEIGHBORS 
*****************************************************************************/

// Finds distances to nearest neighbors and sorts
// V should be pointer to first element of n-by-d array
// D should be a pre-allocated n-length array of edge structs
// Returns a sorted list of edge datastructs
void brute_NN(edge_t * V, struct ALG_EDGE * D, nid_t n, int d, nid_t k, nid_t vi)
{
    nid_t i; int j; //counters (one for nodes, other for dims)
    edge_t s, diff, m=0; //sum for dist calc, difference
    for (i=0; i<n; i++) { //for each other node
        s = 0.0; //sum for euclidian distance calculation
        for (j=0; j<d; j++) {
            diff = *(V+vi*d+j)-*(V+i*d+j); //difference between 2 elements of 2d array
            s += diff*diff; //pow 2
        }
        D[i].a = vi;
        D[i].b = i;
        D[i].w = sqrt(s);
        if (D[i].w > m) { m = D[i].w; } //keep track of max edge
    } //NOTE: this loop calcs distance to self, but that's faster than to check each loop if i==self
    D[vi].a = vi; D[vi].b=vi; D[vi].w = m; //make self the furthest neighbor 
    MergeSort(D, n); //sort the distances
}


// Find the k-nearest neighbors of node vi in a n-by-d matrix V of points
// n is number of nodes, d is number of dimensions
// NN is nearest neighbor storage array (length K)
// D is distances to those nearest neighbors (length K
void ALG_KNNS(edge_t * V, nid_t * NN, edge_t * D, nid_t n, int d, nid_t k, nid_t vi)
{
    nid_t i;
    struct ALG_EDGE * DS; //array of distances
    DS = malloc(n*sizeof(struct ALG_EDGE)); //allocate the array
    brute_NN(V, DS, n, d, k, vi); //find distances and sort for this point
    for (i=0; i<k; i++) { //assign k nearest values
        NN[i] = DS[i].b; //i-th nearest neighbor's id
        D[i] = DS[i].w; //distance to i-th nearest neighbor
    }
}


// Find the k-nearest neighbors of all points in a n-by-d matrix V of points
// n is number of nodes, d is number of dimensions
// NN is nearest neighbor storage array (n-by-K)
// D is distances to those nearest neighbors (n-by-K)
void ALG_KNN(edge_t * V, nid_t * NN, edge_t * D, nid_t n, int d, nid_t k)
{
    nid_t i, j;
    struct ALG_EDGE * DS; //array of distances
    DS = malloc(n*sizeof(struct ALG_EDGE)); //allocate the array
    for (i=0; i<n; i++) { //for each node,
        brute_NN(V, DS, n, d, k, i); //find distances and sort for this point
        for (j=0; j<k; j++) { //assign k nearest values
            *(NN+i*k+j) = DS[j].b; //i-th nearest neighbor's id
            *(D+i*k+j) = DS[j].w; //distance to i-th nearest neighbor
        }
    }
    free(DS); //fre the array
}


// Finds the graph connecting mutual K-nearest neighbors
struct ALGRAPH * ALG_KNNM(edge_t * V, nid_t n, int d, nid_t k)
{
    // Declare variables
    nid_t * NN; //ids of nearest neighbors
    edge_t * D; //distances to nearest neighbors
    nid_t tnn; //store this nearest neighbor id
    nid_t i, j, m; //counters
    struct ALGRAPH * G;
    
    // Find nearest neighbors
    NN = malloc(n*k*sizeof(nid_t)); //allocate the arrays
    D = malloc(n*k*sizeof(edge_t));
    ALG_KNN(V, NN, D, n, d, k); //find nearest neighbors + distances

    // Make the graph
    G = NEW_ALGRAPH(n, 0); //allocate new graph
    for (i=0; i<n; i++) { //for each node
        for (j=0; j<k; j++) { //for each nearest neighbor of that node
            tnn = *(NN+i*k+j);
            if ( i < tnn ) { //only do lesser id -> larger id 
                for (m=0; m<k; m++) { //check if they're mutual neighbors
                    if ( *(NN+tnn*k+m) == i ) { //if they are in each other's lists,
                        ALG_ADD_EDGE(G, i, tnn, *(D+i*k+j)); //add edge with dist as weight
                        break;
                    }
                }
            }
        }
    }
    free(D); //free the distances array
    free(NN); //free the ids array
    return G; //return the graph of mutual nearest neighbors
}








/*****************************************************************************
                             AND / OR / ADD / THRESH
*****************************************************************************/

// AND two graphs together (only keep edges which exist in both graphs
struct ALGRAPH * ALG_AND(struct ALGRAPH * G1, struct ALGRAPH * G2)
{
    nid_t i; len_t j; //counters
    nid_t tid; edge_t tw; //this id and edge weight
    struct ALGRAPH * Go;
    Go = NEW_ALGRAPH(G1->n, G1->dir); //new graph of same type as input
    for (i=0; i<G1->n; i++) { //for each node
        for (j=0; j<LIST_LENGTH(G1->nodes[i]); j++) { //for each edge coming from that node
            tid = LIST_GET_ID(G1->nodes[i], j); //id of this edge 
            if ( G1->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                tw = ALG_WEIGHT(G2, i, tid); //weight of this edge in g2
                if ( tw != NO_EDGE_VAL ) { //if they're also connected in g2
                    ALG_ADD_EDGE(Go, i, tid, LIST_GET_VAL(G1->nodes[i], j)); //add the edge with the weight from g1
                }
            }
        }
    }
    return Go; //return the output graph
}


// OR 2 graphs together ("superimpose" them)
struct ALGRAPH * ALG_OR(struct ALGRAPH * G1, struct ALGRAPH * G2)
{
    nid_t i; len_t j; //counters
    nid_t tid; //this id
    struct ALGRAPH * Go;
    Go = NEW_ALGRAPH(G1->n, G1->dir); //new graph of same type as input
    for (i=0; i<G1->n; i++) { //for each node; add all G1 edges
        for (j=0; j<LIST_LENGTH(G1->nodes[i]); j++) { //for each edge coming from that node
            tid = LIST_GET_ID(G1->nodes[i], j); //id of this edge's endpoint
            if ( G1->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                ALG_ADD_EDGE(Go, i, tid, LIST_GET_VAL(G1->nodes[i], j)); //add the edge
            }
        }
    }
    for (i=0; i<G2->n; i++) { //add all G2 edges which haven't already been added
        for (j=0; j<LIST_LENGTH(G2->nodes[i]); j++) { //for each edge coming from that node
            tid = LIST_GET_ID(G2->nodes[i], j); //id of this edge's endpoint
            if ( G2->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                if ( !ALG_CONNECTED(G1, i, tid) ) { //if they're not connected in g1
                    ALG_ADD_EDGE(Go, i, tid, LIST_GET_VAL(G2->nodes[i], j)); //add the edge
                }
            }
        }
    }
    return Go; //return the OR'd graph
}


// OR 2 graphs together in place ("superimpose" them), that is, adds edges from G2
// to G1 which aren't already in G1
void ALG_ORIP(struct ALGRAPH * G1, struct ALGRAPH * G2)
{
    nid_t i; len_t j; //counters
    nid_t tid; //this id
    for (i=0; i<G2->n; i++) { //add all G2 edges which aren't in G1
        for (j=0; j<LIST_LENGTH(G2->nodes[i]); j++) { //for each edge coming from each of G2's nodes
            tid = LIST_GET_ID(G2->nodes[i], j); //id of this edge's endpoint
            if ( G2->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                if ( !ALG_CONNECTED(G1, i, tid) ) { //if they're not connected in g1
                    ALG_ADD_EDGE(G1, i, tid, LIST_GET_VAL(G2->nodes[i], j)); //add the edge
                }
            }
        }
    }
}


// Add the edge weights of two graphs together
void ALG_ADD(struct ALGRAPH * G1, struct ALGRAPH * G2)
{
    nid_t i; len_t j; //counters
    nid_t tid; //this id
    for (i=0; i<G2->n; i++) { //add all G2 edge values to G1 edges
        for (j=0; j<LIST_LENGTH(G2->nodes[i]); j++) { //for each edge coming from that node
            tid = LIST_GET_ID(G2->nodes[i], j); //id of this edge's endpoint
            if ( G2->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                ALG_ADD_VAL_TO_EDGE(G1, i, tid, LIST_GET_VAL(G2->nodes[i], j)); //add the value
            }
        }
    }
}



// Thresholds a graph.  That is, only keeps edges which have weights above some threshold th
void ALG_THRESH(struct ALGRAPH * G, struct ALGRAPH * Go, edge_t th)
{
    nid_t i; len_t j; //counters
    nid_t tid; //this id
    edge_t tw; //this weight
    for (i=0; i<G->n; i++) { //for each node
        for (j=0; j<LIST_LENGTH(G->nodes[i]); j++) { //for each edge coming from that node
            tid = LIST_GET_ID(G->nodes[i], j); //id of this edge's endpoint
            if ( G->dir || i<tid ) { //only add lesser id -> larger id in undirected graphs
                tw = LIST_GET_VAL(G->nodes[i], j); //weight of this edge in g
                if ( tw > th ) { //if this edge is above the threshold
                    ALG_ADD_EDGE(Go, i, tid, tw); //add the edge to output graph
                }
            }
        }
    }
}


