/* SPC.h

    Header for the superparamagnetic clustering functions

    Brendan Hasz
    haszx010@umn.edu
    Aug 2015
    David Redish Lab
    University of Minnesota, Twin Cities
*/

//#include <math.h>
#ifndef DEF_EDGELIST_LEN //if ALGRAPH has not been included yet,
    #include "ALGRAPH.h"
    typedef struct ALGRAPH Graph; //define graph datastruct to use
#endif

// Random number generators
#define randf ((float)rand()/(float)(RAND_MAX)) //random float
#define randi(a) (rand()%a)  //not a perfect random int generator, but good enough
#define randn ( sqrt(-2*log(randf)) * sin(2*M_PI*randf) )

// Rename graph functions so implementation swapping is easier
#define NEW_GRAPH(n) NEW_ALGRAPH(n, 0)
#define CLEAR_GRAPH(G) ALG_CLEAR(G)
#define EMPTY_GRAPH(G) ALG_EMPTY(G)
#define GRAPH_ADD_EDGE(G, i, j, w) ALG_ADD_EDGE(G, i, j, w)

// SPC function
//void SPC_APPROX(edge_t *, nid_t, int, int, nid_t, int, int *, edge_t **, nid_t **, nid_t **, nid_t **, int *, nid_t **, nid_t ***, nid_t **);

void SPC_HIERARCHICAL(edge_t *, nid_t, int, int, nid_t, int, int *, nid_t **, nid_t ***, nid_t **);

