/*  FAIRSPLIT.c

    Fair-split decomposition and application functions

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FAIRSPLIT.h"

/*****************************************************************************
                                SORTING
*****************************************************************************/

// Get value of a point with ind=ni in dimension di
val_t GET_VAL(FS_LNODE ** DL, val_t * V, int d, int di, fsp_t ni)
{
    return *(V+d*DL[di][ni].data+di); //value of point ni in dimension di
}

// Set DL[di][ni] = B, updating cross-pointers in adjacent lists
void SET_EQUAL_CPTRS(FS_LNODE ** DL, FS_LNODE * B, int d, int di, fsp_t ni)
{
    // set DL[di][ni] = B
    DL[di][ni].data = B->data;
    DL[di][ni].adj_n = B->adj_n; //keep cross refs
    DL[di][ni].adj_p = B->adj_p;

    // set last dim's 'next adj' cross pointer to here
    DL[(d+di-1)%d][B->adj_p].adj_n = ni;

    // set next dim's 'previous adj' cross pointer to here
    DL[(d+di+1)%d][B->adj_n].adj_p = ni;
}

// Set B = DL[di][ni]
void SET_EQUAL(FS_LNODE * B, FS_LNODE ** DL, int di, fsp_t ni)
{
    B->data = DL[di][ni].data; //pointer to data in V
    B->adj_p = DL[di][ni].adj_p; //cross-pointer to point in previous dimension's list
    B->adj_n = DL[di][ni].adj_n; //cross-pointer to point in next dimension's list
}

// Set objs in range A = B
void FS_CopyBack(FS_LNODE ** DL, FS_LNODE * B, int d, int di, fsp_t a, fsp_t b)
{
    int i;
    for (i=a; i<b; i++) {
        SET_EQUAL_CPTRS(DL, &B[i], d, di, i);
    }
}

void FS_TopDownMerge(FS_LNODE ** DL, val_t * V, FS_LNODE * B, int d, int di, fsp_t a, fsp_t m, fsp_t b)
{
    int i0=a, i1=m, j;
    for (j=a; j<b; j++) { //while we haven't sorted the whole run,
        // if left run is nonexhausted, and the right run IS, or smallest L run item is LEQ smallest R run item
        if (i0 < m && (i1 >= b || GET_VAL(DL, V, d, di, i0) <= GET_VAL(DL, V, d, di, i1))) {
            SET_EQUAL(&B[j], DL, di, i0); //B[j] = A[i0]
            i0++; //move onto next item of Left run
        } else { //right run has next smallest element
            SET_EQUAL(&B[j], DL, di, i1); //B[j] = A[i1]
            i1++; //move onto next item of right run
        }
    }
}

void FS_TopDownMergeSort(FS_LNODE ** DL, val_t * V, FS_LNODE * B, int d, int di, fsp_t a, fsp_t b)
{
    fsp_t m;
    if (b-a < 2) { return; } //only one element
    m = (a+b)/2; //find midpoint
    FS_TopDownMergeSort(DL, V, B, d, di, a, m); //sort first half recursively
    FS_TopDownMergeSort(DL, V, B, d, di, m, b); //sort second half recursively
    FS_TopDownMerge(DL, V, B, d, di, a, m, b); //merge the two lists
    FS_CopyBack(DL, B, d, di, a, b); //copy merged back into A
}

// Sort objects in dl array, setting cross pointers
void FS_MergeSort(FS_LNODE ** DL, val_t * V, int d, int di, fsp_t n)
{
    FS_LNODE * B;
    B = malloc(n*sizeof(FS_LNODE));
    FS_TopDownMergeSort(DL, V, B, d, di, 0, n); //call the recursive mergesort algo
    free(B); //free from mem
}



/*****************************************************************************
                                PRINTING FUNCTIONS
*****************************************************************************/

// Prints a point given the tree and an ind
void to_string_tree_node(FS_TREE * T, fsp_t id)
{
    int di;
    if (T->tree[id].c2 < 0) { //if leaf,
        printf("%d: P = %d: ( ", id, T->tree[id].c1);
        for (di=0; di<T->d; di++) { printf("%f, ", *(T->V+T->d*T->tree[id].c1+di)); }
        printf(" )\n");

    } else { //not leaf
        printf("%d: R = [ ", id); //placeholder corresponding to internal nodes
        for (di=0; di<T->d; di++) { printf("[%f, %f], ", T->tree[id].Ra[di], T->tree[id].Rb[di]); }
        printf(" ]\n");
    }
}


// Worker function for the to_string method
void to_string_recurse(FS_TNODE * tree, fsp_t id, val_t * V, int d, int rl)
{
    int i, di;
    for (i=0; i<rl; i++) printf("  "); //indentation
    if (tree[id].c2 < 0) { //if leaf,
        printf("P = %d: ( ", tree[id].c1);
        for (di=0; di<d; di++) { printf("%f, ", *(V+d*tree[id].c1+di)); }
        printf(" )\n");
    } else { //not leaf
        printf("R = [ "); //placeholder corresponding to internal nodes
        for (di=0; di<d; di++) { printf("[%f, %f], ", tree[id].Ra[di], tree[id].Rb[di]); }
        printf(" ]\n");
        to_string_recurse(tree, tree[id].c1, V, d, rl+1); //recursively print children
        to_string_recurse(tree, tree[id].c2, V, d, rl+1); //recursively print children
    }
}

// Depth-first traverse, assigning parent pointers (s=self ind, p=parent ind)
void FS_TREE_TO_STRING(FS_TREE * T)
{
    to_string_recurse(T->tree, T->root, T->V, T->d, 0); //recursively print the tree
}



/*****************************************************************************
                                STRUCTURE CLEARING
*****************************************************************************/

// Recursively clear node bounds
void clear_bounds(FS_TNODE * T, fsp_t ni)
{
    if (T[ni].c2 >= 0) { //if internal node,
        clear_bounds(T, T[ni].c1);
        clear_bounds(T, T[ni].c2);
    }
    free(T[ni].Ra); //clear lower bounds
    free(T[ni].Rb); //clear upper bounds
}

// Clear a wsp tree from memory
void CLEAR_FS_TREE(FS_TREE * T)
{
    clear_bounds(T->tree, T->root); //clear the Rs
    free(T->tree);
    free(T);
}


/********************************************************************
                              FS DECOMPOSITION
********************************************************************/

/* Recursively create FS tree
   Inputs:
    T - pointer to first node in tree's node array
    ti - pointer to ind of first empty spot in T
    DL - the doubly-linked lists
    V - N-by-D array of values
    d - number of dimensions
    first - vector of inds of first element in list for each dimension
    last - vector of inds of last element in list for each dimension
*/
fsp_t FS_RECURSE(FS_TNODE * T, fsp_t * ti, FS_LNODE ** DL, val_t * V, int d, fsp_t * first, fsp_t * last)
{
    // Declare variables
    int di; //dimension counter
    float mdw, tdw; //max dimension width, this dimension's width
    int mdi; //max dimension ind
    fsp_t fbi, fei, i; //index of search from beginning, and end, total step counter
    fsp_t adp; //pointer to this point in adjacent dimension
    fsp_t * firsts1, * firsts2, * lasts1, * lasts2; //firsts and lasts for new lists
    int td; //this dimension
    
    // Check to see if base case (list of length 1)
    if (first[0] == last[0]) { //if we've got only a single point, (first==last in any dimension)
        T[*ti].c1 = DL[0][first[0]].data; //store pointer to this point's data
        T[*ti].c2 = -1; //store null value to signify leaf
        //T[*ti].wspl = NULL; //no wsp list yet
        T[*ti].Ra = NULL; //points dont have bounds, and re-storing data would be waste of O(n) space; it's already in V
        T[*ti].Rb = NULL;
        (*ti)++; //we used the spot at ti, so increment next avail spot (NOTE: increment operator binds tighter than the f@%&ing dereference operator, so can't do "*ti++"!!!)
        return *ti-1; //return pointer to this tree node
    }

    // 1. Find widest dimension
    mdw = 0; //max dimension width,
    mdi = 0; //max dimension ind
    for (di=0; di<d; di++) { // for each dimension
        tdw = GET_VAL(DL, V, d, di, last[di]) - GET_VAL(DL, V, d, di, first[di]); //width of dimension di
        if ( tdw > mdw ) { //if widest dimension so far,
            mdw = tdw; //store new max dim width
            mdi = di; //this dim is widest so far
        }
    }

    // 2. Walk thru widest dim's list from each end to find midpoint
    //assigning .half=0 or 1 as you go to points in adjacent dims
    i=0; //index of search from beginning, and end, total step counter
    fbi = first[mdi]; //start from beginning of this list in widest dim
    fei = last[mdi]; //and last
    while ( 1 ) { //while we haven't met in the middle
        if (i%2==0) { //step from beginning half the time
            adp = DL[mdi][fbi].adj_n; //pointer to this point in adjacent dimension
            for (di=1; di<d; di++) { // set flags for this point in adjacent dimension's lists
                DL[(mdi+di)%d][adp].half = (char) 0; //set flag to first half
                adp = DL[(mdi+di)%d][adp].adj_n; //pointer to this point in dim after that
            }
            if ( fei==fbi ) break; //quit if we've met in middle
            fbi = DL[mdi][fbi].next; //step forward
        } else { //step from end other half of the time
            adp = DL[mdi][fei].adj_n; //pointer to this point in adjacent dimension
            for (di=1; di<d; di++) { // set flags for this point in adjacent dimension's lists
                DL[(mdi+di)%d][adp].half = (char) 1; //set flag to first half
                adp = DL[(mdi+di)%d][adp].adj_n; //pointer to this point in dim after that
            }
            if ( fei==fbi ) break; //quit if we've met in middle
            fei = DL[mdi][fei].prev; //step backwards
        }
        i++; //increment step counter
    }

    // 3. Split widest dim's list (node we met at is last element in 1st half's list)
    // Create new firsts and lasts lists for both lists
    // Then set this list's midpoint to have null next pointer, and mid+1 to have null prev pointer
    //fsp_t firsts1[d], firsts2[d], lasts1[d], lasts2[d]; //firsts and lasts for new lists
    firsts1 = malloc(d*sizeof(fsp_t)); //firsts and lasts for new lists
    firsts2 = malloc(d*sizeof(fsp_t));
    lasts1 = malloc(d*sizeof(fsp_t));
    lasts2 = malloc(d*sizeof(fsp_t));
    if (i%2==0) { // we met in the exact middle of the list
        fei = DL[mdi][fbi].next; //point to first element of 2nd half's list
    } else { // list had even #, met after middle
        fbi = DL[mdi][fbi].prev; //point to last element of 1st half's list
    }
    DL[mdi][fbi].next = -1; //make fbi last element of 1st half's list
    DL[mdi][fei].prev = -1; //make fei first element of 2nd half's list
    firsts1[mdi] = first[mdi]; //same first element for this dim
    firsts2[mdi] = fei; //2nd list now starts at fei
    lasts1[mdi] = fbi; //first half's list now ends at fbi
    lasts2[mdi] = last[mdi]; //second half's list ends at same place

    // 4. Spit other dim's lists according to widest dim's split
    //Walk through each adjacent dim's list, re-setting pointers according to the 0/1 flags, to create 2 separate lists
    for (di=1; di<d; di++) { // for each dimension
        td = (mdi+di)%d; //this dimension
        fbi = -1; // now is pointer to most recent element in first half
        fei = -1; // now is pointer to most recent element in second half
        i = first[td]; //first in this list overall
        while ( fbi!=last[td] && fei!=last[td] ) { //while we haven't hit the end of the list
            if (DL[td][i].half==0) { //this node to go in first half's list
                DL[td][i].prev = fbi; //point back to most recent element in list 1
                if (fbi >= 0) { //if this isn't the first element in the list,
                    DL[td][fbi].next = i; //previous element of list 1 points next to here
                } else { //this is the first element of list 1
                    firsts1[td] = i; //store pointer to first element of list 1
                }
                fbi = i; //store most recent entry in list 1
            } else { //this node to go in 2nd half's list
                DL[td][i].prev = fei; //point back to most recent element in list 2
                if (fei >= 0) { //if this isn't the first element in the list,
                    DL[td][fei].next = i; //previous element of list 2 points next to here
                } else { //this is the first element of list 1
                    firsts2[td] = i; //store pointer to first element of list 1
                }
                fei = i; //store most recent entry in list 1
            }
            i = DL[td][i].next; //move on to next node in parent list
        }
        DL[td][fbi].next = -1; //last element of first list points next to null
        DL[td][fei].next = -1; //last element of 2nd list points next to null
        lasts1[td] = fbi; //store last element of 1st list
        lasts2[td] = fei; //store last element of 2nd list
    }

    // 5. Recursively split each list, getting pointers to Tree nodes for each
    fbi = FS_RECURSE(T, ti, DL, V, d, firsts1, lasts1); //get pointer to tree node for 1st half
    fei = FS_RECURSE(T, ti, DL, V, d, firsts2, lasts2); //get pointer to tree node for 2nd half

    // 6. Create new tree node, and link to its children, store bounding box if internal node, increment value at ti
    T[*ti].c1 = fbi; //store pointer to child 1 (1st half of list)
    T[*ti].c2 = fei; //store pointer to child 2 (2nd half of list)
    //T[*ti].wspl = NULL; //no wsp list yet
    T[*ti].Ra = malloc(d*sizeof(val_t)); //store lower bounds for dimensions
    T[*ti].Rb = malloc(d*sizeof(val_t)); //store upper bounds for dimensions
    for (di=0; di<d; di++){
        T[*ti].Ra[di] = GET_VAL(DL, V, d, di, first[di]); //lower bound for this dim
        T[*ti].Rb[di] = GET_VAL(DL, V, d, di, last[di]); //lower bound for this dim
    }
    (*ti)++; //we used the spot at ti, so increment next avail spot (NOTE: increment operator binds tighter than the damn dereference operator, so can't do "*ti++")
    free(firsts1); free(firsts2); free(lasts1); free(lasts2); //free the arrays
    return *ti-1; //return pointer to this tree node
}

// Get value of a point with ind=ni in dimension di
val_t GET_VAL_FROM_IND(val_t * V, int d, int di, fsp_t ni)
{
    return *(V+d*ni+di); //value of point ni in dimension di
}

// Depth-first traverse, assigning parent pointers (s=self ind, p=parent ind)
// Returns the number of points contained in the bounding box for that node
fsp_t assign_parents(FS_TNODE * T, fsp_t s, fsp_t p)
{
    T[s].p = p; //assign pointer to parent for this node
    if (T[s].c2 >= 0) { //if internal node (not leaf)
        assign_parents(T, T[s].c1, s); //make child point to me
        assign_parents(T, T[s].c2, s); //make 2nd child point to me, summing points
    }
    return 1;
}


/* Find the fair-split tree for some data, returning the binary tree
   Inputs:
    V - N-by-D array of values
    n - number of points
    d - number of dimensions
    root - will return the index of the root of the tree (also is the number of elements in the tree)
*/
FS_TREE * FAIR_SPLIT(val_t * V, fsp_t n, int d)
{
    int di; //counter thru dimensions
    fsp_t ni; //counter thru datapoints
    FS_LNODE ** DL;
    FS_TREE * T;
    fsp_t * firsts, * lasts; //pointers to first and last elements of lists
    fsp_t ti; //first available element in FST

    // 1. Construct doubly-linked list for each dimension of data
    DL = malloc(d*sizeof(FS_LNODE *)); //a d-length array of pointers to doubly-linked lists
    for (di=0; di<d; di++) { // for each dimension, (NOTE: this is 630mb for 3m spikes and 10dims, in 10 63mb chunks)
        if (!(DL[di] = malloc(n*sizeof(FS_LNODE)))) { //allocate that dim's doubly-linked list
            printf("ERROR in FS_FAIR_SPLIT: cannot allocate memory for doubly-linked lists\n");
        }
        for (ni=0; ni<n; ni++) { //for each node in that new list
            DL[di][ni].data = ni; //pointer to data
            DL[di][ni].adj_n = ni; //pointers to same point in adjacent dimension's lists
            DL[di][ni].adj_p = ni;
        }
    }

    // 2. Sort each dimension's list
    for (di=0; di<d; di++) { // for each dimension,
        FS_MergeSort(DL, V, d, di, n); //sort this dimension's list, keeping track of cross-pointers
    }

    // Set forward/back pointers within each dimension's list (so list pointers go in sorted dir)
    for (di=0; di<d; di++) { // for each dimension,
        for (ni=0; ni<n; ni++) { //for each node
            DL[di][ni].prev = ni-1; //pointer to previous node in list
            DL[di][ni].next = ni+1; //pointer to next node in list
        }
        DL[di][n-1].next = -1; //set last node's next pointer to null (1st's prev will be set in loop)
    }

    // 3. Allocate space for the binary tree
    T = malloc(sizeof(FS_TREE));
    if (!(T->tree = malloc((n+n-1)*sizeof(FS_TNODE)))) { //allocate space, should only be 72mb for 3m spikes
        printf("ERROR in FS_FAIR_SPLIT: cannot allocate memory for FS_TREE's list\n");
    }
    T->d = d; //save the number of dimensions in the tree struct
    T->n = n+n-1; //number of nodes in the tree
    T->V = V; //save pointer to data array in the tree struct
    ti=0; //first available element in FST

    firsts = malloc(d * sizeof(fsp_t));
    lasts = malloc(d * sizeof(fsp_t));

    for (di=0; di<d; di++) { // for each dimension,
        firsts[di] = 0; //sorted, so set first to first
        lasts[di] = n-1; //likewise
    }

    // 4. Recursively form the FS tree
    T->root = FS_RECURSE(T->tree, &ti, DL, V, d, firsts, lasts); //form the tree and save pointer to root

    // 5. Free the linked lists from memory
    for (di=0; di<d; di++) { // for each dimension,
        free(DL[di]); //free that malloc'd block
    }

    // 6. Depth-first traverse the tree, assigning parent backpointers
    assign_parents(T->tree, T->root, -1);

    // free the doubly-linked lists
    free(DL);
    free(firsts);
    free(lasts);

    // Return pointer to the tree
    return T;

}


/*********************************************************************
                CHEAP MINIMUM SPANNING TREE USING ONLY THE TREE
*********************************************************************/

// Calculates the distance between 2 nodes (given the TREE ids)
val_t calcdist(FS_TREE * T, fsp_t a, fsp_t b)
{
    fsp_t di;
    val_t av, bv, ddiff, dsum=0;
    if (T->tree[a].c2 < 0) { //if a is leaf
        if (T->tree[b].c2 < 0) { // both are leaf
            for (di=0; di<T->d; di++) {
                av = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[a].c1);
                bv = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[b].c1);
                ddiff = bv-av;
                dsum += ddiff*ddiff;
            }
        } else { //b is internal node, a is leaf
            for (di=0; di<T->d; di++) {
                av = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[a].c1);
                bv = 0.5*(T->tree[b].Ra[di]+T->tree[b].Rb[di]); //midpoint of bounds
                ddiff = bv-av;
                dsum += ddiff*ddiff;
            }
        }
    } else { //a is internal node
        if (T->tree[b].c2 < 0) { // b is leaf, a is internal node
            for (di=0; di<T->d; di++) {
                av = 0.5*(T->tree[a].Ra[di]+T->tree[a].Rb[di]); //midpoint of bounds
                bv = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[b].c1);
                ddiff = bv-av;
                dsum += ddiff*ddiff;
            }
        } else { //both are internal nodes
            for (di=0; di<T->d; di++) {
                av = 0.5*(T->tree[a].Ra[di]+T->tree[a].Rb[di]); //midpoint of bounds
                bv = 0.5*(T->tree[b].Ra[di]+T->tree[b].Rb[di]); //midpoint of bounds
                ddiff = bv-av;
                dsum += ddiff*ddiff;
            }
        }
    }
    return sqrt(dsum); //return distance between points/node centers
}


// Connect the closest children of two internal nodes
void connect_closest_children(FS_TREE * T, Graph * G, fsp_t a, fsp_t b)
{
    while (T->tree[a].c2>0 || T->tree[b].c2>0) { //find closest points: while we still have an internal node,
        if ( T->tree[a].c2 > 0 ) { //if a isn't leaf
            if (calcdist(T, T->tree[a].c1, b) < calcdist(T, T->tree[a].c2, b)) { //if c1 closer than c2
                a = T->tree[a].c1; //move down the tree
            } else { //c2 is closer to B
                a = T->tree[a].c2; //move down the tree
            }
        }
        if ( T->tree[b].c2 > 0 ) { //if b isn't leaf
            if (calcdist(T, T->tree[b].c1, a) < calcdist(T, T->tree[b].c2, a)) { //if c1 closer than c2
                b = T->tree[b].c1; //move down the tree
            } else { //c2 is closer to B
                b = T->tree[b].c2; //move down the tree
            }
        }
    }
    ALG_ADD_EDGE(G, T->tree[a].c1, T->tree[b].c1, calcdist(T, a, b)); //add edge to algraph
}


// recursion function for FS_AEMST
void aemst_recurse(FS_TREE * T, Graph * G, fsp_t nid)
{
    if (T->tree[nid].c2 > 0) { //as long as I am not a leaf,
        connect_closest_children(T, G, T->tree[nid].c1, T->tree[nid].c2); //recursively connect closest two points of my children
        aemst_recurse(T, G, T->tree[nid].c1); //recurse on children
        aemst_recurse(T, G, T->tree[nid].c2);
    }
}


// Calculates the (very) Approximate Euclidian Spanning Graph using a cheap algorithm
void FS_AEMSG(FS_TREE * T, Graph * G)
{
    aemst_recurse(T, G, T->root); //recursively connect shortest edge separating all sibling internal nodes
}


/*****************************************************************************
                           CHEAP KNN USING ONLY THE TREE
*****************************************************************************/

/*  Adds an element to the KNN list
    T - the FS tree
    NNd - list of knn dists
    NNid - list of knn ids
    cfn - current distance of the furthest neighbor
    ck - current length of the list
    k - number of nearest neighbors desired
    nid - id of new node
    ndist - distance from p to current node
*/
void add_ann_to_list(val_t * NNd, fsp_t * NNid, val_t * cfn, int * ck, int k, fsp_t nid, val_t ndist)
{
    int i=0;
    val_t odist; //dist to last point in list
    fsp_t oid; //id of last point in list
    if ( *ck < 1 ) { //if the list is empty
        NNd[0] = ndist; //add this point as first element
        NNid[0] = nid;
        *cfn = ndist; //record the dist as the current furthest element of the list
        *ck = 1; //list now contains 1 element
    } else { //there are already elements in the list
        if ( *ck < k) { //if not full
            NNd[*ck] = 100*ndist; //last element is larger than new one
            NNid[*ck] = -1; //null id value for last element
            (*ck)++; //then the list is getting bigger
        }
        while ( i < *ck && ndist > NNd[i] ) i++; //find where to insert new point
        if ( i < k && NNid[i] == nid ) { //if this point has already been added to the list,
            (*ck)--; //jk no we didn't make the list longer
        } else {
            while ( i < *ck ) { //while we haven't reached the end of the list
                odist = NNd[i]; //swap 'registers' with list elements
                oid = NNid[i];
                NNd[i] = ndist;
                NNid[i] = nid;
                ndist = odist;
                nid = oid;
                i++; //move on to next point
            }
            *cfn = NNd[--i]; //save distance of new furthest node
        }
    }
}


/* Calculates the minimum possible distance between a d-box or point and a point
    T - fair split tree
    a - the d-box or point
    b - the point for which to calculate the min dist from a
*/
val_t calcmindist(FS_TREE * T, fsp_t a, fsp_t b)
{
    fsp_t di;
    val_t av, bv, ddiff, dsum=0;
    if (T->tree[a].c2 < 0) { //if a is a point (leaf)
        for (di=0; di<T->d; di++) {
            av = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[a].c1);
            bv = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[b].c1);
            ddiff = bv-av;
            dsum += ddiff*ddiff;
        }
    } else { //a is internal node
        for (di=0; di<T->d; di++) {
            bv = GET_VAL_FROM_IND(T->V, T->d, di, T->tree[b].c1); //get b value in this dim
            if ( bv > T->tree[a].Rb[di] ) { //b is above upper bound
                ddiff = bv - T->tree[a].Rb[di]; //dist from bound to point
            } else if ( bv < T->tree[a].Ra[di] ) { //b is below lower bound
                ddiff = T->tree[a].Ra[di] - bv; //dist from bound to point
            } else { //b is within the bounds
                ddiff = 0; //so dist could be zero
            }
            dsum += ddiff*ddiff;
        }
    }
    return sqrt(dsum); //return distance between points/node centers
}


/* Returns the TREE ids and VALUE ids in dfs-order for a fair-split tree

*/
void dfs_list(FS_TREE * T, fsp_t nid, fsp_t * L, fsp_t * VL, fsp_t * li)
{
    if (T->tree[nid].c2 < 0) { //if leaf,
        *(VL+*li) = T->tree[nid].c1; //save this node's VALUE id
        *(L+*li) = nid; //save this node's tree id
        (*li)++; //increment the index of where we're at in the list
    } else { //internal node
        dfs_list(T, T->tree[nid].c1, L, VL, li); //recurse on children in DFS order
        dfs_list(T, T->tree[nid].c2, L, VL, li); //recurse on children in DFS order
    }
}


/*  Recurses on the tree to find the points which are KNN of p.
    T - the FS tree
    KNNv - list of knn dists
    KNNid - list of knn ids
    cfn - current distance of the furthest neighbor
    ck - current length of the list
    k - number of nearest neighbors desired
    nid - id of the node for whom we're finding the KNN list (TREE id)
    cid - id of the current node we're on in the tree
    cdist - mindist possible between cid and nid
*/
void aknn_recurse(FS_TREE * T, val_t * NNd, fsp_t * NNid, val_t * cfn, int * ck, int k, fsp_t nid, fsp_t cid, val_t cdist)
{
    val_t a, b; //to store dists to R and L child
    if (T->tree[cid].c2 < 0) { //if current node is leaf,
        if (T->tree[nid].c1 != T->tree[cid].c1) { //don't insert me into my own list
            add_ann_to_list(NNd, NNid, cfn, ck, k, T->tree[cid].c1, cdist); //insert (value/"true") id and dist to sorted KNN list
        }
    } else { //cid is internal node
        a = calcmindist(T, T->tree[cid].c1, nid); //get dists to children
        b = calcmindist(T, T->tree[cid].c2, nid);
        if ( a < b ) { //do a first if closer
            if (a < *cfn) { //if a is closer than any in the list,
                aknn_recurse(T, NNd, NNid, cfn, ck, k, nid, T->tree[cid].c1, a);
                if (b < *cfn) { //then check other child
                    aknn_recurse(T, NNd, NNid, cfn, ck, k, nid, T->tree[cid].c2, b);
                }
            }
        } else { //b is closer
            if (b < *cfn) { //if a is closer than any in the list,
                aknn_recurse(T, NNd, NNid, cfn, ck, k, nid, T->tree[cid].c2, b);
                if (a < *cfn) { //then check other child
                    aknn_recurse(T, NNd, NNid, cfn, ck, k, nid, T->tree[cid].c1, a);
                }
            }
        }
    }
}


/*  Computes the KNN using the fair split tree
    T - the fair split tree
    n - number of points
    NNid - list of ids of nearest neighbors for each point (n*k-length vector)
    NNd - list of dists to nearest neighbors for each point (n*k-length vector)
    k - number of nearest neighbors to find
*/
void FS_AKNN(FS_TREE * T, fsp_t n, fsp_t * NNid, val_t * NNd, int k)
{
    // Declare variables
    fsp_t * dfsl, * dfsl_v;
    fsp_t li;
    fsp_t i, a, b; int j; //counters
    val_t td;
    val_t * tvals; //pointer to this point's dimension values
    int k2;
    val_t * MDs;
    int * CK;
    fsp_t print_div;

    // Get the DFS list
    dfsl = malloc(T->n * sizeof(fsp_t));
    dfsl_v = malloc(T->n * sizeof(fsp_t));
    li = 0; //index for the list
    dfs_list(T, T->root, dfsl, dfsl_v, &li);  //get the dfs list

    // Initialize temp KNN lists
    k2 = ceil(k/2.0);
    MDs = calloc(n, sizeof(val_t));
    CK = calloc(n, sizeof(int));
    for (i=0; i<n; i++) { //for each point
        *(MDs+dfsl_v[i]) = 0; //init to 0
        *(CK+dfsl_v[i]) = 0; //init to 0
        if (i-k2 < 0) { //if off lower edge of list
            a = 0; b = k+1;
        } else if (i+k2 >= n) { //if off upper end of list
            a = n-k-1; b = n-1;
        } else { //within the bounds
            a = i-k2;  b = a+k2+1;
        }
        for (j=a; j<=b; j++) { //for each of its nearest neighbors,
            if (j != i) { //don't add me to my own list
                td = calcmindist(T, dfsl[j], dfsl[i]); //dist from this node to this element of its list
                add_ann_to_list(NNd+k*dfsl_v[i], NNid+k*dfsl_v[i], MDs+dfsl_v[i], CK+dfsl_v[i], k, dfsl_v[j], td); //insert (value/"true") id and dist to sorted KNN list
            }
        }
    }

    // Run search to find better options for each point
    print_div = n/10;
    for (i=0; i<n; i++) { //for each point, update list with better options than dfs-adjacents
        //if (i%print_div==0) { //print every once in a while
        //    mexPrintf("  %.1f percent complete\n", 100 * (float) i / (float) n);
        //    mexEvalString("drawnow;");
        //}
        aknn_recurse(T, NNd+k*dfsl_v[i], NNid+k*dfsl_v[i], MDs+dfsl_v[i], CK+dfsl_v[i], k, dfsl[i], T->root, 0);
    }

    // Garbage
    free(dfsl);
    free(dfsl_v);
    free(MDs);
    free(CK);
}


/* Computes the mutual K-nearest neighbors and adds the edges to graph G
*/
void FS_AMKNN(FS_TREE * T, Graph * G, fsp_t n, int k)
{
    // Find the k nearest neighbors for each point
    fsp_t i; int j, m; //counters
    //fsp_t IDs[n]; //ids of nearest neighbors
    fsp_t * IDs = malloc(n*sizeof(fsp_t)); //ids of nearest neighbors
    fsp_t tnn; //current nearest neighbor
    fsp_t * NNid = calloc(n*k, sizeof(fsp_t)); //nearest neighbor ids
    val_t * NNd = calloc(n*k, sizeof(val_t)); //nearest neighbor distances
    if (NNid == NULL || NNd == NULL) {
        printf("ERROR in FS_AMKNN: Could not allocate enough space for nearest neighbors");
        exit(1);
    }

    // Find the K-nearest neighbors
    FS_AKNN(T, n, NNid, NNd, k);

    // Now find nearest neighbors who are mutual
    for (i=0; i<n; i++) { //for each point, check for mutual KNNs
        for (j=0; j<k; j++) { //for each element of their KNN list
            tnn = *(NNid+i*k+j); //value id of this neighbor
            if (i < tnn) { //only do lesser -> greater id edges
                for (m=0; m<k; m++) { //check if they're mutual neighbors
                    if ( *(NNid+tnn*k+m) == i ) { //if they are in each other's lists
                        ALG_ADD_EDGE(G, i, tnn, *(NNd+i*k+j)); //add the edge to MKNN graph
                        m += k; //exit loop since we found the neighbor
                    }
                }
            }
        }
    }

    // Garbage
    free(NNid);
    free(NNd);
    free(IDs);
}
