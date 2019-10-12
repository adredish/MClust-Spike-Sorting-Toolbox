/*  Mergesort.c

    Sort an array of structs

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
*/

#include <stdlib.h>
#include <stdio.h>

// What structs/functions to use, so we can swap what type we're sorting easily
#include "ALGRAPH.h"
typedef struct ALG_EDGE Obj;
#define SET_EQUAL(O1, O2) es_set_edges_equal(O1, O2)
#define GET_VAL(O) es_get_edge_val(O)


// Set obj A = B
void CopyBack(Obj * A, Obj * B, int a, int b)
{
    int i;
    for (i=a; i<b; i++) {
        SET_EQUAL(&A[i], &B[i]);
    }
}

void TopDownMerge(Obj * A, Obj * B, int a, int m, int b)
{
    int i0, i1, j;
    i0=a; i1=m;
    for (j=a; j<b; j++) { //while we haven't sorted the whole run,
        // if left run is nonexhausted, and the right run IS, or smallest L run item is LEQ smallest R run item
        if (i0 < m && (i1 >= b || GET_VAL(&A[i0]) <= GET_VAL(&A[i1]))) {
            SET_EQUAL(&B[j], &A[i0]); //B[j] = A[i0]
            i0++; //move onto next item of Left run
        } else { //right run has next smallest element
            SET_EQUAL(&B[j], &A[i1]); //B[j] = A[i1]
            i1++; //move onto next item of right run
        }
    }
}

void TopDownMergeSort(Obj * A, Obj * B, int a, int b)
{
    int m;
    if (b-a < 2) { return; } //only one element
    m = (a+b)/2; //find midpoint
    TopDownMergeSort(A, B, a, m); //sort first half recursively
    TopDownMergeSort(A, B, m, b); //sort second half recursively
    TopDownMerge(A, B, a, m, b); //merge the two lists
    CopyBack(A, B, a, b); //copy merged back into A
}

// Sort objects in array A of length n, in-place
void MergeSort(Obj * A, int n)
{
    Obj * B;
    B = malloc(n*sizeof(Obj));
    TopDownMergeSort(A, B, 0, n); //call the recursive mergesort algo
    free(B); //free the extra list from mem
}


