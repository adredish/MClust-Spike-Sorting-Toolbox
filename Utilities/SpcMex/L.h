/*  L.h
 
    List data structure and utility functions.  
    Elements can only be added or deleted from end
  
    To create a new linked list:
       L = NEW_L(N);
       where N is the number of nodes initially allocated

    Functions:
    NEW_L - make a new linked list
    L_GET - get the value in the i-th node
    L_FIND - get the value in a node with specified id
    L_FIND_ID - get the id of a node with specified value
    L_PUSH - add a node to the end of the list
    L_POP - remove last node from list (get id + val)
    L_TO_STRING - print out a list
    L_CLEAR - clear a linked list data structure from memory

    Internal Functions:
    double_size - double the size of the allocated array

    Notes:
    float is used for the value type because double operations take 3-24 times
        as long on CUDA cores.  However, this may be fixed in future CUDA 
        systems, in which case you can just change the definition of val_t to 
        double.
    long int is used for the id datatype because normal ints have only a range
        up to 65k, and we might be working with graphs with millions of nodes.
    unsigned short int is used for the length type for the list because for 
        our purposes it is unlikely that any list will have length >65k

    Brendan Hasz
    haszx010@umn.edu
    Feb 2015
    David Redish Lab
    University of Minnesota, Twin Cities
 */

// Datatypes of struct elements
//typedef float val_t; //data type of list element values
typedef float val_t; //data type of list element values
typedef int nid_t; //data type of identifier values (depends, you'll need at least 32 bits, this will be 'int' on some machines/compilers, and 'long int' on others.  Start with long int, make an array of len 1000, and save it to file.  If the filesize is 8kb, then 'long int'=64 bits and you can use 'int'
typedef int len_t; //datatype of list length 

// Default null values of struct elements
#define NULL_ID -1 //value of an unassigned id
#define NULL_VAL 0.0 //value of an unassigned value

struct L //list head data structure
{
    struct LNODE * nodes; //array of allocated nodes
    len_t n; //number of nodes (not less than 0, and shouldn't be longer than 65k, so unsigned short int)
    len_t na; //number of allocated slots
};

struct LNODE //list node element data structure
{
    nid_t id; //id of the node
    val_t v; //value stored by this node
};

void double_size(struct L *);
struct L * NEW_L(len_t);
nid_t L_GET_ID(struct L *, len_t);
val_t L_GET_VAL(struct L *, len_t);
void L_SET_ID(struct L *, len_t, nid_t);
void L_SET_VAL(struct L *, len_t, val_t);
val_t L_FIND(struct L *, nid_t);
nid_t L_FIND_ID(struct L *, val_t);
struct LNODE * L_ADD(struct L *, nid_t, val_t);
int L_ADD_VAL_BY_ID(struct L *, nid_t, val_t);
int L_POP(struct L *, nid_t *, val_t *);
len_t L_LENGTH(struct L *);
void L_TO_STRING(struct L *);
void L_EMPTY(struct L *);
void L_CLEAR(struct L *);
struct LNODE * L_PUSH_ID(struct L *, nid_t, nid_t);
int L_POP_ID(struct L *, nid_t *, nid_t *);
val_t L_SUM_VALS(struct L *);

