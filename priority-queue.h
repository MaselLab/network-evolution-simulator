/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
//#define PriorityQueue FixedEvent

/* add an element to the queue with an associated priority */
extern int insert_with_priority(FixedEvent **, FixedEvent **, int, float);

/* remove the element from the queue that has the highest priority,
   and return it (also known as "PopElement(Off)", or "GetMinimum") */
extern float get_next(FixedEvent **, FixedEvent **, int *);

extern void delete_element(FixedEvent **, FixedEvent **, int, float);

/* look at the element with highest priority without removing it */
extern float peek_at_next(FixedEvent *, FixedEvent *, int *);

/* This implementation stores the binary heap in a 1 dimensional array. */

/*  Author: Shane Saunders 
    http://www.cosc.canterbury.ac.nz/tad.takaoka/alg/heaps/heaps.html 
    Modified by: Alex Lancaster (2008)
*/

/*** Structure Definitions ***/

/* Frontier set items in Dijkstra's algorithm stored in the binary heap
 * structure.
 */
typedef struct {
    int item;  /* vertex number is used for the item. */
    float key;  /* distance is used as the key. */
} bheap_item_t;

/* Binary heap structure for frontier set in Dijkstra's algorithm.
 * a[] - stores (distance, vertex) pairs of the binary heap.
 * p[] - stores the positions of vertices in the binary heap a[].
 * n - is the size of the binary heap.
 */
typedef struct {
  bheap_item_t *a;
  int *p;
  int n;
  long key_comps;
  
  int count_ops;

} bheap_t;

extern int insert_with_priority_heap(bheap_t *, int, float);

/* remove the element from the queue that has the highest priority,
   and return it (also known as "PopElement(Off)", or "GetMinimum") */
extern float get_next_heap(bheap_t *, int *, int *);

extern void delete_element_heap(bheap_t *, int);

/*** Function prototypes. ***/

/* Binary heap functions. */

/* bh_alloc() allocates space for a binary heap of size n and initialises it.
 * Returns a pointer to the binary heap.
 */
bheap_t *bh_alloc(int n);

/* bh_free() frees the space taken up by the binary heap pointed to by h.
 */
void bh_free(bheap_t *h);

/* bh_min() returns the item with the minimum key in the binary heap pointed to
 * by h.
 */
int bh_min(bheap_t *h);

/* bh_insert() inserts an item and its key value into the binary heap pointed
 * to by h.
 */
void bh_insert(bheap_t *h, int item, float key);

/* bh_delete() deletes an item from the binary heap pointed to by h.
 */
void bh_delete(bheap_t *h, int item);

/* bh_decrease_key() decreases the value of 'item's key and then performs
 * sift-down until 'item' has been relocated to the correct position in the
 * binary heap.
 */
void bh_decrease_key(bheap_t *h, int item, long new_key);

void bh_dump(bheap_t *h);
