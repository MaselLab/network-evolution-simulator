#include "netsim.h"
#include "priority-queue.h"

/* add an element to the queue with an associated priority */
void insert_with_priority(FixedEvent **queue, FixedEvent **end, int cell, float time) {
  add_fixed_event(cell, time, queue, end);
  if (queue == NULL) 
    printf("queue is empty after adding an event\n");
}

/* remove the element from the queue that has the highest priority,
   and return it (also known as "PopElement(Off)", or "GetMinimum") */
float get_next(FixedEvent **queue, FixedEvent **end, int *cell) {
  float retval = 999.0;

  if (queue != NULL) {
    retval = (*queue)->time;
    *cell = (*queue)->geneID;
    delete_fixed_event_start(queue, end);
  } else {
    printf("queue is empty!\n");
  }
  return retval;
}

/* look at the element with highest priority without removing it */
float peek_at_next(FixedEvent *queue, FixedEvent *end, int *cell) {
  //
}
