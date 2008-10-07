#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "netsim.h"
#include "random.h"
#include "priority-queue.h"

int main(int argc, char *argv[])
{
  FixedEvent *queue;
  FixedEvent *end;
  int cell;
  int i;
  long seed  = 28121;
  float next;

  queue = NULL;
  end = NULL;

  for (i=0; i < 10; i++) {
    insert_with_priority(&(queue), &(end), trunc(10*ran1(&seed)), ran1(&seed));
  }

  next = get_next(&queue, &end, &cell);
  printf("get next time=%.3f, cell =%d\n", next, cell);

  for (i=0; i < 2; i++) {
    insert_with_priority(&(queue), &(end), trunc(10*ran1(&seed)), ran1(&seed));
  }


  while (queue != NULL) {
  //for (i=0; i < 5; i++) {
    next = get_next(&queue, &end, &cell);
    printf("get next time=%.3f, cell =%d\n", next, cell);
  }

}
