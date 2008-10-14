/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
//#define PriorityQueue FixedEvent

/* add an element to the queue with an associated priority */
extern void insert_with_priority(FixedEvent **, FixedEvent **, int, float);

/* remove the element from the queue that has the highest priority,
   and return it (also known as "PopElement(Off)", or "GetMinimum") */
extern float get_next(FixedEvent **, FixedEvent **, int *);

extern void delete_element(FixedEvent **, FixedEvent **, int, float);

/* look at the element with highest priority without removing it */
extern float peek_at_next(FixedEvent *, FixedEvent *, int *);
