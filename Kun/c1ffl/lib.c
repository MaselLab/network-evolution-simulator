/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 *
 * lib.c
 * 
 * This file mainly contains helper/support functions which aren't
 * necessarily specific to the model, such as
 * maintaining the linked list data structures or file I/O
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdlib.h>
#include <math.h>
#include <sys/stat.h>   /* for file/directory permissions used by mkdir() */
#include <errno.h>      /* for error codes */
#include "lib.h"

int sls_store(FixedEvent *i, 
	      FixedEvent **start, 
	      FixedEvent **last)
{
  FixedEvent *old, *p;

  int pos = 0;
  
  p = *start;
  if (!*last) { /* first element in list */
    i->next = NULL;
    *last = i;
    *start = i;
    return pos;
  }
  old=NULL;
  while (p) {
    if (p->time < i->time) {
      old = p;
      p = p->next;
      pos++;
    }
    else {
      if (old) { /* goes in the middle */
        old->next = i;
        i->next = p;
        return pos;
      } else {
	i->next = p; /* new first element */
	*start = i;
	return pos;
      }
    }
  }
  (*last)->next = i; /* put on end */
  i->next = NULL;
  *last = i;
  return pos;
}


//void delete_time_course(TimeCourse *start2)
//{
//  TimeCourse *info, *start;
//  
//  start = start2; 
//   
//  while (start){
//    info = start;
//    start = start->next;
//    free(info);
//  }
//}

/*append to the end*/      
void sls_store_end(FixedEvent *i, 
                   FixedEvent **start, 
                   FixedEvent **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

//void sls_store_end2(TimeCourse *i, 
//                    TimeCourse **start, 
//                    TimeCourse **last)
//{
//  i->next = NULL;
//  if (!*last) *start = i;
//  else (*last)->next = i;
//  *last = i;
//}
