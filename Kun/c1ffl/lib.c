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


// System panics need to be handled better, but this will do for the moment.
// --blr
//static void
//panic(char *file, int line, char *msg){
//	// Notes:  
//	//
//	// The fflush is (strictly speaking) redundant as abort() will also 
//	//   cause all open files to be flushed.  It has been retained for 
//	//   readability.
//	//
//	// abort() bypasses calls to atexit() and onexit() handlers.  The 
//	//   assumption here is that things have gone sufficiently wrong
//	//   that calling those functions would make a bigger mess than 
//	//   avoiding them.  Perhaps "fatal()" could be implemented as an
//	//   additional logging function with a call to "exit()" instead  
//	//   of "abort()" if this proves problematic.
//	//
//	// Example:
//	//    panic( __FILE__, __LINE__, "Can't open kdis.txt." );
//	//
//
//	// Do a sanity check on the two strings.
//	if(file==NULL){
//		file="[Unspecified File]";
//	}
//	if(msg==NULL){
//		msg="[No error msg given.]";
//	}
//
//	// Send panic msg to stderr.  
//        fprintf(stderr, "%s::%d  PANIC:  program aborting. '%s'\n",
//                file, line, msg);
//
//	// If we have an error file, send the message there as well.
//	//   (Not guaranteed to work, which is why we send to stderr first.)
//	if(fperrors!=NULL){
//		fprintf(fperrors, "%s::%d  PANIC:  program aborting. '%s'\n",
//			file, line, msg);
//	}
//
//	// Cleanup code should go here.
//	fflush(NULL);
//
//	// Bye!
//	abort();
//}

// TODO: remove, keep track of comparisons only for debugging
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


void delete_time_course(TimeCourse *start2)
{
  TimeCourse *info, *start;
  
  start = start2; 
   
  while (start){
    info = start;
    start = start->next;
    free(info);
  }
}

//void display2(TimeCourse *start)
//{
//  TimeCourse *info;
//
//  info = start;
//  while (info){
//    fprintf(fperrors, "time %g conc %g\n", info->time, info->concentration);
//    info = info->next;
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

void sls_store_end2(TimeCourse *i, 
                    TimeCourse **start, 
                    TimeCourse **last)
{
  i->next = NULL;
  if (!*last) *start = i;
  else (*last)->next = i;
  *last = i;
}

//void display(FixedEvent *start)
//{
//  FixedEvent *info;
//
//  info = start;
//  while (info){
//    fprintf(fperrors,"gene %d time %f\n",info->event_id,info->time);
//    info = info->next;
//  }
//}