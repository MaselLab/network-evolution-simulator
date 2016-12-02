/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include "netsim.h" 


extern void delete_queues(CellState *state);
extern void free_mem_CellState(CellState *state);
extern int sls_store(FixedEvent *i, 
		     FixedEvent **start, 
		      FixedEvent **last);
extern void delete_time_course(TimeCourse *start2);
extern void display(FixedEvent *start);
extern void display2(TimeCourse *start);
extern void sls_store_end(FixedEvent *i, 
			  FixedEvent **start, 
			  FixedEvent **last);
extern void sls_store_end2(TimeCourse *i, 
			   TimeCourse **start, 
			   TimeCourse **last);
extern void remove_from_array(int,
                              int,
                              int [],
                              int *,
                              int );
extern void create_output_directory(char *);
extern void create_output_file(char [80], char *, FILE **, int);
extern void read_kdisassembly(float [NUM_K_DISASSEMBLY]);
