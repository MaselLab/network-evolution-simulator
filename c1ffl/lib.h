/* 
 * Yeast transcriptional network simulator
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007-2018 Arizona Board of Regents (University of Arizona)
 */
#include <stdio.h>
#include "netsim.h" 

extern int sls_store(FixedEvent *i, 
		     FixedEvent **start, 
		      FixedEvent **last);
extern void delete_time_course(TimeCourse *start2);

extern void sls_store_end(FixedEvent *i, 
			  FixedEvent **start, 
			  FixedEvent **last);
extern void sls_store_end2(TimeCourse *i, 
			   TimeCourse **start, 
			   TimeCourse **last);
