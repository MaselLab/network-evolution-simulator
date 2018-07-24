/* 
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents (University of Arizona)
 */
#ifndef LIB_H
#define LIB_H

#include <stdio.h>
#include "netsim.h" 
#include "cellular_activity.h"
#include "netsim.h"

/* 
 * define macros for logging warning/errors 
 */
#ifndef LOG_OFF
#define LOG(...) { FILE *fperror; fperror=fopen("error.txt","a+"); fprintf(fperror, "%s: ", __func__); fprintf (fperror, __VA_ARGS__) ; fflush(fperror); fclose(fperror);} 
#endif

int add_fixed_event(int, float, FixedEvent **, FixedEvent **);

void delete_fixed_event(int, int, FixedEvent **, FixedEvent **);

void delete_fixed_event_from_head(FixedEvent **, FixedEvent **);

int check_concurrence(CellState *, float);

void free_fixedevent(CellState *);

void release_memory(Genotype*, Genotype *, RngStream *, RngStream[N_THREADS]);

#endif