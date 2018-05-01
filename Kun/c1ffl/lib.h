/* -*- Mode: C; tab-width: 2; indent-tabs-mode: nil; c-basic-offset: 2 -*- */
/* 
 * Yeast transcriptional network simulator
 *
 * Authors: Joanna Masel, Alex Lancaster
 * Copyright (c) 2007, 2008, 2009 Arizona Board of Regents (University of Arizona)
 */
#ifndef LIB_H
#define LIB_H

#include <stdio.h>
#include "netsim.h" 
#include "cellular_activity.h"
#include "netsim.h"

int add_fixed_event(int, float, FixedEvent **, FixedEvent **);

void delete_fixed_event(int, int, FixedEvent **, FixedEvent **);

void delete_fixed_event_from_head(FixedEvent **, FixedEvent **);

int check_concurrence(CellState *, float);

void free_fixedevent(CellState *);

void release_memory(Genotype*, Genotype *, RngStream *, RngStream[N_THREADS]);

#endif