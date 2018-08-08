/* 
 * Authors: Joanna Masel, Alex Lancaster, Kun Xiong
 * Copyright (c) 2018 Arizona Board of Regents on behalf of the University of Arizona
 
 * This file is part of network-evolution-simulator.
 * network-evolution-simulator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * network-evolution-simulator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 * You should have received a copy of the GNU Affero General Public License
 * along with network-evolution-simulator. If not, see <https://www.gnu.org/licenses/>.
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
