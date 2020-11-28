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

#include "RngStream.h" 

/* Normally distributed random number with mean zero and variance 1 */
extern float gasdev(RngStream);

/* Returns an exponentially distributed, positive, random deviate of unit mean */
extern float expdev(RngStream);

/*Newton-Raphson root-finding method*/
extern float rtsafe(void (*funcd)(float, int, float, float*, float*, float*, float*, float*, int), 
		    int, float, float *, float*, float*, float, float, int);

/*find max or min*/
extern void find_max(float*, int, int, float*, int*, float);
extern void find_x(float*, int, int, float, float*, int);
