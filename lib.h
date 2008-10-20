#include "netsim.h" 

extern float rtsafe(void (*funcd)(float, float, GillespieRates *, KonStates *, float *, float *), 
		    float x, GillespieRates *rates, KonStates *konStates, float x1, float x2, float xacc);
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
