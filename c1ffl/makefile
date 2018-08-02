#compiling option
CFLAGS = -O3 #the default is performance release 

#simulation mode
#CPPFLAGS #the default is full simulation.   

#objects
objects=main.o numerical.o netsim.o lib.o RngStream.o mutation.o cellular_activity.o

#target
simulator: $(objects)	
ifeq ($(CC),icc)
	$(CC) $(CFLAGS) $(CPPFLAGS) -o simulator $(objects) -qopenmp
else
	$(CC) $(CFLAGS) $(CPPFLAGS) -o simulator $(objects) -fopenmp -lm
endif

main.o: netsim.h RngStream.h lib.h	

numerical.o: netsim.h RngStream.h

RngStream.o: RngStream.h

mutation.o: netsim.h RngStream.h numerical.h mutation.h

cellular_activity.o: cellular_activity.h lib.h numerical.h

lib.o: lib.h

netsim.o: netsim.c
ifeq ($(CC),icc) 
	$(CC) $(CFLAGS) $(CPPFLAGS) -c netsim.c -qopenmp
else
	$(CC) $(CFLAGS) $(CPPFLAGS) -c netsim.c -fopenmp
endif

#clean
.PHONY:clean

clean:
	-rm simulator $(objects)
