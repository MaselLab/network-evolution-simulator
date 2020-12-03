#compiling option
CFLAGS = -O3 # the default is performance release
CFLAGS_INTEL = -fp-model source -fp-model precise # enable accurate arithematics on floating numbers in intel C compiler

#other flags go here
#CPPFLAGS =   

#objects
objects=main.o numerical.o netsim.o lib.o RngStream.o mutation.o cellular_activity.o

%.o: %.c
ifeq ($(CC),icc) 
	$(CC) $(CFLAGS) $(CFLAGS_INTEL) $(CPPFLAGS) -c $< -o $@ -qopenmp
else
	$(CC) $(CFLAGS) $(CPPFLAGS) -c $< -o $@ -fopenmp
endif

#target
simulator: $(objects)	
ifeq ($(CC),icc)
	$(CC) $(CFLAGS) $(CFLAGS_INTEL) $(CPPFLAGS) -o simulator $(objects) -qopenmp
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

#clean
.PHONY:clean

clean:
	-rm simulator $(objects)
