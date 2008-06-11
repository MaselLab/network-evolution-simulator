all: netsim

CC = gcc
CFLAGS = -g -std=c99 $(EXTRACFLAGS)
#  -O3
OBJS = random.o lib.o
LIBS = -lm 
#-lefence
OTHER = Makefile

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-check: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-gprof: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg $<

netsim: $(OBJS)
netsim-check: $(OBJS)
netsim-gprof: $(OBJS)

## check specific directory
check:	netsim-check
	./netsim-check -r 4 -d output
	@diff -r --exclude=.svn --exclude=NOTES output regression-tests/after-kon-change-r-4 && echo -e "************\nPassed regression\n***********"

profiling:	netsim-gprof

clean:
	rm -f netsim netsim-check netsim-gprof $(OBJS)
