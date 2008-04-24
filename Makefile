all: netsim

CC = gcc
CFLAGS = -g -std=c99
OBJS = random.o lib.o
LIBS = -lm 
OTHER = Makefile

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%: %.c $(OBJS) $(OTHER)
	$(CC) -O3 $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-check: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-gprof: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg $<

netsim: $(OBJS)
netsim-check: $(OBJS)
netsim-gprof: $(OBJS)

check:	netsim-check
	./netsim-check
	@diff -r --exclude=.svn --exclude=NOTES output test && echo -e "************\nPassed regression\n***********"

profiling:	netsim-gprof

clean:
	rm -f netsim netsim-check netsim-gprof $(OBJS)