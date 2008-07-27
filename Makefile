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

%: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-check: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-gprof: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg $<

netsim: $(OBJS)
netsim-check: $(OBJS)
netsim-gprof: $(OBJS)
netsim-bigtf: $(OBJS)

## check specific directory
check-haploid:	clean
	make EXTRACFLAGS="-m32 -DPLOIDY=1" netsim-check
	./netsim-check -r 4 -d output
	@diff -r --exclude=.svn --exclude=NOTES --exclude=cellsize.dat --exclude=growthrate.dat --exclude=netsimerrors.txt output regression-tests/after-kon-change-haploid-after-plus1-r-4 && echo -e "************\nPassed regression\n***********"

check-diploid:	clean
	make EXTRACFLAGS="-m32 -DPLOIDY=2" netsim-check
	./netsim-check -r 4 -d output
	@diff -r --exclude=.svn --exclude=NOTES --exclude=cellsize.dat --exclude=growthrate.dat --exclude=netsimerrors.txt output regression-tests/after-kon-change-diploid-after-plus1-r-4 && echo -e "************\nPassed regression\n***********"


profiling:	netsim-gprof

clean:
	rm -f netsim netsim-check netsim-gprof $(OBJS)
