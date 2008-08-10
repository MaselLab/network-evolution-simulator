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

%-selection: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-gprof: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg $<

netsim: $(OBJS)
netsim-check: $(OBJS)
netsim-selection: $(OBJS)
netsim-gprof: $(OBJS)
netsim-bigtf: $(OBJS)


## check specific directory
check-haploid:	clean
	make EXTRACFLAGS="-m32 -DPLOIDY=1 -DHIND_LENGTH=6" netsim-check
	./netsim-check -r 4 -d output
	@diff -r --exclude=.svn --exclude=NOTES --exclude=cellsize.dat --exclude=growthrate.dat --exclude=netsimerrors.txt output regression-tests/after-remove-lopt-haploid-r-4 && echo -e "************\nPassed regression\n***********"

check-diploid:	clean
	make EXTRACFLAGS="-m32 -DPLOIDY=2 -DHIND_LENGTH=6" netsim-check
	./netsim-check -r 4 -d output
	@diff -r --exclude=.svn --exclude=NOTES --exclude=cellsize.dat --exclude=growthrate.dat --exclude=netsimerrors.txt output regression-tests/after-remove-lopt-diploid-r-4 && echo -e "************\nPassed regression\n***********"

check-selection:	clean
	make EXTRACFLAGS="-m32 -DPLOIDY=2 -DHIND_LENGTH=15 -DSELECTION_GENE=10" netsim-selection
	./netsim-selection -r 4 -d selection

profiling:	netsim-gprof

clean:
	rm -f netsim netsim-check netsim-gprof $(OBJS)
