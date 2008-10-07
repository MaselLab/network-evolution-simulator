all: netsim

CC = gcc
CFLAGS = -g -std=c99 $(EXTRACFLAGS)
#  -O3
OBJS = random.o lib.o netsim.o priority-queue.o
LIBS = -lm
#-lefence
OTHER = Makefile

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

%: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

#%-check: %.c %.h $(OBJS) $(OTHER)
#	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) $<

%-check: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main.c

main-matrix: main-matrix.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main-matrix.c

%-selection: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main.c

%-gprof: %.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg $<

qtest: qtest.c $(OBJS) priority-queue.o 
	$(CC) $(CFLAGS) $(OBJS) priority-queue.o -o $@ $(LIBS) qtest.c

netsim: $(OBJS)
netsim-check: $(OBJS)
netsim-selection: $(OBJS)
netsim-gprof: $(OBJS)
netsim-bigtf: $(OBJS)

## common command for doing regression test diff
DIFF_CMD := @diff -r  --exclude=tfsbound.dat --exclude=.svn --exclude=NOTES --exclude=cellsize.dat --exclude=growthrate.dat --exclude=netsimerrors.txt RUN regression-tests/ORIG && echo -e "************\nPassed regression\n***********"

## check specific directory
check-haploid:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DNO_SEPARATE_GENE" netsim-check
	./netsim-check -r 4 -p 1 -d output -c -1.0
	$(subst RUN,output,$(subst ORIG,2008-08-29-haploid-dilution-hind-15-r-4,$(DIFF_CMD)))

check-diploid:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15  -DNO_SEPARATE_GENE" netsim-check
	./netsim-check -r 4 -p 2 -d output -c -1.0
	$(subst RUN,output,$(subst ORIG,2008-08-29-diploid-dilution-hind-15-r-4,$(DIFF_CMD)))

check-replication:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15  -DNO_SEPARATE_GENE" netsim-check
	./netsim-check -r 4 -p 2 -d output -c 0.55
	$(subst RUN,output,$(subst ORIG,2008-08-29-replication-dilution-hind-15-r-4,$(DIFF_CMD)))

check-selection:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15" netsim-selection
	./netsim-selection -r 4 -p 2 -d selection -c -1.0
	$(subst RUN,selection,$(subst ORIG,2008-09-25-selection-r-4,$(DIFF_CMD)))

check-sample-output:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15" netsim-selection
	./netsim-selection -r 4 -p 2 -d selection -c -1.0 -t 150
	$(subst RUN,selection,$(subst ORIG,2008-08-29-sample-output-11genes-diploid-r-4,$(DIFF_CMD)))

profiling:	netsim-gprof

clean:
	rm -f netsim netsim-check netsim-gprof $(OBJS)
