all: netsim

CC = gcc
#Currently enabled for Linux needs to be enbled for Windows
#CFLAGS = -g -std=c99 $(EXTRACFLAGS) -Wall
CFLAGS = -g $(EXTRACFLAGS) -Wall
#  -O3 -pg
OBJS = random.o lib.o netsim.o priority-queue.o
LIBS = -lm
#-lefence
OTHER = Makefile

%.o: %.c %.h
	$(CC) $(CFLAGS) -c $<

netsim netsim-gs netsim-full-10 netsim-full-500: main.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main.c

%-check: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main.c

%-selection: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) main.c

%-gprof: %.c %.h $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) -pg main.c

main-matrix: main-matrix.c $(OBJS) $(OTHER)
	$(CC) $(CFLAGS) -I/usr/local/MathSoft/extern/include/ $(OBJS) -L/usr/local/MathSoft/bin/glnxa64  -Wl,-rpath-link /usr/local/MathSoft/bin/glnxa64 -o $@ $(LIBS) -leng -lmx -lut -lmat -lexpat -lmwfl main-matrix.c

qtest: qtest.c $(OBJS) priority-queue.o 
	$(CC) $(CFLAGS) $(OBJS) priority-queue.o -o $@ $(LIBS) qtest.c

rnd_test: rnd_test.c random.o
	$(CC) $(CFLAGS) random.o -o $@ $(LIBS) rnd_test.c

tf: tfonly.c $(OBJS) $(OTHER) 
	$(CC) $(CFLAGS) $(OBJS) -o $@ $(LIBS) tfonly.c

netsim: $(OBJS)
netsim-gs: $(OBJS)
netsim-full-500: $(OBJS)
netsim-full-10: $(OBJS)
netsim-check: $(OBJS)
netsim-selection: $(OBJS)
netsim-gprof: $(OBJS)
netsim-bigtf: $(OBJS)

main-matrix-32bit: clean
	make EXTRACFLAGS="-m32" main-matrix

main-matrix-64bit: clean
	make EXTRACFLAGS="-m64" main-matrix


tf-clean: clean
	make EXTRACFLAGS="-DTFGENES=10 -DNGENES=1 -DNPROTEINS=10" tf

## common command for doing regression test diff
DIFF_CMD := @diff -r --exclude=rounding*.dat --exclude=.svn --exclude=NOTES --exclude=koff*.dat  --exclude=tfsbound*.dat  --exclude=cellsize*.dat --exclude=growthrate*.dat --exclude=netsimerrors.txt RUN regression-tests/ORIG && echo -e "************\nPassed regression\n***********"

## check specific directory
check-haploid:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DNO_SEPARATE_GENE -DPOP_SIZE=1" netsim-check
	./netsim-check -r 4 -p 1 -d output -c -1.0
	$(subst RUN,output,$(subst ORIG,2009-05-26-haploid-with-ran1-r-4,$(DIFF_CMD)))

check-diploid:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15  -DNO_SEPARATE_GENE -DPOP_SIZE=1" netsim-check
	./netsim-check -r 4 -p 2 -d output -c -1.0
	$(subst RUN,output,$(subst ORIG,2009-05-26-diploid-with-ran1-r-4,$(DIFF_CMD)))

check-replication:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15  -DNO_SEPARATE_GENE -DPOP_SIZE=1" netsim-check
	./netsim-check -r 4 -p 2 -d output -c 0.505
	$(subst RUN,output,$(subst ORIG,2009-05-26-replication-with-ran1-r-4,$(DIFF_CMD)))

check-selection:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=1" netsim-selection
	./netsim-selection -r 4 -p 2 -d selection -c -1.0
	$(subst RUN,selection,$(subst ORIG,2009-05-26-selection-with-ran1-r-4,$(DIFF_CMD)))

check-burn-in:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=1" netsim-selection
	./netsim-selection -r 4 -p 2 -d selection -c -1.0 --kon 0.225 --konafter 1e-4
	$(subst RUN,selection,$(subst ORIG,2009-05-26-burn-in-diploid-r-4,$(DIFF_CMD)))


## TODO: not recently tested
check-sample-output:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=1" netsim-selection
	./netsim-selection -r 4 -p 2 -d selection -c -1.0 -t 150
	$(subst RUN,selection,$(subst ORIG,2008-08-29-sample-output-11genes-diploid-r-4,$(DIFF_CMD)))

## run for 4 divisions with reduced S and G2 phases
check-multiple-pops:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=4" netsim
	./netsim -r 5 -p 2 -d multiple-pops -c 0.505 -n -s 4 --timesphase 1.0 --timeg2phase 0.0
	$(subst RUN,multiple-pops,$(subst ORIG,2009-05-26-division-time-multiple-pops-r-5,$(DIFF_CMD)))

## all known working regression tests
check-all: check-haploid check-diploid check-replication check-burn-in check-multiple-pops check-selection 

run-full-pops: clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=500" netsim-full-500
	./netsim-full-500 -r 8 -p 2 -d multiple-pops-500 -c 1.0 -n -s 500 --timesphase 30.0 --timeg2phase 30.0 --growthscaling 10.0  --kon 0.2225 --konafter 1e-4

run-10-pops: clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=10" netsim-full-10
	./netsim-full-10 -r 3 -p 2 -d multiple-pops-10 -c 1.0 -n -s 60 --timesphase 30.0 --timeg2phase 30.0 --growthscaling 2.0  --kon 0.2225 --konafter 1e-4

run-gs:	clean
	make EXTRACFLAGS="-m32 -DHIND_LENGTH=15 -DPOP_SIZE=1" netsim-gs
	for ((gs=10;gs>=0;gs-=1)); do  \
	   g=$$(qalc -t "round(10^(2*$${gs}/10))"); \
	   for ((rand=04;rand<64;rand+=1)); do  \
	     echo "random number=" $${rand} "growth rate scaling=" $${g} ; \
	     dir=$$(printf "%03g" $${g})-$$(printf "%02d" $${rand}); \
	     echo $${dir} ; \
	     time ./netsim-gs -r $${rand} -p 2 -d growthscaling/$${dir} -c 1.0 -n -s 1 --timemax=1000.0 --timesphase 0.0 --timeg2phase 0.0 --growth $${g} --kon 0.2225 --konafter 1e-4 ; \
	   done \
	done	


profiling:	netsim-gprof

clean:
	rm -f netsim netsim-gs netsim-full-500 netsim-check netsim-gprof netsim-selection tfonly $(OBJS)
