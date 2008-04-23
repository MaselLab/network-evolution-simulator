all: netsim

%: %.c Makefile
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -O3 -g -std=c99 -o $@ -lm $<

%-check: %.c
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -g -std=c99 -o $@ -lm -pg $<

%-gprof: %.c
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -g -std=c99 -o $@ -lm -pg $<

check:	netsim-check
	./netsim-check
	@diff -r --exclude=.svn --exclude=NOTES output test && echo -e "************\nPassed regression\n***********"

profiling:	netsim-gprof

clean:
	rm netsim