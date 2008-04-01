all: netsim

%: %.c
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -g -std=c99 -o $@ -lm $<

%-gprof: %.c
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -g -std=c99 -o $@ -lm -pg $<


check:	netsim
	./netsim
	@diff -r --exclude=.svn --exclude=NOTES output test && echo -e "************\nPassed regression\n***********"

profiling:	netsim-gprof

clean:
	rm netsim