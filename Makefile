all: netsim

%: %.c
#	gcc34 -g -std=c99 -o $@ -lm $<
	gcc -g -std=c99 -o $@ -lm $<

check:	netsim
	./netsim
	@diff -r --exclude=.svn --exclude=NOTES output check && echo -e "************\nPassed regression\n***********"

clean:
	rm netsim