all: netsim

%: %.c
	gcc34 -g -std=c99 -o $@ -lm $<