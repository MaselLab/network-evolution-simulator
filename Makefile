all: netsim

%: %.c
	gcc34 -std=c99 -o $@ -lm $<