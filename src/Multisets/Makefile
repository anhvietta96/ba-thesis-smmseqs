CC?=gcc
CFLAGS=-g -m64 -Wall -Werror -Wunused-parameter -Wunused-variable -O3

all:multisets.o

%.o: %.c
	$(CC) -c $< -o $@ $(CFLAGS) -MT $@ -MMD -MP -MF $(@:.o=.d)
