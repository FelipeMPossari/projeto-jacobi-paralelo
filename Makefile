SEQ_EXEC = jacobi_sequencial
PAR_EXEC = jacobi_paralelo

CC = gcc

CFLAGS_SEQ = -Wall

CFLAGS_PAR = -Wall -fopenmp

SRCDIR = src

all: $(SEQ_EXEC) $(PAR_EXEC)

$(SEQ_EXEC): $(SRCDIR)/$(SEQ_EXEC).c
	$(CC) $(CFLAGS_SEQ) $< -o $@

$(PAR_EXEC): $(SRCDIR)/$(PAR_EXEC).c
	$(CC) $(CFLAGS_PAR) $< -o $@

clean:
	rm -f $(SEQ_EXEC) $(PAR_EXEC)