# Define os nomes dos executáveis
SEQ_EXEC = jacobi_sequencial
PAR_EXEC = jacobi_paralelo

# Define o compilador
CC = gcc

# Flags para a versão sequencial
CFLAGS_SEQ = -Wall

# Flags para a versão paralela (inclui -fopenmp)
CFLAGS_PAR = -Wall -fopenmp

# Diretório dos arquivos fonte
SRCDIR = src

all: $(SEQ_EXEC) $(PAR_EXEC)

# Regra para compilar a versão sequencial
$(SEQ_EXEC): $(SRCDIR)/$(SEQ_EXEC).c
	$(CC) $(CFLAGS_SEQ) $< -o $@

# Regra para compilar a versão paralela
$(PAR_EXEC): $(SRCDIR)/$(PAR_EXEC).c
	$(CC) $(CFLAGS_PAR) $< -o $@

# Regra para limpar os executáveis
clean:
	rm -f $(SEQ_EXEC) $(PAR_EXEC)