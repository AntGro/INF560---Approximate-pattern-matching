SRC_DIR=src
HEADER_DIR=include
OBJ_DIR=obj

CC=mpicc
CFLAGS=-O3 -I$(HEADER_DIR) -Wall
LDFLAGS=

SRC= apmMPI1.c

OBJ= $(OBJ_DIR)/apmMPI1.o

all: $(OBJ_DIR) apmMPI1

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $^

apmMPI1:$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f apmMPI1 $(OBJ) ; rmdir $(OBJ_DIR)
