SRC_DIR=src
HEADER_DIR=include
OBJ_DIR=obj

CC=mpicc
CFLAGS=-O3 -I$(HEADER_DIR) -Wall -fopenmp
LDFLAGS=

SRC= apmAG1.c

OBJ= $(OBJ_DIR)/apmAG1.o

all: $(OBJ_DIR) apmAG1

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $^

apmAG1:$(OBJ)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $^

clean:
	rm -f apmAG1 $(OBJ) ; rmdir $(OBJ_DIR)
