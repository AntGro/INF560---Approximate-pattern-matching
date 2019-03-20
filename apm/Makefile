SRC_DIR=src
HEADER_DIR=include
OBJ_DIR=obj
OBJ_CUDA_DIR=obj_cuda

CC=mpicc
CFLAGS=-O3 -I$(HEADER_DIR) -Wall -fopenmp
LDFLAGS=
CUDA_LINK=-L/usr/local/cuda/bin/..//lib64 -lcudadevrt -lcudart_static -lrt -lpthread -ldl

SRC= apm.c

OBJ= $(OBJ_DIR)/apm.o $(OBJ_CUDA_DIR)/cuda.o

all: $(OBJ_DIR) $(OBJ_CUDA_DIR) apm

$(OBJ_DIR):
	mkdir $(OBJ_DIR)

$(OBJ_CUDA_DIR):
	mkdir $(OBJ_CUDA_DIR)
    
$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) $(CFLAGS) -c -o $@ $^

cuda_run: $(SRC_DIR)/cuda.cu
	nvcc -I. -o cuda $(SRC_DIR)/cuda.cu

$(OBJ_CUDA_DIR)/%.o : $(SRC_DIR)/%.cu
	nvcc -I. -c -o $@ $^

cuda:
	nvcc -I. -c -o $(OBJ_DIR)/cuda.o $(SRC_DIR)/cuda.cu

apm:$(OBJ)
	mpic++ $(CFLAGS) $(LDFLAGS) -o $@ $^ $(CUDA_LINK)

clean:
	rm -f apm $(OBJ) ; rmdir $(OBJ_DIR) ; rmdir $(OBJ_CUDA_DIR)
