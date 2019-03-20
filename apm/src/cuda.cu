// CUDA runtime
#include <cuda_runtime.h>
#include <cuda.h>

#include "cuda.h"

int __device__ min3(int a, int b, int c) {
    return ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)));
}


int __device__ levenshtein_cuda(char *s1, char *s2, int len, int *column) {
    unsigned int x, y, lastdiag, olddiag;

    for (y = 1; y <= len; y++) {
        column[y] = y;
    }
    for (x = 1; x <= len; x++) {
        column[0] = x;
        lastdiag = x - 1;
        for (y = 1; y <= len; y++) {
            olddiag = column[y];
            column[y] = min3(
                    column[y] + 1,
                    column[y - 1] + 1,
                    lastdiag + (s1[y - 1] == s2[x - 1] ? 0 : 1)
            );
            lastdiag = olddiag;

        }
    }
    return (column[len]);
}

void __global__ matchesKernel(int* d_n_matches, char * d_buf, char * d_pattern, int i, int size_pattern, int offset, int n_bytes, int approx_factor){

    /* Traverse the input data up to the end of the file */
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;
    int distance = 0 ;
    int size ;

    size = size_pattern ;
    int* columns = (int *) malloc((size_pattern + 1) * sizeof(int));
    while (j < n_bytes) {
        if (n_bytes - j < size_pattern ){
            size = n_bytes - j ;
        }

        distance = levenshtein_cuda(d_pattern + offset, &d_buf[j], size, columns ) ;
        if ( distance <= approx_factor) {
            atomicAdd(&d_n_matches[i], 1);
        }

        j += stride;
    }
    free(columns);

}

int __host__ gpu_find_matches (int nb_patterns, char** pattern, char * buf, int n_bytes, int* n_matches, int approx_factor) {
/* Check each pattern one by one */
    int i;
    int* d_n_matches;
    char * d_pattern;
    char* d_buf;
    int* offset = (int *)malloc( nb_patterns * sizeof( int ) ) ;
    int* lens = (int *)malloc( nb_patterns * sizeof( int ) ) ;
    int sum_lens;
    lens[0] = strlen(pattern[0]);
    offset[0] = 0;
    sum_lens = lens[0];
    for (i = 1; i < nb_patterns; i++) {
        offset[i] = offset[i-1] + lens[i-1];
        lens[i] = strlen(pattern[i]);
        sum_lens += lens[i];
    }
    char* concat_patterns = (char*) malloc( sum_lens * sizeof( char ) ) ;
    for (i = 0; i < nb_patterns; i++) {
        strcpy (concat_patterns + offset[i], pattern[i]);
    }

    cudaMalloc((void **)&d_n_matches, nb_patterns*sizeof(int));
    cudaMalloc((void **)&d_pattern, sum_lens*sizeof(char));
    cudaMalloc((void **)&d_buf, n_bytes);
    cudaMemcpy(d_pattern, concat_patterns, sum_lens*sizeof(char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_buf, buf, n_bytes, cudaMemcpyHostToDevice);
    cudaMemcpy(d_n_matches, n_matches, nb_patterns*sizeof(int), cudaMemcpyHostToDevice);

    int Dg = 4;
    int Db = 256;
    for (i = 0; i < nb_patterns; i++) {
        matchesKernel<<<Dg,Db>>>(d_n_matches, d_buf, d_pattern, i, lens[i], offset[i], n_bytes, approx_factor);
    }

    cudaMemcpy(n_matches, d_n_matches, nb_patterns*sizeof(int), cudaMemcpyDeviceToHost);

    return 0;
}
