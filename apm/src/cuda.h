#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>

int gpu_find_matches (int nb_patterns, char** pattern, char * buf, int n_bytes, int* n_matches, int approx_factor);

