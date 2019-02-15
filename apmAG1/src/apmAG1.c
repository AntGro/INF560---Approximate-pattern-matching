/**
 * APPROXIMATE PATTERN MATCHING
 *
 * INF560
 */
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/time.h>
#include <mpi.h>

#define APM_DEBUG 0

char * read_input_file(char *filename, int *size) {
    char *buf;
    off_t fsize;
    int fd = 0;
    int n_bytes = 1;

    /* Open the text file */
    fd = open(filename, O_RDONLY);
    if (fd == -1) {
        fprintf(stderr, "Unable to open the text file <%s>\n", filename);
        return NULL;
    }


    /* Get the number of characters in the textfile */
    fsize = lseek(fd, 0, SEEK_END);
    if (fsize == -1) {
        fprintf(stderr, "Unable to lseek to the end\n");
        return NULL;
    }

#if APM_DEBUG
    printf( "File length: %lld\n", fsize ) ;
#endif

    /* Go back to the beginning of the input file */
    if (lseek(fd, 0, SEEK_SET) == -1) {
        fprintf(stderr, "Unable to lseek to start\n");
        return NULL;
    }

    /* Allocate data to copy the target text */
    buf = (char *) malloc(fsize * sizeof(char));
    if (buf == NULL) {
        fprintf(stderr, "Unable to allocate %lld byte(s) for main array\n",
                fsize);
        return NULL;
    }

    n_bytes = read(fd, buf, fsize);
    if (n_bytes != fsize) {
        fprintf(stderr,
                "Unable to copy %lld byte(s) from text file (%d byte(s) copied)\n",
                fsize, n_bytes);
        return NULL;
    }

#if APM_DEBUG
    printf( "Number of read bytes: %d\n", n_bytes ) ;
#endif

    *size = n_bytes;


    close(fd);


    return buf;
}


#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))
#define max(a,b) (a>=b?a:b)
#define min(a,b) (a<=b?a:b)

int levenshtein(char *s1, char *s2, int len, int *column) {
    unsigned int x, y, lastdiag, olddiag;

    for (y = 1; y <= len; y++) {
        column[y] = y;
    }
    for (x = 1; x <= len; x++) {
        column[0] = x;
        lastdiag = x - 1;
        for (y = 1; y <= len; y++) {
            olddiag = column[y];
            column[y] = MIN3(
                    column[y] + 1,
                    column[y - 1] + 1,
                    lastdiag + (s1[y - 1] == s2[x - 1] ? 0 : 1)
            );
            lastdiag = olddiag;

        }
    }
    return (column[len]);
}


int main(int argc, char **argv) {
    char **pattern;
    char *filename;
    int approx_factor = 0;
    int nb_patterns = 0;
    int i, j;
    char *buf;
    struct timeval t1, t2;
    double duration;
    int n_bytes;
    int *n_matches;

    /* Check number of arguments */
    if (argc < 4) {
        printf("Usage: %s approximation_factor "
               "dna_database pattern1 pattern2 ...\n",
               argv[0]);
        return 1;
    }

    MPI_Init(&argc, &argv);

    int n, rk;

    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);
    /* Get the distance factor */
    approx_factor = atoi(argv[1]);

    /* Grab the filename containing the target text */
    filename = argv[2];

    /* Get the number of patterns that the user wants to search for */
    nb_patterns = argc - 3;

    /* Fill the pattern array */
    pattern = (char **) malloc(nb_patterns * sizeof(char *));
    if (pattern == NULL) {
        fprintf(stderr,
                "Unable to allocate array of pattern of size %d\n",
                nb_patterns);
        return 1;
    }

    /* len of the longest pattern */
    int maxLen = 0;

    /* Grab the patterns */
    for (i = 0; i < nb_patterns; i++) {
        int l;

        l = strlen(argv[i + 3]);
        maxLen = max(maxLen, l);
        if (l <= 0) {
            fprintf(stderr, "Error while parsing argument %d\n", i + 3);
            return 1;
        }

        pattern[i] = (char *) malloc((l + 1) * sizeof(char));
        if (pattern[i] == NULL) {
            fprintf(stderr, "Unable to allocate string of size %d\n", l);
            return 1;
        }

        strncpy(pattern[i], argv[i + 3], (l + 1));
    }

    /* Timer start */
    gettimeofday(&t1, NULL);


    /* Dispatch the data */
    int *sendcounts;
    int *displs;
    char *bufRec;
    int recvcount;
    if (rk == 0) {
        printf("Approximate Pattern Mathing: "
               "looking for %d pattern(s) in file %s w/ distance of %d\n",
               nb_patterns, filename, approx_factor);

        buf = read_input_file(filename, &n_bytes);
        if (buf == NULL) {
            MPI_Abort(MPI_COMM_WORLD, 0);
            return 1;
        }
        MPI_Bcast(&n_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Bcast(&n_bytes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    }

    /* Create communicator */

    int color = (rk < n_bytes);
    int newN, newRk;


    MPI_Comm MPI_COMM_NEW;
    MPI_Comm_split(MPI_COMM_WORLD, color, rk, &MPI_COMM_NEW);
    MPI_Comm_size(MPI_COMM_NEW, &newN);
    MPI_Comm_size(MPI_COMM_NEW, &newRk);

    if (rk == 0) {
        sendcounts = (int *) malloc(sizeof(int)*newN);
        displs = (int *) malloc(sizeof(int)*newN);
        int commonSize = n_bytes / newN;
        int aux = 0;
        for (i = 0; i < newN; i++) {
            /* Compute what is to be sent */
            sendcounts[i] = commonSize + min(maxLen, n_bytes - i * commonSizemin);
            displs[i] = i * commonSize
        }
`       sendcounts[newN - 1] = n_bytes - ((newN - 1) * commonSize);
        MPI_Scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_NEW);
        bufRec = (char *) malloc(sizeof(char)*recvcount);
        MPI_Scatterv(buf, sendcounts, displs, MPI_CHAR, bufRec, recvcount, MPI_CHAR, 0, MPI_COMM_NEW);
    }
        
    else {
        if (rk < n_bytes) {
            MPI_Scatter(sendcounts, 1, MPI_INT, &recvcount, 1, MPI_INT, 0, MPI_COMM_NEW);
	        bufRec = (char *) malloc(sizeof(char)*recvcount);
            MPI_Scatterv(buf, sendcounts, displs, MPI_CHAR, bufRec, recvcount, MPI_CHAR, 0, MPI_COMM_NEW);
        }
    }

    /* Allocate the array of matches */
    int *n_matches_loc = (int *) malloc(nb_patterns * sizeof(int));
    if (n_matches_loc == NULL) {
        fprintf(stderr, "Error: unable to allocate memory for %ldB\n",
                nb_patterns * sizeof(int));
        return 1;
    }
    
    /*****
     * BEGIN MAIN LOOP
     ******/
    if (rk < n_bytes) {
        n_bytes = recvcount;

    /* Check each pattern one by one */
    for (i = 0; i < nb_patterns; i++) {
        int size_pattern = strlen(pattern[i]);
        int *column;

        /* Initialize the number of matches to 0 */
        n_matches_loc[i] = 0;

        column = (int *) malloc((size_pattern + 1) * sizeof(int));
        if (column == NULL) {
            fprintf(stderr, "Error: unable to allocate memory for column (%ldB)\n",
                    (size_pattern + 1) * sizeof(int));
            return 1;
        }

        /* Traverse the input data up to the end of the file */
        for (j = 0; j < ((rk != newN -1) ? n_bytes - maxLen: n_bytes); j++) {
            int distance = 0;
            int size;

#if APM_DEBUG
            if ( j % 100 == 0 )
            {
            printf( "Procesing byte %d (out of %d)\n", j, n_bytes ) ;
            }
#endif

            size = size_pattern;
            if (n_bytes - j < size_pattern) {
                size = n_bytes - j;
            }

            distance = levenshtein(pattern[i], &bufRec[j], size, column);

            if (distance <= approx_factor) {
                n_matches_loc[i]++;
            }
        }
        free(column);
    }

        /*****
        * END MAIN LOOP
        ******/
        if (rk == 0) {
            n_matches = (int *) malloc(nb_patterns * sizeof(int));
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_NEW);
            /* Timer stop */
            gettimeofday(&t2, NULL);

            duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

            printf("APM done in %lf s\n", duration);

            for (i = 0; i < nb_patterns; i++) {
                printf("Number of matches for pattern <%s>: %d\n",
                    pattern[i], n_matches[i]);
            }
        }
        else {
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_NEW);
        }
    }
    

    return 0;
}
