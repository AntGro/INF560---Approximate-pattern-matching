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
#include <omp.h>

#define APM_DEBUG 0

char *
read_input_file(char *filename, int *size) {
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
#define max(a, b) (a>=b?a:b)
#define min(a, b) (a<=b?a:b)

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


double mpi_data_split(int argc, char **argv, int n, int rk) {
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
        return 1.0;
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
            return 1.0;
        }

        pattern[i] = (char *) malloc((l + 1) * sizeof(char));
        if (pattern[i] == NULL) {
            fprintf(stderr, "Unable to allocate string of size %d\n", l);
            return 1.0;
        }

        strncpy(pattern[i], argv[i + 3], (l + 1));
    }

    /* Timer start */
    gettimeofday(&t1, NULL);

    if (rk == 0) {
        printf("Approximate Pattern Mathing: "
               "looking for %d pattern(s) in file %s w/ distance of %d (function called: mpi_data_split)\n",
               nb_patterns, filename, approx_factor);
    }

    /* Every process loads the data */
    buf = read_input_file(filename, &n_bytes);
    if (buf == NULL) {
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 1.0;
    }

    int rkMax = min(n, n_bytes) - 1;
    int start_check = rk * max (n_bytes / n, 1);
    int end_check = (rk + 1) * max (n_bytes / n, 1);

    if (rk > rkMax) {
        start_check = 0;
        end_check = 0;
    }

    if (rk == rkMax) {
        end_check = n_bytes;
    }
    /* Create communicator */
    /*
    int color = (rk < n_bytes);
    int newN, newRk;


    MPI_Comm MPI_COMM_NEW;
    MPI_Comm_split(MPI_COMM_WORLD, color, rk, &MPI_COMM_NEW);
    MPI_Comm_size(MPI_COMM_NEW, &newN);
    MPI_Comm_size(MPI_COMM_NEW, &newRk);
*/
    int *n_matches_loc = (int *) malloc((nb_patterns) * sizeof(int));
    /*****
     * BEGIN MAIN LOOP
     ******/
    if (rk <= rkMax) {
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
                return 1.0;
            }

            /* Traverse the input data up to the end of the file */
            for (j = start_check; j < end_check; j++) {
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

                distance = levenshtein(pattern[i], &buf[j], size, column);

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
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            /* Timer stop */
            gettimeofday(&t2, NULL);

            duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

            printf("APM done in %lf s\n", duration);

            for (i = 0; i < nb_patterns; i++) {
                printf("Number of matches for pattern <%s>: %d\n",
                       pattern[i], n_matches[i]);
            }
        } else {
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }
    }
    return duration;
}

double mpi_omp_data_split(int argc, char **argv, int n, int rk) {
    char **pattern;
    char *filename;
    int approx_factor = 0;
    int nb_patterns = 0;
    int i;
    char *buf;
    struct timeval t1, t2;
    double duration;
    int n_bytes;
    int *n_matches;

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
        return 1.0;
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
            return 1.0;
        }

        pattern[i] = (char *) malloc((l + 1) * sizeof(char));
        if (pattern[i] == NULL) {
            fprintf(stderr, "Unable to allocate string of size %d\n", l);
            return 1.0;
        }

        strncpy(pattern[i], argv[i + 3], (l + 1));
    }

    /* Timer start */
    gettimeofday(&t1, NULL);

    if (rk == 0) {
        printf("Approximate Pattern Mathing: "
               "looking for %d pattern(s) in file %s w/ distance of %d (function called: mpi_omp_data_split)\n",
               nb_patterns, filename, approx_factor);
    }

    /* Every process loads the data */
    buf = read_input_file(filename, &n_bytes);
    if (buf == NULL) {
        MPI_Abort(MPI_COMM_WORLD, 0);
        return 1.0;
    }

    int rkMax = min(n, n_bytes) - 1;
    int start_check = rk * max (n_bytes / n, 1);
    int end_check = (rk + 1) * max (n_bytes / n, 1);

    if (rk > rkMax) {
        start_check = 0;
        end_check = 0;
    }

    if (rk == rkMax) {
        end_check = n_bytes;
    }

    int *n_matches_loc = (int *) malloc((nb_patterns) * sizeof(int));
    /*****
     * BEGIN MAIN LOOP
     ******/
    if (rk <= rkMax) {
        /* Check each pattern one by one */
#pragma omp parallel default(shared)
        {
#pragma omp for schedule (static)
            for (i = 0; i < nb_patterns; i++) {
                int size_pattern = strlen(pattern[i]);
                int *column;
                int j;

                /* Initialize the number of matches to 0 */
                n_matches_loc[i] = 0;

                column = (int *) malloc((size_pattern + 1) * sizeof(int));
                if (column == NULL) {
                    fprintf(stderr, "Error: unable to allocate memory for column (%ldB)\n",
                            (size_pattern + 1) * sizeof(int));
                    //return 1; <- need to handle that!
                }

                /* Traverse the input data up to the end of the file */
                for (j = start_check; j < end_check; j++) {
                    int distance = 0;
                    int size = size_pattern;

#if APM_DEBUG
                    if ( j % 100 == 0 )
            {
            printf( "Procesing byte %d (out of %d)\n", j, n_bytes ) ;
            }
#endif

                    if (n_bytes - j < size_pattern) {
                        size = n_bytes - j;
                    }

                    distance = levenshtein(pattern[i], &buf[j], size, column);

                    if (distance <= approx_factor) {
                        n_matches_loc[i]++;
                    }
                }
                //printf ("Pattern %s handled by thread %d/%d of the #%d process out of %d: %d\n", pattern[i], omp_get_thread_num(), omp_get_num_threads(), rk, n, n_matches_loc[i]);
                free(column);
            }
        }


        /*****
        * END MAIN LOOP
        ******/
        if (rk == 0) {
            n_matches = (int *) malloc(nb_patterns * sizeof(int));
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
            /* Timer stop */
            gettimeofday(&t2, NULL);

            duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

            printf("APM done in %lf s\n", duration);

            for (i = 0; i < nb_patterns; i++) {
                printf("Number of matches for pattern <%s>: %d\n",
                       pattern[i], n_matches[i]);
            }
        } else {
            MPI_Reduce(n_matches_loc, n_matches, nb_patterns, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        }


    }
    return duration;
}

double mpi_pattern_split(int argc, char **argv) {

    int my_rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &my_rank);
    int comm_size;
    MPI_Comm_size (MPI_COMM_WORLD, &comm_size);

    char ** pattern ;
    char * filename ;
    int approx_factor = 0 ;
    int nb_patterns = 0 ;
    int i, j ;
    char * buf ;
    struct timeval t1, t2;
    double duration ;
    int n_bytes ;
    int * n_matches_rank ;

    /* Check number of arguments */
    if ( argc < 4 )
    {
        printf( "Usage: %s approximation_factor "
                "dna_database pattern1 pattern2 ...\n",
                argv[0] ) ;
        return 1 ;
    }

    /* Get the distance factor */
    approx_factor = atoi( argv[1] ) ;

    /* Grab the filename containing the target text */
    filename = argv[2] ;

    /* Get the number of patterns that the user wants to search for */
    nb_patterns = argc - 3 ;
    int modulo = nb_patterns%comm_size;
    int countData = nb_patterns/comm_size;
    int nb_patterns_rank;

    /* Fill the pattern array */

    if (my_rank>=modulo) {
        pattern = (char **)malloc( countData * sizeof( char * ) ) ;
        nb_patterns_rank = countData;
        if ( pattern == NULL )
        {
            fprintf( stderr,
                     "Unable to allocate array of pattern of size %d\n",
                     countData ) ;
            return 1 ;
        }
    }
    else {
        pattern = (char **)malloc( countData+1 * sizeof( char * ) ) ;
        nb_patterns_rank = countData + 1;
        if ( pattern == NULL )
        {
            fprintf( stderr,
                     "Unable to allocate array of pattern of size %d\n",
                     countData+1 ) ;
            return 1 ;
        }
    }


    /* Grab the patterns */
    for ( i = 0 ; i < nb_patterns ; i++ )
    {

        if (i%comm_size == my_rank){

            int l ;
            l = strlen(argv[i+3]) ;
            if ( l <= 0 )
            {
                fprintf( stderr, "Error while parsing argument %d\n", i+3 ) ;
                return 1 ;
            }

            pattern[i/comm_size] = (char *)malloc( (l+1) * sizeof( char ) ) ;
            if ( pattern[i/comm_size] == NULL )
            {
                fprintf( stderr, "Unable to allocate string of size %d\n", l ) ;
                return 1 ;
            }

            strncpy( pattern[i/comm_size], argv[i+3], (l+1) ) ;
        }
    }

    gettimeofday(&t1, NULL);

    buf = read_input_file( filename, &n_bytes ) ;
    if ( buf == NULL )
    {
        return 1 ;
    }

    /* Allocate the array of matches */
    if (my_rank>=modulo) {
        n_matches_rank = (int *)malloc( countData * sizeof( int ) ) ;
        if ( n_matches_rank == NULL )
        {
            fprintf( stderr, "Error: unable to allocate memory for %ldB\n",
                     countData * sizeof( int ) ) ;
            return 1 ;
        }
    }
    else {
        n_matches_rank = (int *)malloc( countData+1 * sizeof( int ) ) ;
        if ( n_matches_rank == NULL )
        {
            fprintf( stderr, "Error: unable to allocate memory for %ldB\n",
                     (countData+1) * sizeof( int ) ) ;
            return 1 ;
        }
    }



    /*****
     * BEGIN MAIN LOOP
     ******/

    /* Check each pattern one by one */
    for ( i = 0 ; i < nb_patterns_rank ; i++ )
    {
        int size_pattern = strlen(pattern[i]) ;
        int * column ;

        /* Initialize the number of matches to 0 */
        n_matches_rank[i] = 0 ;

        column = (int *)malloc( (size_pattern+1) * sizeof( int ) ) ;
        if ( column == NULL )
        {
            fprintf( stderr, "Error: unable to allocate memory for column (%ldB)\n",
                     (size_pattern+1) * sizeof( int ) ) ;
            return 1 ;
        }

        /* Traverse the input data up to the end of the file */
        for ( j = 0 ; j < n_bytes ; j++ )
        {
            int distance = 0 ;
            int size ;

#if APM_DEBUG
            if ( j % 100 == 0 )
          {
          printf( "Procesing byte %d (out of %d)\n", j, n_bytes ) ;
          }
#endif

            size = size_pattern ;
            if ( n_bytes - j < size_pattern )
            {
                size = n_bytes - j ;
            }

            distance = levenshtein( pattern[i], &buf[j], size, column ) ;

            if ( distance <= approx_factor ) {
                n_matches_rank[i]++ ;
            }
        }

        free( column );
    }


    if (my_rank>0){
        int j;
        for (j = 0 ; j < nb_patterns_rank ; j++){
            MPI_Request req;
            MPI_Status sta;
            MPI_Send(n_matches_rank+j, 1, MPI_INT, 0, j, MPI_COMM_WORLD);
        }

    }

    else{
        /* Timer start */
        printf("Approximate Pattern Mathing: "
               "looking for %d pattern(s) in file %s w/ distance of %d (function called: mpi_pattern_split)\n",
               nb_patterns, filename, approx_factor);
        char ** all_patterns;
        all_patterns = (char **)malloc( nb_patterns * sizeof( char * ) ) ;
        if ( all_patterns == NULL )
        {
            fprintf( stderr,
                     "Unable to allocate array of pattern of size %d\n",
                     nb_patterns ) ;
            return 1 ;
        }

        for ( i = 0 ; i < nb_patterns ; i++ )
        {
            int l ;
            l = strlen(argv[i+3]) ;
            if ( l <= 0 )
            {
                fprintf( stderr, "Error while parsing argument %d\n", i+3 ) ;
                return 1 ;
            }

            all_patterns[i] = (char *)malloc( (l+1) * sizeof( char ) ) ;
            if ( all_patterns[i] == NULL )
            {
                fprintf( stderr, "Unable to allocate string of size %d\n", l ) ;
                return 1 ;
            }

            strncpy( all_patterns[i], argv[i+3], (l+1) ) ;
        }

        int * n_matches ;
        n_matches = (int *)malloc( nb_patterns * sizeof( int ) );
        int sender;
        for (sender = 1 ; sender < comm_size ; sender++){
            int limit;
            if (sender >= modulo){
                limit = countData;
            }
            else {
                limit = countData + 1;
            }

            int j;
            for (j = 0 ; j < limit ; j++){
                MPI_Request req;
                MPI_Status sta;
                MPI_Recv(n_matches+j*comm_size+sender, 1, MPI_INT, sender, j, MPI_COMM_WORLD, &sta);
            }

        }

        for (j = 0 ; j < nb_patterns_rank ; j++){
            n_matches[comm_size*j] = n_matches_rank[j];
        }

        /* Timer stop */
        gettimeofday(&t2, NULL);
        duration = (t2.tv_sec -t1.tv_sec)+((t2.tv_usec-t1.tv_usec)/1e6);
        printf( "APM done in %lf s\n", duration ) ;

        for ( i = 0 ; i < nb_patterns ; i++ ) {
            printf( "Number of matches for pattern <%s>: %d\n",
                    all_patterns[i], n_matches[i] ) ;
        }

    }


    /*****
     * END MAIN LOOP
     ******/

    return duration ;
}

double sequential(int argc, char **argv, int *n_bytes_copy) {
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
        return 1.;
    }

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
        return 1.;
    }

    /* Grab the patterns */
    for (i = 0; i < nb_patterns; i++) {
        int l;

        l = strlen(argv[i + 3]);
        if (l <= 0) {
            fprintf(stderr, "Error while parsing argument %d\n", i + 3);
            return 1.;
        }

        pattern[i] = (char *) malloc((l + 1) * sizeof(char));
        if (pattern[i] == NULL) {
            fprintf(stderr, "Unable to allocate string of size %d\n", l);
            return 1.;
        }

        strncpy(pattern[i], argv[i + 3], (l + 1));
    }


    printf("Approximate Pattern Mathing: "
           "looking for %d pattern(s) in file %s w/ distance of %d (function called: sequential)\n",
           nb_patterns, filename, approx_factor);

    buf = read_input_file(filename, &n_bytes);
    *n_bytes_copy = n_bytes;
    if (buf == NULL) {
        return 1.;
    }

    /* Allocate the array of matches */
    n_matches = (int *) malloc(nb_patterns * sizeof(int));
    if (n_matches == NULL) {
        fprintf(stderr, "Error: unable to allocate memory for %ldB\n",
                nb_patterns * sizeof(int));
        return 1.;
    }

    /*****
     * BEGIN MAIN LOOP
     ******/

    /* Timer start */
    gettimeofday(&t1, NULL);

    /* Check each pattern one by one */
    for (i = 0; i < nb_patterns; i++) {
        int size_pattern = strlen(pattern[i]);
        int *column;

        /* Initialize the number of matches to 0 */
        n_matches[i] = 0;

        column = (int *) malloc((size_pattern + 1) * sizeof(int));
        if (column == NULL) {
            fprintf(stderr, "Error: unable to allocate memory for column (%ldB)\n",
                    (size_pattern + 1) * sizeof(int));
            return 1.;
        }

        /* Traverse the input data up to the end of the file */
        for (j = 0; j < n_bytes; j++) {
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

            distance = levenshtein(pattern[i], &buf[j], size, column);

            if (distance <= approx_factor) {
                n_matches[i]++;
            }
        }

        free(column);
    }

    /* Timer stop */
    gettimeofday(&t2, NULL);

    duration = (t2.tv_sec - t1.tv_sec) + ((t2.tv_usec - t1.tv_usec) / 1e6);

    printf("APM done in %lf s\n", duration);

    /*****
     * END MAIN LOOP
     ******/

    for (i = 0; i < nb_patterns; i++) {
        printf("Number of matches for pattern <%s>: %d\n",
               pattern[i], n_matches[i]);
    }

    return duration;
}

int get_N(int n, char* hostnames, int* lens, int* displs, int sum_lens) {
    char* hostnames_dico = (char*) malloc(sum_lens);
    int lens_dico[n], displs_dico[n];
    hostnames_dico = hostnames;
    lens_dico[0] = lens[0];
    displs_dico[0] = 0;
    int dico_size = 1;
    int i, j, k;
    for (i = 1; i < n; i++) {
        int new_name = 1;
        for (j = 0; j < dico_size; j++) {
            if (lens_dico[j] != lens[i]) continue;
            int same_name = 1;
            for (k = 0; k < lens[i]; k++) {
                if (hostnames_dico[displs_dico[j] + k] != hostnames[displs[i] + k]) {
                    same_name = 0;
                    break;
                }
            }
            if (same_name) {
                new_name = 0;
                break;
            }
        }
        if (new_name) {
            lens_dico[dico_size] = lens[i];
            displs_dico[dico_size] = displs_dico[dico_size-1] + lens_dico[dico_size-1];
            strncpy (hostnames_dico + displs_dico[dico_size], hostnames + displs[i], lens[i]);
            dico_size += 1;
        }
    }

    return dico_size;
}

int main(int argc, char **argv) {

    MPI_Init(&argc, &argv);

    /* Check number of arguments */
    if (argc < 4) {
        printf("Usage: %s approximation_factor "
               "dna_database pattern1 pattern2 ...\n",
               argv[0]);
        return 1;
    }

    int n, rk;

    MPI_Comm_size(MPI_COMM_WORLD, &n);
    MPI_Comm_rank(MPI_COMM_WORLD, &rk);

    //get number of Nodes
    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);

    int lens[n];
    int displs[n];
    int sum_lens;
    char* hostnames;
    if (rk==0) {
        MPI_Gather(&len, 1, MPI_INT, lens, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Reduce(&len, &sum_lens, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        hostnames = (char*) malloc(sum_lens);
        displs[0] = 0;
        int i;
        for (i = 1; i < n; i++) {
            displs[i] = displs[i-1] + lens[i-1];
        }
        MPI_Gatherv(name, len, MPI_CHAR, hostnames, lens, displs, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    else {
        MPI_Gather(&len, 1, MPI_INT, lens, 1, MPI_INT, 0, MPI_COMM_WORLD);
        MPI_Reduce(&len, &sum_lens, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Gatherv(name, len, MPI_CHAR, hostnames, lens, displs, MPI_CHAR, 0, MPI_COMM_WORLD);
    }

    int N;

    if (rk == 0) {
        // Get value of N
        N = get_N(n, hostnames, lens, displs, sum_lens);
    }

    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double duration_mpi_data_split = mpi_data_split(argc, argv, n, rk);

    double duration_mpi_omp_data_split = mpi_omp_data_split(argc, argv, n, rk);

    double duration_mpi_pattern_split = mpi_pattern_split(argc, argv);



    if (rk == 0) {
        int n_bytes;
        double duration_sequential = sequential(argc, argv, &n_bytes);

        char filename[200];
        sprintf(filename, "testsResults/%d %d %d %d.txt", n, N, n_bytes, argc - 3); //(n, N, size of dna, nb patterns)

        FILE *fp = fopen(filename, "wt");
        fprintf(fp, "sequential %f\n", duration_sequential);
        fprintf(fp, "mpi_data_split %f\n", duration_mpi_data_split);
        fprintf(fp, "mpi_omp_data_split %f\n", duration_mpi_omp_data_split);
        fprintf(fp, "mpi_pattern_split %f\n", duration_mpi_pattern_split);
        fclose(fp);
    }

    MPI_Finalize();

}