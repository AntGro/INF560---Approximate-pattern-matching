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

char * 
read_input_file( char * filename, int * size )
{
    char * buf ;
    off_t fsize;
    int fd = 0 ;
    int n_bytes = 1 ;

    /* Open the text file */
    fd = open( filename, O_RDONLY ) ;
    if ( fd == -1 ) 
    {
        fprintf( stderr, "Unable to open the text file <%s>\n", filename ) ;
        return NULL ;
    }


    /* Get the number of characters in the textfile */
    fsize = lseek(fd, 0, SEEK_END);
    if ( fsize == -1 )
    {
        fprintf( stderr, "Unable to lseek to the end\n" ) ;
        return NULL ;
    }

#if APM_DEBUG
    printf( "File length: %lld\n", fsize ) ;
#endif

    /* Go back to the beginning of the input file */
    if ( lseek(fd, 0, SEEK_SET) == -1 ) 
    {
        fprintf( stderr, "Unable to lseek to start\n" ) ;
        return NULL ;
    }

    /* Allocate data to copy the target text */
    buf = (char *)malloc( fsize * sizeof ( char ) ) ;
    if ( buf == NULL ) 
    {
        fprintf( stderr, "Unable to allocate %lld byte(s) for main array\n",
                fsize ) ;
        return NULL ;
    }

    n_bytes = read( fd, buf, fsize ) ;
    if ( n_bytes != fsize ) 
    {
        fprintf( stderr, 
                "Unable to copy %lld byte(s) from text file (%d byte(s) copied)\n",
                fsize, n_bytes) ;
        return NULL ;
    }

#if APM_DEBUG
    printf( "Number of read bytes: %d\n", n_bytes ) ;
#endif

    *size = n_bytes ;


    close( fd ) ;


    return buf ;
}


#define MIN3(a, b, c) ((a) < (b) ? ((a) < (c) ? (a) : (c)) : ((b) < (c) ? (b) : (c)))

int levenshtein(char *s1, char *s2, int len, int * column) {
    unsigned int x, y, lastdiag, olddiag;

    for (y = 1; y <= len; y++)
    {
        column[y] = y;
    }
    for (x = 1; x <= len; x++) {
        column[0] = x;
        lastdiag = x-1 ;
        for (y = 1; y <= len; y++) {
            olddiag = column[y];
            column[y] = MIN3(
                    column[y] + 1, 
                    column[y-1] + 1, 
                    lastdiag + (s1[y-1] == s2[x-1] ? 0 : 1)
                    );
            lastdiag = olddiag;

        }
    }
    return(column[len]);
}




int 
main( int argc, char ** argv )
{

  MPI_Init(&argc, &argv);
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
  int modulo = comm_size%nb_patterns;
  int countData = comm_size/nb_patterns; 
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


  printf( "Approximate Pattern Mathing: "
          "looking for %d pattern(s) in file %s w/ distance of %d\n", 
          nb_patterns, filename, approx_factor ) ;

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
      MPI_Send(n_matches_rank[j], 1, MPI_INT, 0, j, MPI_COMM_WORLD);
    }
    
  }

  else{
     /* Timer start */
    gettimeofday(&t1, NULL);

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
	MPI_Status sta;
	MPI_Recv(n_matches[j*comm_size+sender], 1, MPI_INT, sender, j, MPI_COMM_WORLD, &sta);
      }
    }
      
    for (j = 0 ; j < nb_patterns_rank ; j++){
      n_matches[comm_size*j] = n_matches_rank[j];
    }
    
    /* Timer stop */
    gettimeofday(&t2, NULL);
    duration = (t2.tv_sec -t1.tv_sec)+((t2.tv_usec-t1.tv_usec)/1e6);
    printf( "APM done in %lf s\n", duration ) ;
    
    for ( i = 0 ; i < nb_patterns_rank ; i++ ) {
      printf( "Number of matches for pattern <%s>: %d\n", 
              pattern[i], n_matches[i] ) ;
    }
    
  }


  /*****
   * END MAIN LOOP
   ******/

  MPI_Finalize();
  return 0 ;
}
