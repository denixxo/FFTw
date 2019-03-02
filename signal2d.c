# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h> 

# include "fftw3.h"

int main ( void );
void signal2d ( void );
double frand ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

Signal 2d: retrouve la transformée d'un champ de valeurs en 2D. 

Version décembre 2012
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "Signal2d\n" );
  printf ( " version C \n" );
  printf ( "  Requiert la biblio FFTW3.\n" );

      signal2d ( ); // test @ 2d 

/*
  fin
*/
  printf ( "\n" );
  printf ( "Signal 2d - fin \n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void signal2d ( void )

/******************************************************************************/
/*
  Purpose:

    Signal2d: apply FFT to complex 2D data.

  Discussion: voir la discussion pour le code Poisson.c (aussi 2D). 

*/
{
  int i;
  fftw_complex *in;
  fftw_complex *in2;
  int j;
  int nx = 100;
  int ny = 100;
  fftw_complex *out;
  fftw_plan plan_backward;
  fftw_plan plan_forward;
  unsigned int seed = 123456789;

  printf ( "\n" );
  printf ( "fonction signal2d - début \n" );
/*
Données de départ
*/
  in = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );

  srand ( seed );

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      in[i*ny+j][0] = sin(2*PI/(nx-1)*i)*cos(2*PI/(ny-1)*j) + frand()/10;
      in[i*ny+j][1] = 0.0;
    }
  }

  printf ( "\n" );
  printf ( "  En entrées:\n" );
  printf ( "\n" );

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      printf ( "  %4d  %4d  %12f  %12f\n", i, j, in[i*ny+j][0], in[i*ny+j][1] );
    }
  }
/*
  Transformée FFT : vers l'espace complex
*/
  out = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );

  plan_forward = fftw_plan_dft_2d ( nx, ny, in, out, FFTW_FORWARD, 
    FFTW_ESTIMATE );

  fftw_execute ( plan_forward );

  printf ( "\n" );
  printf ( "  FFT - les Coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      printf ( "  %4d  %4d  %12f  %12f\n", i, j, out[i*ny+j][0], out[i*ny+j][1] );
    }
  }
/*
Recouvre les données de départ : 
*/
  in2 = fftw_malloc ( sizeof ( fftw_complex ) * nx * ny );

  plan_backward = fftw_plan_dft_2d ( nx, ny, out, in2, FFTW_BACKWARD, 
    FFTW_ESTIMATE );

  fftw_execute ( plan_backward );

  printf ( "\n" );
  printf ( "  Données recouvrées: voir fichier Signal2d.dat\n" );
  printf ( "\n" );

  FILE *fp ; 

  fp = fopen( "Signal2d.dat", "w+") ; 

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      fprintf ( fp, "  %4d  %4d  %12f  %12f\n",  i, j,
	       in[i*ny+j][0], in2[i*ny+j][0] / ( double ) ( nx * ny ) );
      //      in2[i*ny+j][1] / ( double ) ( nx * ny ) );
    }
  }

  fclose(fp); 

/*
Libere l'espace utilise: 
*/
  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );

  fftw_free ( in );
  fftw_free ( in2 );
  fftw_free ( out );

  return;
}
//*****************************************************************************/

double frand ( void )

//*****************************************************************************/
/*
  Purpose:

    FRAND returns random values between 0 and 1.
*/
{
  double value;

  value = ( ( double ) rand ( ) / ( RAND_MAX ) );

  return value;
}
//*****************************************************************************/

void timestamp ( void )

/******************************************************************************/
/*
  Purpose:

    TIMESTAMP prints the current YMDHMS date as a time stamp.

  Example:

    31 May 2001 09:45:54 AM
*/
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  printf ( "%s\n", time_buffer );

  return;
# undef TIME_SIZE
}
