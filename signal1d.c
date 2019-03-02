# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h> 

# include "fftw3.h"

int main ( void );
void signal1d ( void );
double frand ( void );
void timestamp ( void );

/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Purpose:

Signal 1d: retrouve la transformée d'une fonction sinusoidale ..

version C : décembre 2012

*/
{
  timestamp ( );

      signal1d ( ); // 1d test
/*
  fin
*/
  printf ( "\n" );
  printf ( "fin - Signal 1d\n" );
  timestamp ( );

  return 0;
}
/******************************************************************************/

void signal1d ( void )

/******************************************************************************/
/*
  Purpose:

    Démonstration de la bibliothèque  FFTW pour donées (complexes) @  1D. 

  Discussion:

*/
{
  int i;
  fftw_complex *in;
  fftw_complex *in2;
  int n = 100;
  fftw_complex *out;
  fftw_plan plan_backward;
  fftw_plan plan_forward;
  unsigned int seed = 123456789;

  printf ( "\n" );
  printf ( "Signa1d: début \n" );

/*
Genere les donnees (1d) 
*/
  in = fftw_malloc ( sizeof ( fftw_complex ) * n );

  for ( i = 0; i < n; i++ )
  {
    in[i][0] = sin( 2*PI/n * i) ;
    in[i][1] = 0.0;
  }

  printf ( "\n" );
  printf ( " Rendu des chiffres:\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %12f  %12f\n", i, in[i][0], in[i][1] );
  }
/*
Genere la transformee de Fourier
*/
  out = fftw_malloc ( sizeof ( fftw_complex ) * n );

  plan_forward = fftw_plan_dft_1d ( n, in, out, FFTW_FORWARD, FFTW_ESTIMATE );

  fftw_execute ( plan_forward );

  printf ( "\n" );
  printf ( " FFT - les coefficients:\n" );
  printf ( "\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %12f  %12f\n", i, out[i][0], out[i][1] );
  }
/*
Transformée inverse: retrouve le signal de départ 
*/
  in2 = fftw_malloc ( sizeof ( fftw_complex ) * n );

  plan_backward = fftw_plan_dft_1d ( n, out, in2, FFTW_BACKWARD, FFTW_ESTIMATE );

  fftw_execute ( plan_backward );

  printf ( "\n" );
  printf ( "  Nouveau vecteur pour le signal:\n" );
  printf ( "\n" );
  printf ( "  Normalisation retrouvé en divisant par  n:\n" );

  for ( i = 0; i < n; i++ )
  {
    printf ( "  %3d  %12f  %12f\n", i, 
      in2[i][0] / ( double ) ( n ), in2[i][1] / ( double ) ( n ) );
  }
/*
Libère l'espace utilisé : 
*/
  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );

  fftw_free ( in );
  fftw_free ( in2 );
  fftw_free ( out );

  return;
}

double frand ( void )

//*****************************************************************************/
/*
  Purpose:

    FRAND returns random values between 0 and 1.

  D
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
