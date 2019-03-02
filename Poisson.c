# include <stdlib.h>
# include <stdio.h>
# include <time.h>
# include <math.h> 

# include "fftw3.h"

# define PI 3.14159

int main ( void );

void poisson ( void );
void test04 ( void );
double frand ( void );
void timestamp ( void );


/******************************************************************************/

int main ( void )

/******************************************************************************/
/*
  Poisson: retrouve le potentiel d'une distribution de masse dans l'espace avec FFTW. 

  Modified:  [ ] 

  Author: John Burkardt, Christian Boily. 
*/
{
  timestamp ( );

  printf ( "\n" );
  printf ( "Code Poisson \n" );
  printf ( "  version C \n" );
  printf ( "  Requiert la bibliothèque  FFTW3 .\n" );

  poisson ( ); /* 2d test */
/*
  Fin. 
*/
  printf ( "\n" );
  printf ( "Poisson - fin\n" );
  timestamp ( );

  return 0;
}

void poisson ( void )

/******************************************************************************/
/*
  Discussion:

L'espace reel est décompose en nx x ny mailles où sont distribuées les sources (S) 
amenant le potentiel Phi.

La résolution de l'équation de Poisson se fait par transformée de Fourier et évoque 
le théorème de la convolution: 

Si Phi(x) = int( G(x-x') S(x') dx'  

alors la transformée de Fourier du potentiel est 

   Phik(k) = Gk(k) * Sk(k) 

avec le souscripte "k" indiquant une représentation dans l'espace des vecteurs d'onde 
k. Le potentiel Phi(x) s'obtient d'une transformée de Fourier "inverse", de l'espace 
de nombre d'onde "k", vers l'espace réel x. 

Le programme fourni ici se limite au cas à deux dimensions libres, x et y. Les indices 
utilisés sont alors i et j (pour le numéro de la maille dans chacune des dimensions). 

Attention - nous utilisons la bibliothèque de référence FFTW versions 3.3.3 (2011). 
Les fonctions C font référence aux indices de nombres complexes (toutes les quantités 
sont complexes par défault), par exemple pour A :  

      a[i*ny+j][0]      partie réeelle de A(I,J).
      a[i*ny+j][1]      composante imaginaire de  A(I,J). 

et ainsi de suite. 

Version: décembre 2012 / suivant un patron de John Burkardt. 

*/
{

  int i;
  fftw_complex *S;
  fftw_complex *G;
  fftw_complex *Phi;

  int j;
  // Fixe la dimension de l'espace réel en terme de nx x ny mailles: 
  //
  int nx = 64;
  int ny = 64;

  fftw_complex *Sk;
  fftw_complex *Gk;  
  fftw_complex *Phik;

  fftw_plan plan_backward;
  fftw_plan plan_forward;

  unsigned int seed = 123456789;

  double h = 1.0; 

  printf ( "\n" );
  printf ( " fonction Poisson\n" );
  printf ( "  Demontre les transformée de Fourier par  FFTW3 sur une maille de  %d x %d données complexes.\n",
    nx, ny );
  printf ( "\n" );
/*
  Crée les tableaux des Sources (S) et de la fonciton de Green (G)
  Dimensions de l'espace "physique" : nx x ny mais dimension de l'espace de Fourier = 2nx x 2ny. 

*/
  G = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );
  S = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );

  srand ( seed );

  /* D'abord les sources (forme de densité quelconque) 
  */ 

  double r; 

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      r = h * sqrt( (i-nx/2)*(i-nx/2) + (j-ny/2)*(j-ny/2) ) ; 
      S[i*2*ny+j][0] =  0.0 * exp( -r/h ) ; // Choix 1: Profil de densité arbitraire mais axial
      S[i*2*ny+j][1] = 0.0;
    }
  }

   S[nx/2*2*ny+ny/2][0] = 4*PI*1.0/PI/h/h; // Choix 2: Point-masse = 1 à l'origine des coordonnées

  // Double l'espace (éliminer le biais de la FFT): sources nulles partout ailleurs 
  //
  for ( i = nx; i < 2*nx; i++ )
  {
    for ( j = ny; j < 2*ny; j++ )
    {
      r = h * sqrt( (i-nx/2)*(i-nx/2) + (j-ny/2)*(j-ny/2) ) ; 
      S[i*2*ny+j][0] = 0.0 ; 
      S[i*2*ny+j][1] = 0.0;
    }
  }

  /* Ensuite: la fonction de Green 
  */ 
  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
      {   r = h * fmax(sqrt( (double)(i)*(i) + (double)(j)*(j) ), 1.0) ; 
      	G[i*2*ny+j][0] = - h*h / (4*PI) / r ;
        G[i*2*ny+j][1] = 0.0; 

   G[(2*nx-1-i)*2*ny+j][0] = G[i*2*ny+(2*ny-1-j)][0] = G[(2*nx-1-i)*2*ny+(2*ny-1-j)][0] = G[i*2*ny+j][0] ; 
   G[(2*nx-1-i)*2*ny+j][1] = G[i*2*ny+(2*ny-1-j)][1] = G[(2*nx-1-i)*2*ny+(2*ny-1-j)][1] = G[i*2*ny+j][1];
    }
  }

/*
  Transformées de Fourier, sources et fonction de Green.
*/
  Sk = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );
  Gk = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );

  plan_forward = fftw_plan_dft_2d ( 2*nx, 2*ny, S, Sk, FFTW_FORWARD, FFTW_ESTIMATE );
  fftw_execute ( plan_forward );
  plan_forward = fftw_plan_dft_2d ( 2*nx, 2*ny, G, Gk, FFTW_FORWARD, FFTW_ESTIMATE );
  fftw_execute ( plan_forward );

  /* Construction du potentiel dans l'espace de Fourier: application du Th. de la convolution
     Phik = Gk * Sk (quantités complexes en principe). 
  */

 Phik = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );
 
  for ( i = 0; i < 2*nx; i++ )
  {
    for ( j = 0; j < 2*ny; j++ )
    {
      Phik[i*2*ny+j][0] = Gk[i*2*ny+j][0] * Sk[i*2*ny+j][0] -  Gk[i*2*ny+j][1] * Sk[i*2*ny+j][1] ;
      Phik[i*2*ny+j][1] = Gk[i*2*ny+j][1] * Sk[i*2*ny+j][0] +  Gk[i*2*ny+j][0] * Sk[i*2*ny+j][1]; 
 
    } }

/*
  Recouvre le potentiel dans l'espace physique: FFT inverse (= backward) 
*/
  Phi = fftw_malloc ( sizeof ( fftw_complex ) * 2*nx * 2*ny );

  plan_backward = fftw_plan_dft_2d ( 2*nx, 2*ny, Phik, Phi, FFTW_BACKWARD, FFTW_ESTIMATE );
  fftw_execute ( plan_backward );

  printf ( "\n" );
  printf ( "  Normalisation des résultats (diviser par nx * ny): voir Poisson.dat\n" );
  printf ( "\n" );

  FILE *fp ; fp = fopen( "Poisson.dat", "w+") ; 

  /* Note: pour la normalisation "finale" on pose que la grille nx h = 3 unités de longeur, ce qui 
     change la norme du potentiel 
  */ 

  for ( i = 0; i < nx; i++ )
  {
    for ( j = 0; j < ny; j++ )
    {
      fprintf ( fp, "  %4d  %4d  %12f  %12f %12f \n",  i, j, \
		S[i*2*ny+j][0], Phi[i*2*ny+j][0]/(double)(nx*ny) * (nx*h/2/3.),  Phi[i*2*ny+j][1]/(double)(nx*ny)* (nx*h/3.) );
    }
  }

  fclose(fp); 

/*
  Libere l'espace mémoire: 
*/
  fftw_destroy_plan ( plan_forward );
  fftw_destroy_plan ( plan_backward );

  fftw_free (S); fftw_free (Sk );
  fftw_free (G); fftw_free (Gk);
  fftw_free (Phi); fftw_free (Phik);

  return;
}
/******************************************************************************/

//*****************************************************************************/

double frand ( void )

//*****************************************************************************/
/*
  Purpose:

    FRAND returns random values between 0 and 1.

  Modified:

    23 October 2005

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
