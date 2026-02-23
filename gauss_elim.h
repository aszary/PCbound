#ifndef PI
#define PI 3.14159265359
#endif
#ifndef Boolean
#define Boolean char
#endif
#ifndef SUCCESS
#define SUCCESS 1
#endif
#ifndef ERROR
#define ERROR 0
#endif
#ifndef True
#define True 1
#define False 0
#endif
 
 
#define ACCURACY 1.0e-12
 
 
 
/*
***********************************************************
       SOLVE LINEAR SYSTEM VIA GAUSSIAN ELIMINATION
 
 Source code taken from "Numerical Recipes" in Java at
    http://introcs.cs.princeton.edu/java/95linear/GaussianElimination.java.html
 
 and ported to C by Jean Souviron 2016/01/06
 
 
***********************************************************
*/
int LinearSystemByGaussian (double *A, double *X, int Order )
{
  int    N=Order, Length=Order+1, max, i, p, j ;
  double t, sum, alpha ;
 
 
  for (p = 0; p < N; p++) {
 
    /* find pivot row and swap */
    max = p;
    for (i = p + 1; i < N; i++) {
      if (fabs(A[i*Length+p]) > fabs(A[max*Length+p])) {
	max = i;
      }
    }
 
    /* Exchange rows */
    for ( i = 0 ; i < Length ; i++ )
      {
	t = A[p*Length+i] ;
	A[p*Length+i] = A[max*Length+i] ;
	A[max*Length+i] = t ;
      }
 
 
    /* singular or nearly singular */
    if (fabs(A[p*Length+p]) <= ACCURACY) {
      for ( j = 0 ; j < Order ; j++ )
	X[j] = 0.0 ;
      return 0 ;
    }
 
    /* pivot within A and b */
    for (i = p + 1; i < N; i++) {
      alpha = A[i*Length+p] / A[p*Length+p];
      A[i*Length+N] -= alpha * A[p*Length+N];
      for (j = p; j < N; j++) {
	A[i*Length+j] -= alpha * A[p*Length+j];
      }
    }
  }
 
  /* back substitution */
  for (i = N - 1; i >= 0; i--) {
    sum = 0.0;
    for (j = i + 1; j < N; j++) {
      sum += A[i*Length+j] * X[j];
    }
    X[i] = (A[i*Length+N] - sum) / A[i*Length+i];
  }
 
  return 1;
}
 
 
/*
***********************************************************************
      GENERALIZED LEAST-SQUARE FOR ANY NUMBER OF PARAMETERS
 
                      (Jean SOUVIRON - DAO Victoria - 1985)
 
         (Method of R. SEDGEWICK in "Algorithms"
                                    addison-wesley publishing co. 1984)
 
      Modified 2015/12/12 Jean Souviron
         Porting to C
	 Upgrading to take any number of parameters
 
-----------------------------------------------------------------------
 
      INPUT :
                 NDims = number of different parameters+1 (for the data to be 
		         fitted to)
                  NPts = total number of points
           InputMatrix = array of parameters for each point
                         (last vector has to be the data to be fitted)
      OUTPUT :
          OutputCoeffs = array of coefficients (dimension NDims-1 )
                 Sigma
 
      RETURN VALUE :
           1 if everything is OK
           0 if an error occurred (memory or singular array)
 
 
  Example : if we want to fit an ellipse with a group of points according to
            equation x2/a2 + y2/b2 = 1
 
            the input matrix will be N pts long, and 3 columns wide, with for
            each point the x2, the y2, and 1.
            the output coeffs will be 2 items, 1/a2 and 1/b2
            the sigma will be computed by applying these coeffs to the first 2 
	    columns for each point, and comparing with the result (1).
            
 
***********************************************************************
*/
int LeastSquares ( double *InputMatrix, int NDims, int NPts, 
		   double *OutputCoeffs, double *Sigma )
{
  double *Coeffs=NULL ;
  double  T, Z1, Z2 ;
  int     i,j, k, io = 0, ko, NParams ;
 
 
  NParams = NDims-1 ;
 
/*
  Allocations
*/
  if ( (Coeffs = (double *)calloc(NParams*NDims, sizeof(double))) == NULL )
    {
      return 0 ;
    }
 
/*
  Initializations
*/
  io = 0 ;
  for ( i = 0 ; i < NParams ; i++ )
    {
      OutputCoeffs[i] = 0.0 ;
 
      for ( j = 0 ; j < NDims ; j++ )
	{
	  T = 0.0 ;
 
	  for ( k = 0 ; k < NPts ; k++ )
	    {
	      ko = k * NDims ;
 
	      Z1 = InputMatrix[ko+i] ;
	      Z2 = InputMatrix[ko+j] ;
 
	      T = T + (Z1*Z2) ;
	    }
 
	  Coeffs[io+j]= T ;
	}
      io = io + NDims ;
    }
 
/*
-----------------
   Solving by gaussian elimination
-----------------
*/
  i = LinearSystemByGaussian ( Coeffs, OutputCoeffs, NParams ) ;
 
  free ( Coeffs );
  if ( i != 1 )
    return 0 ;
 
/*
-----------------
  Computes sigma
-----------------
*/
  *Sigma = 0.0 ;
  for ( i = 0 ; i < NPts ; i++ )
    { 
      Z1 = InputMatrix[i*NDims+NParams]; /* The value to be compared with */
 
      for ( j = 0 ; j < NParams ; j++ )
	{
	  Z1= Z1 - (OutputCoeffs[j]*InputMatrix[i*NDims+j]) ;
	}
 
      if( Z1 < 0.0 ) 
	Z1 = -Z1;
 
      *Sigma = *Sigma + Z1 ;
    }
 
  *Sigma = (*Sigma/(double)NPts)*sqrt(PI/2.0) ;
 
  return 1 ;
}
