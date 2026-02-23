#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double urandom(long *idum)
{
   int j;
   long k;
   static long iy=0;
   static long iv[NTAB];
   double temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1)
         *idum=1;
      else
         *idum = -(*idum);
      for (j=NTAB+7;j>=0;j--) {
         k=(*idum)/IQ;
         *idum=IA*(*idum-k*IQ)-IR*k;
         if (*idum < 0)
            *idum += IM;
         if (j < NTAB)
            iv[j] = *idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if (*idum < 0)
      *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j] = *idum;
   if ((temp=AM*iy) > RNMX)
      return (RNMX);
   else
      return temp;
}

#undef IA
#undef IM
#undef AM
#undef IQ
#undef IR
#undef NTAB
#undef NDIV

void cmdlineError (char *message) {
   fprintf (stderr, "---------------");
   fprintf (stderr, "Error on Command Line. Exiting-----------\n");
   fprintf (stderr, "%s", message);
   exit (1);
}


int splitstr(char *iline, char *strarray[], char *sep, int nsep)
{
   int          ii, jj, kk, nword;
   int          llen;
   int          inword_f, if_sep;
   char         *newline;

   newline = (char *)calloc(10000, sizeof(char));

   llen = strlen(iline);
   if (llen > 10000) {
      fprintf (stderr, "Input line longer than 10000 characters. ");
      fprintf (stderr, "Exiting\n");
      exit (1);
   }

   ii = 0;
   jj = 0;
   nword = 0;
   inword_f = 0;
   while (ii < llen) {
      if_sep = 0;
      for (kk=0;kk<nsep;kk++) {
         if (*(iline+ii) == *(sep+kk))
            if_sep = 1;
      }
      if (if_sep == 0) {
         *(newline+jj) = *(iline+ii);
         if (inword_f == 0) {
            inword_f = 1;
            strarray[nword++] = newline+jj;
         }
         jj++;
      }
      else {
         if (inword_f == 1) {
            inword_f = 0;
            *(newline+jj) = '\0';
            jj++;
         }
      }
      ii++;
   }

   return (nword);
}


// arrstat[0] = mean; arrstat[1] = rms;
void meanrms(double *inarray, int arrsize, double *arrstat)
{
   int            ii;
   double         x2;


   if (arrsize <= 1 ) {
      fprintf(stderr,"array size <=1 in stat calc....exiting\n");
      exit (1);
   }

   arrstat[0] = 0.0;
   arrstat[1] = 0.0;
   x2 = 0.0;

   for (ii=0;ii<arrsize;ii++)
      arrstat[0]+= inarray[ii];

   arrstat[0]/=(double)(arrsize);

   for (ii=0;ii<arrsize;ii++)
      x2+= (inarray[ii]-arrstat[0])*(inarray[ii]-arrstat[0]);

   arrstat[1] = sqrt(x2/(double)(arrsize-1));
}


double findedge (double x1, double y1, double x2, double y2, double val)
{
   double       m_ln, b_ln;

   m_ln = (y2-y1)/(x2-x1);
   b_ln = y2-m_ln*x2;

   return (val-b_ln)/m_ln;
}


double findbet (double xi_1, double xi_2, double x_i, double xf_1, double xf_2)
{
   double       del_x;

   
   del_x = (x_i - xi_1)/(xi_2 - xi_1);

   return del_x*(xf_2 - xf_1) + xf_1;
}


int ifbound (double theta, double phi, double *outfield, double angsep,
              int nfield)
{
   double       angle;
   int          ii, if_bound;

   if_bound = 0;
   for (ii=0;ii<nfield;ii++) {
      angle = acos(cos(theta)*cos(outfield[6*ii+4])+
                   sin(theta)*sin(outfield[6*ii+4])*cos(phi-outfield[6*ii+5]));

      if (angle < 2.0*angsep)
         if_bound = 1;
   }

   return if_bound;
}


void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out)
{
   *x_out =  x_in*cos(theta) + y_in*sin(theta);

   *y_out = -x_in*sin(theta) + y_in*cos(theta);

}

// Rectangular to polar representation
void rec2pol (double x, double y, double *r, double *theta)
{
  *r = sqrt(x*x+y*y);
  *theta = atan2(y,x)*360/2.0/M_PI;
}


void findbase(double *inarray, int nsize, int basefact, int *spos, int *epos,
               double *minstat)
{

   int          ii, jj, base_length, l1, ra1, ra2;
   double       pass_arr[nsize], arrstat[2];


   base_length = nsize/basefact;
   minstat[1] = 1e+32;

   for (ii=0;ii<nsize-base_length;ii++) {
      ra1 = ii;
      ra2 = ii + base_length - 1;

      for (jj=ra1;jj<=ra2;jj++)
         pass_arr[jj-ra1] = inarray[jj];

      meanrms(pass_arr, ra2-ra1+1, arrstat);

      if (arrstat[1] < minstat[1]){
         minstat[0] = arrstat[0];
         minstat[1] = arrstat[1];
         l1 = ra1;
      }
   }

   *spos = l1;
   *epos = l1 + base_length - 1;
}


double noise (double sigma, long *seed)
{
   double     v1, v2, gauss_rand, Rad_sqr;
   int        ii;

   do
   {
      //obtaining the uniform random numbers between -1 and 1
      v1 = 2.0*urandom(seed) - 1.0;
      v2 = 2.0*urandom(seed) - 1.0;

      Rad_sqr = v1*v1 + v2*v2 ;

   } while(Rad_sqr >= 1.0 || Rad_sqr == 0.0);

   //Generating gaussian random numbers
   gauss_rand = sigma*(sqrt(-2.0*log(Rad_sqr)/Rad_sqr)*v1);

   return  gauss_rand;
}


