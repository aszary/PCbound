//Estimating the ellipse parameters for the polar cap boundary
//
//To compile: gcc polcapfit.c -o ~/bin/polcapfit -lm 
//


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include "ellipse_fit.h"

# define LINELEN  (1000)      // Size of line.
# define WORDSIZE (100)      // Size of filename.


#define USAGE "\
   usage: polcapfit   pol_out\n"

void cmdlineError (char *message);

int splitstr (char *iline, char *strarray[], char *sep, int nsep);

void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out);


int main (int argc, char *argv[])
{
   FILE        *fpt;
   char        *capfile;
   char        *readline, *readval[20];
   double      *x_cap, *y_cap;
   double       x_cent, y_cent, R_maj, R_min, R_avg;
   double       el_fact, el_angl, el_sig;
   double       x_in, y_in, del_x, x_out, y_out;
   int          FixedCenter, FixedAngle, FixedEll;
   int          ii, ncap, nplt, ncnt, if_fit;


   capfile  = (char *)calloc(WORDSIZE, sizeof(char));
   readline = (char *)calloc(LINELEN, sizeof(char));

   if (argc < 2)
      cmdlineError (USAGE);

   sprintf (capfile, "%s", argv[1]);;

   //Reading the polar cap configuration
   if ((fpt = fopen(capfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", capfile);
      exit(1);
   }

   //Reading file size and initializing arrays
   fgets(readline, LINELEN, fpt);

   splitstr (readline, readval, " #=:", 4);

   ncap = atoi(readval[1]);

   fgets(readline, LINELEN, fpt);

   x_cap = (double *)calloc(ncap+1, sizeof(double));
   y_cap = (double *)calloc(ncap+1, sizeof(double));


   //Reading the Polar Cap boundary
   for (ii=0;ii<ncap;ii++) {
      fgets(readline, LINELEN, fpt);

      splitstr (readline, readval, " ", 1);

      x_cap[ii] = atof(readval[9]);
      y_cap[ii] = atof(readval[10]);
   }
   fclose(fpt);


   //Estimating the ellipticity fits
   FixedCenter = 0;
   FixedAngle = 0;
   FixedEll = 0;

   if_fit = FitEllipse (x_cap, y_cap, ncap, &x_cent, &y_cent, FixedCenter, 
		         FixedAngle, FixedEll, &R_maj, &R_min, &R_avg, 
			 &el_fact, &el_angl, &el_sig);


   fprintf (stderr, "x_cent %lf  y_cent %lf  R_maj %lf  R_min %lf  el_fact %lf  el_angl %lf\n", x_cent, y_cent, R_maj, R_min, el_fact, el_angl);

   //Writing out ellipse to plot
   fpt = fopen("out_elsp.dat", "w");

   del_x = R_maj/1000;

   x_in = -R_maj;

   while (x_in <  R_maj) {

      y_in = R_min*sqrt(1.0 - pow(x_in/R_maj, 2.0));

      coortrans (x_in, y_in, el_angl*M_PI/180, &x_out, &y_out);

      fprintf(fpt, "%lf  %lf\n", x_out+x_cent, y_out+y_cent);

      x_in+=del_x;
   }

   x_in = R_maj;

   while (x_in > -R_maj) {

      y_in = -R_min*sqrt(1.0 - pow(x_in/R_maj, 2.0));

      coortrans (x_in, y_in, el_angl*M_PI/180, &x_out, &y_out);

      fprintf(fpt, "%lf  %lf\n", x_out+x_cent, y_out+y_cent);

      x_in-=del_x;
   }

   fclose(fpt);


   return 0;
}


void cmdlineError (char *message) {
   fprintf (stderr, "---------------");
   fprintf (stderr, "Error on Command Line. Exiting-----------\n");
   fprintf (stderr, "%s", message);
   exit (1);
}


void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out)
{
   *x_out =  x_in*cos(theta) + y_in*sin(theta);

   *y_out = -x_in*sin(theta) + y_in*cos(theta);

}


int splitstr(char *iline, char *strarray[], char *sep, int nsep)
{
   int          ii, jj, kk, nword;
   int          llen;
   int          inword_f, if_sep;
   char         *newline;

   newline = (char *)calloc(LINELEN, sizeof(char));

   llen = strlen(iline);
   if (llen > LINELEN) {
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


