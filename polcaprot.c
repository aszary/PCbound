//Estimating the ellipse parameters for the polar cap boundary
//
//To compile: gcc polcaprot.c -o ~/bin/polcaprot -lm 
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

double dotprod (double *vect1, double *vect2);

void crossprod (double *vect1, double *vect2, double *crosvect);

void rotvect (double *invect, double *rotaxis, double *transvect, 
		double th_rot);

void unitvect (double *invect, double *transvect);


int main (int argc, char *argv[])
{
   FILE        *fpt;
   char        *capfile;
   char        *readline, *readval[20];
   double      *in_cap, *rot_cap;
   double      *r_cap, *th_cap, *ph_cap;
   double       inarr[3], rotaxis[3], trnsvect[3];
   double       x_cent, y_cent, z_cent;
   double       RS, th_c, ph_c, th_min, th_max, ph_min, ph_max;
   int          ii, ncap;


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

   //Reading file size
   fgets(readline, LINELEN, fpt);

   splitstr (readline, readval, " #=:", 4);

   ncap = atoi(readval[1]);

   fgets(readline, LINELEN, fpt);


   //Reading centre location
   r_cap = (double *)calloc(ncap+1, sizeof(double));
   th_cap = (double *)calloc(ncap+1, sizeof(double));
   ph_cap = (double *)calloc(ncap+1, sizeof(double));
   in_cap = (double *)calloc(3*ncap+1, sizeof(double));
   rot_cap = (double *)calloc(3*ncap+1, sizeof(double));


   //Reading the Polar Cap boundary
   for (ii=0;ii<ncap;ii++) {
      fgets(readline, LINELEN, fpt);

      splitstr (readline, readval, " ", 1);

      r_cap[ii] = atof(readval[6]);
      th_cap[ii] = atof(readval[7]);
      ph_cap[ii] = atof(readval[8]);

      in_cap[3*ii] = atof(readval[9]);
      in_cap[3*ii+1] = atof(readval[10]);
      in_cap[3*ii+2] = atof(readval[11]);
 
   }
   fclose(fpt);

   //Finding minimum maximum and central angles
   th_min = th_cap[0];
   th_max = th_cap[0];
   ph_min = ph_cap[0]; 
   ph_max = ph_cap[0];

   for (ii=1;ii<ncap;ii++) {

      if (th_min > th_cap[ii])
         th_min = th_cap[ii];
      if (th_max < th_cap[ii])
         th_max = th_cap[ii];

      if (ph_min > ph_cap[ii])
         ph_min = ph_cap[ii];
      if (ph_max < ph_cap[ii])
         ph_max = ph_cap[ii];
   }

   th_c = (th_min+th_max)/2;
   ph_c = (ph_min+ph_max)/2;

   fprintf (stderr, "%lf  %lf  %lf  %lf\n", th_min, th_max, th_c, ph_c);


   //Estimating rotation axis   
   RS = 10000;  //stellar radius in meters

   x_cent = RS*sin(th_c)*cos(ph_c);
   y_cent = RS*sin(th_c)*sin(ph_c);
   z_cent = RS*cos(th_c);

   inarr[0] = sin(th_c)*cos(ph_max) - sin(th_c)*cos(ph_min);
   inarr[1] = sin(th_c)*sin(ph_max) - sin(th_c)*sin(ph_min);
   inarr[2] = 0.0; 

   unitvect (inarr, rotaxis);

   //Rotating the polar cap and writing out

   fpt = fopen("rot_cap.dat","w");
 
   for (ii=0;ii<ncap;ii++) {
      
      inarr[0] = in_cap[3*ii] - x_cent;
      inarr[1] = in_cap[3*ii+1] - y_cent;
      inarr[2] = in_cap[3*ii+2] - z_cent;

      rotvect (inarr, rotaxis, trnsvect, -th_c);

      fprintf (fpt, "%lf  %lf  %lf    %lf  %lf  %lf\n", in_cap[3*ii], in_cap[3*ii+1], in_cap[3*ii+2], trnsvect[0] + x_cent, trnsvect[1] +y_cent, trnsvect[2] + z_cent);
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


void rotvect (double *invect, double *rotaxis, double *transvect, double th_rot)
{
   double       crosvect[3], scal_prod;
   int          ii;
   

   crossprod (rotaxis, invect, crosvect);

   scal_prod = dotprod (invect, rotaxis);

   for (ii=0;ii<3;ii++) 
      transvect[ii] = invect[ii]*cos(th_rot) + crosvect[ii]*sin(th_rot) + 
	                scal_prod*rotaxis[ii]*(1.0 - cos(th_rot));

}

double dotprod (double *vect1, double *vect2)
{
   double       scal;
   int          ii;

   scal = 0;

   for (ii=0;ii<3;ii++) 
      scal+=vect1[ii]*vect2[ii];

   return scal;
}

void crossprod (double *vect1, double *vect2, double *crosvect)
{
   crosvect[0] = vect1[1]*vect2[2] - vect1[2]*vect2[1];

   crosvect[1] = vect1[2]*vect2[0] - vect1[0]*vect2[2];

   crosvect[2] = vect1[0]*vect2[1] - vect1[1]*vect2[0];
      
}

void unitvect (double *invect, double *transvect)
{
   double       scal;
   int          ii;


   scal = 0.0;

   for (ii=0;ii<3;ii++)
      scal+=invect[ii]*invect[ii];

   for (ii=0;ii<3;ii++)
      transvect[ii] = invect[ii]/sqrt(scal);

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


