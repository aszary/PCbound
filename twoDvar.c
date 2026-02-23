// Estimating the spark motion in two dimension 
//
// To compile: gcc twoDvar.c -o ~/bin/twoDvar -lcpgplot -lm

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cpgplot.h"


# define ZERO_lf (0.000001) // round off zero error.
# define LINELEN (10000)     // Size of line.
# define WORDSIZE (100)      // Size of filename.

#define USAGE "\
   usage: twoDvar  ntime\n"


void cmdlineError (char *message);

void sparkconfig (float *x_arr, float *y_arr, double *x_sprk, double *y_sprk,
                     double *rad_sprk, double h_sprk, double h_drft, 
		     double D_cap, double theta_sp, double rad_trk, 
		     double x_cent, double y_cent, double del_theta, 
		     double theta_h, double theta_l, int *ncap);

int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       D_cap, h_sprk, h_drft;
   double       x_cent, y_cent;
   double      *x_sprk, *y_sprk, *rad_sprk;
   double       theta_sp, theta_c, del_theta, theta_h, theta_l;
   double       rad_in, rad_out, rad_trk;
   int          N_sprk, tot_sprk, N_trk, N_cell;
   int          ii, jj, ncel, ncnt, ntime, ncap;

   char        *label;
   float       *x_arr, *y_arr;
   float        xmin, xmax, ymin, ymax;


   if (argc < 2)
      cmdlineError (USAGE);

   ntime = atoi(argv[1]);

   label = (char *)calloc(WORDSIZE, sizeof(char));


   //Initializing parameters

   D_cap = 15.0;  // polar cap diameter in meters.
   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.


   x_cent = D_cap/2;
   y_cent = D_cap/2;


   // Sparks moving along outer rim
   rad_out = D_cap/2;	   
   rad_in = rad_out - h_sprk;

   N_sprk = 3.3*(rad_out*rad_out - rad_in*rad_in)/(h_sprk*h_sprk);

   theta_sp = 2*M_PI/N_sprk;

   N_cell = D_cap/h_drft;

   x_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   y_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   rad_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   
   x_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
   y_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
   
   //Finding the time evolution of sparks and plotting them
   cpgbeg(0, "/xs", 1, 1);

   rad_trk = 0.5*(rad_out + rad_in);
   del_theta = h_drft/rad_trk;
   theta_h = theta_sp/2;
   theta_l = theta_sp/2;

   for (ii=0;ii<ntime;ii++) {
    
      ncap = 0;

      sparkconfig (x_arr, y_arr, x_sprk, y_sprk, rad_sprk, h_sprk, h_drft, 
		     D_cap, theta_sp, rad_trk, x_cent, y_cent, del_theta, 
		     theta_h, theta_l, &ncap);


      //Plotting the sparks on the polar cap
      xmin = 0.0;
      xmax = (float)(D_cap);

      ymin = 0.0;
      ymax = (float)(D_cap);

      sprintf(label, "Iteration # %d", ii+1);

      cpgenv (xmin, xmax, ymin, ymax, 1, 0);

      cpglab ("X (m)", "Y (m)", label);

      cpgsfs(2);

      cpgcirc ((float)(x_cent), (float)(y_cent), (float)(D_cap/2));

      cpgpt (ncap, x_arr, y_arr, 1);


      //Initializing spark location for next time
      theta_h+=del_theta;

      if (theta_h > theta_sp) 
         theta_h-=theta_sp;

      theta_l+=del_theta;

      if (theta_l > theta_sp)
         theta_l-=theta_sp;
   }
   cpgclos();


   return 0;
}


void cmdlineError (char *message) {
   fprintf (stderr, "---------------");
   fprintf (stderr, "Error on Command Line. Exiting-----------\n");
   fprintf (stderr, "%s", message);
   exit (1);
}


void sparkconfig (float *x_arr, float *y_arr, double *x_sprk, double *y_sprk,
                     double *rad_sprk, double h_sprk, double h_drft, 
		     double D_cap, double theta_sp, double rad_trk, 
		     double x_cent, double y_cent, double del_theta, 
		     double theta_h, double theta_l, int *ncap)
{
   double      x_val, y_val, del_x, del_y, rad_val;
   double      theta_c, theta_fl, theta_fh;
   double      rad_fh, rad_lh;
   int         ii, ncnt, ncell, nsprk_l, nsprk_h;

	
   //Upper half, clockwise track 
   theta_c = M_PI - theta_h;
   ncnt = 0;
   nsprk_h = 0;

   while (theta_c >= 0) {

      x_sprk[ncnt] = x_cent + rad_trk*cos(theta_c);
      y_sprk[ncnt] = y_cent + rad_trk*sin(theta_c);
      rad_sprk[ncnt] = 0.5*h_sprk;


      //Spark at leading edge
      if (M_PI - theta_c <= theta_sp/2) {
         theta_fh = 0.5*(M_PI-theta_c) + theta_sp/4;
         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fh));
         rad_fh = rad_trk+h_sprk/2-rad_sprk[ncnt];

         x_sprk[ncnt] = x_cent + rad_fh*cos(M_PI-theta_fh);
         y_sprk[ncnt] = y_cent + rad_fh*sin(M_PI-theta_fh);
      }


      //Spark at trailing edge
      if (theta_c <= theta_sp/2) {
	 theta_fh = 0.5*(theta_c+theta_sp/2);
         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fh));
	 rad_fh = rad_trk+h_sprk/2-rad_sprk[ncnt];

         x_sprk[ncnt] = x_cent + rad_fh*cos(theta_fh);
         y_sprk[ncnt] = y_cent + rad_fh*sin(theta_fh);
      }


      theta_c-=theta_sp;
      ncnt++;
      nsprk_h++;
   }

   //lower half, anti-clockwise track
   theta_c = M_PI + theta_l;
   nsprk_l = 0;

   while (theta_c <= 2*M_PI) {

      x_sprk[ncnt] = x_cent + rad_trk*cos(theta_c);
      y_sprk[ncnt] = y_cent + rad_trk*sin(theta_c);
      rad_sprk[ncnt] = 0.5*h_sprk;

      
      //Spark at leading edge
      if (theta_c - M_PI <= theta_sp/2) {
         theta_fl = 0.5*(theta_c-M_PI) + theta_sp/4;
         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fl));
         rad_fh = rad_trk+h_sprk/2-rad_sprk[ncnt];

         x_sprk[ncnt] = x_cent + rad_fh*cos(M_PI+theta_fl);
         y_sprk[ncnt] = y_cent + rad_fh*sin(M_PI+theta_fl);
      }


      //Spark at trailing edge
      if (2*M_PI-theta_c <= theta_sp/2) {
         theta_fl = 0.5*(2*M_PI-theta_c) + theta_sp/4;
         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fl));
	 rad_lh = rad_trk+h_sprk/2-rad_sprk[ncnt];

	 x_sprk[ncnt] = x_cent + rad_lh*cos(2*M_PI-theta_fl);
	 y_sprk[ncnt] = y_cent + rad_lh*sin(2*M_PI-theta_fl);
      }

      theta_c+=theta_sp;
      ncnt++;
      nsprk_l++;
   }

   //Plugging gap in the beginning
   if (theta_l+theta_h > theta_sp) {
      rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(0.5*(theta_l+theta_h)) - h_sprk);
      x_sprk[ncnt] = x_cent - D_cap/2 + rad_sprk[ncnt];
      y_sprk[ncnt] = y_cent;

      ncnt++;
   }

   //Plugging gap in the end
   theta_fh = M_PI - theta_h - (nsprk_h-1)*theta_sp;
   theta_fl = M_PI - theta_l - (nsprk_l-1)*theta_sp;

   if (theta_fl+theta_fh > theta_sp && theta_fl+theta_fh < 2*theta_sp) {
      rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(0.5*(theta_fl+theta_fh)) - h_sprk);
      x_sprk[ncnt] = x_cent + D_cap/2 - rad_sprk[ncnt];
      y_sprk[ncnt] = y_cent;

      ncnt++;
   }


   //Section of the polar cap within each spark
   ncell = *ncap;
   x_val = h_drft/2;

   while (x_val < D_cap) {

      y_val = h_drft/2;

      while (y_val < D_cap) {

         for (ii=0;ii<ncnt;ii++) {
            del_x = x_val - x_sprk[ii];
	    del_y = y_val - y_sprk[ii];

            rad_val = sqrt(del_x*del_x + del_y*del_y);

            if (rad_val < rad_sprk[ii]) {
	       x_arr[ncell] = (float)(x_val);
	       y_arr[ncell] = (float)(y_val);
	       ncell++;
            }
	 }
         y_val+=h_drft;
      }
      x_val+=h_drft;
   }

   *ncap = ncell;
}
