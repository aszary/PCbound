// Estimating the spark motion in two dimension 
//
// To compile: gcc twoDout.c -o ~/bin/twoDout -lcpgplot -lm

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cpgplot.h"
#include <time.h>


# define ZERO_lf (0.000001) // round off zero error.
# define LINELEN (10000)     // Size of line.
# define WORDSIZE (100)      // Size of filename.

#define USAGE "\
   usage: twoDdrift  ntime  D_cap (m)\n"


void cmdlineError (char *message);

void sparkconfig (float *x_arr, float *y_arr, double *x_sprk, double *y_sprk, 
                     double *rad_sprk, double *theta_sp, double *theta_h, 
		     double *theta_l, double h_sprk, double h_drft, 
		     double D_cap, int N_trk, int *ncap);

void delay(int number_of_seconds);


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       D_cap, h_sprk, h_drft;
   double      *x_sprk, *y_sprk, *rad_sprk;
   double      *theta_sp, *theta_h, *theta_l;
   double       theta_c, del_theta;
   double       rad_in, rad_out, rad_trk;
   int          N_sprk, Max_sprk, N_trk, N_cell;
   int          ii, jj, ncel, ncnt, ntime, ncap;

   char        *label;
   float       *x_arr, *y_arr;
   float        xmin, xmax, ymin, ymax;
   int          del_sec;


   if (argc < 3)
      cmdlineError (USAGE);

   ntime = atoi(argv[1]);
   D_cap = atof(argv[2]);

   label = (char *)calloc(WORDSIZE, sizeof(char));


   //Initializing parameters

   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.


   // Initialization of variables
   N_trk = D_cap/(2.0*h_sprk);

   theta_sp = (double *)calloc(N_trk+1, sizeof(double));
   theta_h = (double *)calloc(N_trk+1, sizeof(double));
   theta_l = (double *)calloc(N_trk+1, sizeof(double));


   rad_out = D_cap/2;	   
   rad_in = rad_out - h_sprk;
   Max_sprk = 0;

   for (ii=0;ii<N_trk;ii++) {
      
      N_sprk = 3.3*(rad_out*rad_out - rad_in*rad_in)/(h_sprk*h_sprk);

      theta_sp[ii] = 2*M_PI/N_sprk;

      Max_sprk+=N_sprk;

      theta_h[ii] = theta_sp[ii]/2;
      theta_l[ii] = theta_sp[ii]/2;

      rad_out-=h_sprk;
      rad_in-=h_sprk;
   }

   N_cell = D_cap/h_drft;

   x_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   y_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   rad_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   
   x_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
   y_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
   
   
   //Finding the time evolution of sparks and plotting them
   cpgbeg(0, "/xs", 1, 1);

   del_sec = 50;

   for (ii=0;ii<ntime;ii++) {
    
      //Estimating the Spark configuration 
      sparkconfig (x_arr, y_arr, x_sprk, y_sprk, rad_sprk, theta_sp, theta_h, 
		      theta_l, h_sprk, h_drft, D_cap, N_trk, &ncap);


      //Plotting the sparks on the polar cap
      xmin = 0.0;
      xmax = (float)(D_cap);

      ymin = 0.0;
      ymax = (float)(D_cap);

      sprintf(label, "Iteration # %d", ii+1);

      cpgenv (xmin, xmax, ymin, ymax, 1, 0);

      cpglab ("X (m)", "Y (m)", label);

      cpgsfs(2);

      cpgcirc ((float)(D_cap/2), (float)(D_cap/2), (float)(D_cap/2));

      cpgpt (ncap, x_arr, y_arr, 1);

      delay (del_sec);

      cpgask (0);

      //Initializing spark location for next time
      rad_out = D_cap/2;
      rad_in = rad_out - h_sprk;

      for (jj=0;jj<N_trk;jj++) {

	 rad_trk = 0.5*(rad_out + rad_in);
	 del_theta = h_drft/rad_trk;

         theta_h[jj]+=del_theta;

         if (theta_h[jj] > theta_sp[jj]) 
            theta_h[jj]-=theta_sp[jj];

         theta_l[jj]+=del_theta;

         if (theta_l[jj] > theta_sp[jj])
            theta_l[jj]-=theta_sp[jj];

	 rad_out-=h_sprk;
         rad_in-=h_sprk;
      }
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
		     double *rad_sprk, double *theta_sp, double *theta_h, 
		     double *theta_l, double h_sprk, double h_drft, 
		     double D_cap, int N_trk, int *ncap)
{
   double      trkrad_ul[N_trk], trkrad_ut[N_trk];
   double      trkrad_dl[N_trk], trkrad_dt[N_trk];
   double      x_val, y_val, del_x, del_y, x_cent, y_cent, rad_val;
   double      theta_c, theta_fl, theta_fh;
   double      rad_in, rad_out, rad_trk, rad_fh, rad_fl;
   int         ii, jj, ncnt, ncell, nsprk_l, nsprk_h;


   x_cent = D_cap/2;
   y_cent = D_cap/2;

   rad_out = D_cap/2;
   rad_in = rad_out - h_sprk;

   ncnt = 0;

   for (ii=0;ii<N_trk;ii++) {
   
      rad_trk = 0.5*(rad_out+rad_in);

      //Upper half, clockwise track
      theta_c = M_PI - theta_h[ii];
      nsprk_h = 0;

      trkrad_ul[ii] = rad_trk;
      trkrad_ut[ii] = rad_trk;
      trkrad_dl[ii] = rad_trk;
      trkrad_dt[ii] = rad_trk;

      while (theta_c >= 0) {

         x_sprk[ncnt] = x_cent + rad_trk*cos(theta_c);
         y_sprk[ncnt] = y_cent + rad_trk*sin(theta_c);
         rad_sprk[ncnt] = 0.5*h_sprk;

         //Spark at leading edge
         if (M_PI - theta_c <= theta_sp[ii]/2) {

            theta_fh = 0.5*(M_PI - theta_c) + theta_sp[ii]/4;

            rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fh));

            trkrad_ul[ii] = rad_trk+h_sprk/2-rad_sprk[ncnt];

            x_sprk[ncnt] = x_cent + trkrad_ul[ii]*cos(M_PI-theta_fh);
            y_sprk[ncnt] = y_cent + trkrad_ul[ii]*sin(M_PI-theta_fh);
         }

         //Spark at trailing edge
         if (theta_c <= theta_sp[ii]/2) {

	    theta_fh = 0.5*(theta_c + theta_sp[ii]/2);

            rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fh));

	    trkrad_ut[ii] = rad_trk + h_sprk/2 - rad_sprk[ncnt];

            x_sprk[ncnt] = x_cent + trkrad_ut[ii]*cos(theta_fh);
            y_sprk[ncnt] = y_cent + trkrad_ut[ii]*sin(theta_fh);
         }

         theta_c-=theta_sp[ii];
         ncnt++;
         nsprk_h++;
      }

      //lower half, anti-clockwise track
      theta_c = M_PI + theta_l[ii];
      nsprk_l = 0;

      while (theta_c <= 2*M_PI) {

         x_sprk[ncnt] = x_cent + rad_trk*cos(theta_c);
         y_sprk[ncnt] = y_cent + rad_trk*sin(theta_c);
         rad_sprk[ncnt] = 0.5*h_sprk;
      
         //Spark at leading edge
         if (theta_c - M_PI <= theta_sp[ii]/2) {

            theta_fl = 0.5*(theta_c - M_PI) + theta_sp[ii]/4;

            rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fl));

            trkrad_dl[ii] = rad_trk + h_sprk/2 - rad_sprk[ncnt];

            x_sprk[ncnt] = x_cent + trkrad_dl[ii]*cos(M_PI+theta_fl);
            y_sprk[ncnt] = y_cent + trkrad_dl[ii]*sin(M_PI+theta_fl);
         }

         //Spark at trailing edge
         if (2*M_PI-theta_c <= theta_sp[ii]/2) {

            theta_fl = 0.5*(2*M_PI-theta_c) + theta_sp[ii]/4;

            rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(theta_fl));

	    trkrad_dt[ii] = rad_trk + h_sprk/2 - rad_sprk[ncnt];

	    x_sprk[ncnt] = x_cent + trkrad_dt[ii]*cos(2*M_PI-theta_fl);
	    y_sprk[ncnt] = y_cent + trkrad_dl[ii]*sin(2*M_PI-theta_fl);
         }

         theta_c+=theta_sp[ii];
         ncnt++;
         nsprk_l++;
      }

      //Plugging gap in the beginning
      if (theta_l[ii] + theta_h[ii] >= theta_sp[ii]) {

         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(0.5*(theta_l[ii] + theta_h[ii])) 
			         - h_sprk);

         x_sprk[ncnt] = x_cent - rad_out + rad_sprk[ncnt];
         y_sprk[ncnt] = y_cent;

	 if (ii > 0) 
            x_sprk[ncnt]-=(0.5*(trkrad_ul[ii]+trkrad_dl[ii])-rad_out+h_sprk/2);

         ncnt++;
      }

      //Plugging gap in the end
      theta_fh = M_PI - theta_h[ii] - (nsprk_h-1)*theta_sp[ii];
      theta_fl = M_PI - theta_l[ii] - (nsprk_l-1)*theta_sp[ii];

      if (theta_fl+theta_fh >= theta_sp[ii] && 
		                 theta_fl+theta_fh <= 2*theta_sp[ii]) {

         rad_sprk[ncnt] = 0.5*(2*rad_trk*sin(0.5*(theta_fl+theta_fh)) - h_sprk);
      
         x_sprk[ncnt] = x_cent + rad_out - rad_sprk[ncnt];
         y_sprk[ncnt] = y_cent;

         if (ii > 0)
            x_sprk[ncnt]+=(0.5*(trkrad_ut[ii]+trkrad_dt[ii])-rad_out+h_sprk/2);

         ncnt++;
      }
      rad_out-=h_sprk;
      rad_in-=h_sprk;
   }

   //Core Spark
   x_sprk[ncnt] = x_cent;
   y_sprk[ncnt] = y_cent;
   rad_sprk[ncnt] = D_cap/2 - N_trk*h_sprk - 1.5*h_drft;
   ncnt++;


   //Section of the polar cap within each spark
   ncell = 0;
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


void delay(int number_of_seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds)
        ;
}
