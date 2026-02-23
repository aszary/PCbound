// Estimating the spark motion in two dimensional elliptical polar cap
//
// To compile: gcc elipsDrift.c -o ~/bin/elipsDrift -lcpgplot -lm

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
   usage: elipsDrift  ntime  incl_angle (deg)  maj_axis (m)  ellip (b/a) corot_angle (deg) \n"


void cmdlineError (char *message);

void sparkconfig (float *x_arr, float *y_arr, double *th_sprk_u, 
		     double *th_sprk_d, double *x_sprk, double *y_sprk, 
                     double *a_sprk, double *theta_sp, double h_sprk, 
		     double h_drft, double a_cap, double b_cap, double th_cap, 
		     double co_angl, double x_cent, double y_cent, 
		     double del_theta, int N_trk, int *N_up, int *N_dn, 
		     int trk_max, int *ncap);

void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out);

void delay(int number_of_seconds);


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       a_cap, b_cap, th_cap, co_angl;
   double       h_sprk, h_drft, a_sprk, b_sprk;
   double       x_cent, y_cent;
   double      *x_sprk, *y_sprk, *rad_sprk;
   double      *theta_sp, *th_sprk_u, *th_sprk_d;
   double       theta_c, del_theta;
   double       a_in, a_out, a_trk, b_in, b_out;
   double       x_in, y_in, x_trns, y_trns;
   int         *N_up, *N_dn;
   int          Max_sprk, N_sprk, N_trk, N_cell, ncel, ncnt, ntime, ncap;
   int          ii, jj, indx_u, indx_d, trk_max;

   char        *label;
   float       *x_arr, *y_arr, *x_elips, *y_elips;
   float        xmin, xmax, ymin, ymax;
   int          del_sec, N_elips;


   if (argc < 6)
      cmdlineError (USAGE);

   ntime = atoi(argv[1]);
   th_cap = atof(argv[2])/180*M_PI;
   a_cap = atof(argv[3]);
   b_cap = a_cap*atof(argv[4]);
   co_angl = atof(argv[5])/180*M_PI;

   label = (char *)calloc(WORDSIZE, sizeof(char));


   //Initializing parameters

   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.

   a_sprk = h_sprk;
   b_sprk = a_sprk*b_cap/a_cap;

   x_cent = sqrt(pow(a_cap*cos(th_cap),2) + pow(b_cap*sin(th_cap),2));
   y_cent = sqrt(pow(a_cap*sin(th_cap),2) + pow(b_cap*cos(th_cap),2));


   // Initialization of variables
   N_trk = b_cap/(2*b_sprk);

   theta_sp = (double *)calloc(N_trk+1, sizeof(double));
   N_up = (int *)calloc(N_trk+1, sizeof(int));
   N_dn = (int *)calloc(N_trk+1, sizeof(int));

   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   b_out = b_cap;
   b_in = b_out - 2*b_sprk;

   Max_sprk = 0;

   trk_max = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk)/2+1;

   for (ii=0;ii<N_trk;ii++) {
      
      N_sprk = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk);

      theta_sp[ii] = 2*M_PI/N_sprk;

      Max_sprk+=N_sprk;

      a_out-=2*a_sprk;
      a_in-=2*a_sprk;

      b_out-=2*b_sprk;
      b_in-=2*b_sprk;
   }

   N_cell = 2*a_cap/h_drft;

   fprintf (stderr, "%d  %d %d\n", N_trk, Max_sprk, N_cell);    


   //Initializing the spark configuration
   x_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   y_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   rad_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   th_sprk_u = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));
   th_sprk_d = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));
   
   x_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
   y_arr = (float *)calloc(N_cell*N_cell+1, sizeof(float));
  
   indx_u = 0;
   indx_d = 0;

   for (ii=0;ii<N_trk;ii++) {

      //Upper track
      th_sprk_u[indx_u] = M_PI - theta_sp[ii]/2 - co_angl;
      N_up[ii] = 1;

      while (th_sprk_u[indx_u+N_up[ii]-1] >= theta_sp[ii] - co_angl) {

         th_sprk_u[indx_u+N_up[ii]]=th_sprk_u[indx_u+N_up[ii]-1]-theta_sp[ii];
         N_up[ii]++;
      }
      indx_u+=trk_max;

      //Lower track
      th_sprk_d[indx_d] = M_PI + theta_sp[ii]/2 - co_angl;
      N_dn[ii] = 1;

      while (th_sprk_d[indx_d+N_dn[ii]-1] <= 2*M_PI - theta_sp[ii] - co_angl) {

         th_sprk_d[indx_d+N_dn[ii]]=th_sprk_d[indx_d+N_dn[ii]-1]+theta_sp[ii];
         N_dn[ii]++;
      }
      indx_d+=trk_max;
   }


   //The outline of the elliptical polar cap
   N_elips = 100;

   x_elips = (float *)calloc(N_elips+1, sizeof(float));
   y_elips = (float *)calloc(N_elips+1, sizeof(float));

   theta_c = 0;

   del_theta = 2*M_PI/N_elips;

   ii = 0;

   while (theta_c < 2*M_PI) {
      
      x_in = a_cap*cos(theta_c);
      y_in = b_cap*sin(theta_c);

      coortrans (x_in, y_in, -th_cap, &x_trns, &y_trns);

      x_elips[ii] = (float)(x_trns + x_cent);
      y_elips[ii] = (float)(y_trns + y_cent);
      ii++;

      theta_c+=del_theta;
   }
   x_elips[ii] = x_elips[0];
   y_elips[ii] = y_elips[0];


   //Finding the time evolution of sparks and plotting them
   cpgbeg(0, "/xs", 1, 1);

   del_sec = 50;

   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   for (ii=0;ii<ntime;ii++) {
    
      //Estimating the Spark configuration 
      sparkconfig (x_arr, y_arr, th_sprk_u, th_sprk_d, x_sprk, y_sprk, 
		      rad_sprk, theta_sp, h_sprk, h_drft, a_cap, b_cap, th_cap,
		      co_angl, x_cent, y_cent, del_theta, N_trk, N_up, N_dn, 
		      trk_max, &ncap);


      //Plotting the sparks on the polar cap
      xmin = 0.0;
      xmax = (float)(2*x_cent);

      ymin = 0.0;
      ymax = (float)(2*y_cent);

      sprintf(label, "Iteration # %d", ii+1);

      cpgenv (xmin, xmax, ymin, ymax, 1, 0);

      cpglab ("X (m)", "Y (m)", label);

      cpgsfs(2);

      cpgline(N_elips+1, x_elips, y_elips);

      cpgpt (ncap, x_arr, y_arr, 1);

      delay (del_sec);

      cpgask (0);


      //Initializing spark location for next time
      indx_u = 0;
      indx_d = 0;

      for (jj=0;jj<N_trk;jj++) {

         a_trk = 0.5*(a_out+a_in);

         del_theta = h_drft/(a_trk);

         //Upper track
         th_sprk_u[indx_u]-=del_theta;

         if (th_sprk_u[indx_u] < M_PI - theta_sp[jj] - co_angl) 
            th_sprk_u[indx_u]+=theta_sp[jj];

         N_up[jj] = 1;

         while (th_sprk_u[indx_u+N_up[jj]-1] >= theta_sp[jj] - co_angl) {

            th_sprk_u[indx_u+N_up[jj]] = th_sprk_u[indx_u+N_up[jj]-1]-theta_sp[jj];
            N_up[jj]++;
         }
         indx_u+=trk_max;

         //Lower track
	 th_sprk_d[indx_d]+=del_theta;

         if (th_sprk_d[indx_d] > M_PI + theta_sp[jj] - co_angl)
            th_sprk_d[indx_d]-=theta_sp[jj];

	 N_dn[jj] = 1;

         while (th_sprk_d[indx_d+N_dn[jj]-1] <= 2*M_PI-theta_sp[jj]-co_angl) {

            th_sprk_d[indx_d+N_dn[jj]] = th_sprk_d[indx_d+N_dn[jj]-1]+theta_sp[jj];
            N_dn[jj]++;
         }
         indx_d+=trk_max;

	 a_trk-=2.0*h_sprk;
      }
   }
   cpgclos();


   return 0;
}


void cmdlineError (char *mssage) {
   fprintf (stderr, "---------------");
   fprintf (stderr, "Error on Command Line. Exiting-----------\n");
   fprintf (stderr, "%s", mssage);
   exit (1);
}

void sparkconfig (float *x_arr, float *y_arr, double *th_sprk_u, 
		     double *th_sprk_d, double *x_sprk, double *y_sprk, 
                     double *a_sprk, double *theta_sp, double h_sprk, 
		     double h_drft, double a_cap, double b_cap, double th_cap, 
		     double co_angl, double x_cent, double y_cent, 
		     double del_theta, int N_trk, int *N_up, int *N_dn, 
		     int trk_max, int *ncap)
{
   double      trk_a_ul[N_trk], trk_b_ul[N_trk];
   double      trk_a_ut[N_trk], trk_b_ut[N_trk];
   double      trk_a_dl[N_trk], trk_b_dl[N_trk];
   double      trk_a_dt[N_trk], trk_b_dt[N_trk];
   double      a_in, a_out, b_in, b_out, e_trk, b_sprk;
   double      x_val, y_val, x_in, y_in, x_trns, y_trns, el_val;
   double      theta_c, theta_fl, theta_fh;
   double      a_trk, b_trk, rad_fh, rad_fl;
   int         indx_u, indx_d, indx_us;
   int         ii, jj, ncnt, ncell, nsprk_l, nsprk_h;


   a_out = a_cap;
   a_in = a_out - 2.0*h_sprk;

   b_out = b_cap;
   b_in = b_out - b_cap/a_cap*2.0*h_sprk; 

   ncnt = 0;

   indx_u = 0;
   indx_d = 0;

   indx_us = 0;

   for (ii=0;ii<N_trk;ii++) {
   
      a_trk = 0.5*(a_out+a_in);
      b_trk = 0.5*(b_out+b_in);

      //Upper half, clockwise track
      for (jj=0;jj<N_up[ii];jj++) {

         x_in = a_trk*cos(th_sprk_u[indx_u+jj] - th_cap);
	 y_in = b_trk*sin(th_sprk_u[indx_u+jj] - th_cap);

	 coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

	 x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;
         a_sprk[ncnt] = h_sprk;

         //Spark at leading edge
         if (M_PI - co_angl - th_sprk_u[indx_u+jj] <= theta_sp[ii]/2) {
		 
            theta_fh = 0.5*(M_PI-co_angl-th_sprk_u[indx_u+jj]) + theta_sp[ii]/4;

            a_sprk[ncnt] = a_trk*sin(theta_fh);

            trk_a_ul[ii] = a_trk + h_sprk - a_sprk[ncnt];

            x_in = trk_a_ul[ii]*cos(M_PI - co_angl - theta_fh - th_cap);
            y_in = trk_a_ul[ii]*b_trk/a_trk*sin(M_PI-co_angl-theta_fh-th_cap);

	    coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

	    x_sprk[ncnt]+=x_cent;
	    y_sprk[ncnt]+=y_cent;
         }

         //Spark at trailing edge
         if (th_sprk_u[indx_u+jj] <= theta_sp[ii]/2 - co_angl) {

	    theta_fh = 0.5*(th_sprk_u[indx_u+jj] + theta_sp[ii]/2 + co_angl);

            a_sprk[ncnt] = a_trk*sin(theta_fh);

	    trk_a_ut[ii] = a_trk + h_sprk - a_sprk[ncnt];

            x_in = trk_a_ut[ii]*cos(theta_fh - th_cap - co_angl);
            y_in = trk_a_ut[ii]*b_trk/a_trk*sin(theta_fh - th_cap - co_angl);

	    coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

	    x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

	 ncnt++;
      }


      //lower half, anti-clockwise track
      for (jj=0;jj<N_dn[ii];jj++) {

         x_in = a_trk*cos(th_sprk_d[indx_d+jj] - th_cap);
         y_in = b_trk*sin(th_sprk_d[indx_d+jj] - th_cap);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;
         a_sprk[ncnt] = h_sprk;
      
	 //Spark at leading edge
         if (th_sprk_d[indx_d+jj] - M_PI + co_angl <= theta_sp[ii]/2) {

            theta_fl = 0.5*(th_sprk_d[indx_d+jj]-M_PI+co_angl)+theta_sp[ii]/4;

            a_sprk[ncnt] = 0.5*(2*a_trk*sin(theta_fl));

            trk_a_dl[ii] = a_trk + h_sprk - a_sprk[ncnt];

            x_in = trk_a_dl[ii]*cos(M_PI - co_angl + theta_fl - th_cap);
            y_in = trk_a_dl[ii]*b_trk/a_trk*sin(M_PI-co_angl+theta_fl-th_cap);

	    coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

            x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }
      
	 //Spark at trailing edge
         if (2*M_PI-th_sprk_d[indx_d+jj] - co_angl <= theta_sp[ii]/2) {

            theta_fl = 0.5*(2*M_PI-th_sprk_d[indx_d+jj]-co_angl)+theta_sp[ii]/4;

            a_sprk[ncnt] = 0.5*(2*a_trk*sin(theta_fl));

            trk_a_dt[ii] = a_trk + h_sprk - a_sprk[ncnt];

            x_in = trk_a_dt[ii]*cos(2*M_PI - theta_fl - th_cap - co_angl);
	    y_in = trk_a_dt[ii]*b_trk/a_trk*sin(2*M_PI-theta_fl-th_cap-co_angl);

	    coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

	    x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

         ncnt++;
      }


      //Plugging gap in the beginning
      if (th_sprk_d[indx_d] - th_sprk_u[indx_u] >= theta_sp[ii]) {

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(th_sprk_d[indx_d]-th_sprk_u[indx_u])) - 2.0*h_sprk);

	 e_trk = a_trk + h_sprk - a_sprk[ncnt];

         x_in = e_trk*cos(M_PI - th_cap - co_angl);
         y_in = e_trk*b_trk/a_trk*sin(M_PI - th_cap - co_angl);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;


         ncnt++;
      }


      //Plugging gap in the end
      theta_fh = th_sprk_u[indx_u+N_up[ii]-1];
      theta_fl = 2*M_PI - th_sprk_d[indx_d+N_dn[ii]-1];

      if (theta_fl+theta_fh >= theta_sp[ii] && 
		                 theta_fl+theta_fh <= 2*theta_sp[ii]) {

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(theta_fl+theta_fh))-2.0*h_sprk);
    
	 e_trk = a_trk + h_sprk - a_sprk[ncnt]; 

	 x_in = e_trk*cos(2*M_PI - th_cap - co_angl);
         y_in = e_trk*b_trk/a_trk*sin(2*M_PI - th_cap - co_angl);
 
	 coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+= x_cent;
         y_sprk[ncnt]+=y_cent;

         ncnt++;
      }

      indx_u+=trk_max;
      indx_d+=trk_max;

      a_out-=2.0*h_sprk;
      a_in-=2.0*h_sprk;

      b_out-=b_cap/a_cap*2.0*h_sprk;
      b_in-=b_cap/a_cap*2.0*h_sprk;
   }

   //Core Spark
   x_sprk[ncnt] = x_cent;
   y_sprk[ncnt] = y_cent;
   a_sprk[ncnt] = a_cap - N_trk*2.0*h_sprk - 1.5*h_drft;
   ncnt++;


   //Section of the polar cap within each spark
   ncell = 0;
   x_val = h_drft/2;

   while (x_val < 2*x_cent) {

      y_val = h_drft/2;

      while (y_val < 2*y_cent) {

         for (ii=0;ii<ncnt;ii++) {
            x_in = x_val - x_sprk[ii];
	    y_in = y_val - y_sprk[ii];

	    coortrans(x_in, y_in, th_cap, &x_trns, &y_trns);

	    b_sprk = a_sprk[ii]*b_cap/a_cap;

	    el_val = sqrt(x_trns*x_trns/a_sprk[ii]/a_sprk[ii] +
                            y_trns*y_trns/b_sprk/b_sprk);

            if (el_val < 1) {
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


void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out)
{
   *x_out =  x_in*cos(theta) + y_in*sin(theta);

   *y_out = -x_in*sin(theta) + y_in*cos(theta);

}


void delay(int number_of_seconds)
{
    // Converting time into milli_seconds
    int milli_seconds = 1000 * number_of_seconds;

    // Storing start time
    clock_t start_time = clock();

    // looping till required time is not achieved
    while (clock() < start_time + milli_seconds);
}
