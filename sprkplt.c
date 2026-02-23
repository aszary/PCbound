// Reading the polar cap and LOS from file and plot spark configuration
//
//
// To compile: gcc sprkplt.c -o ~/bin/sprkplt -lcpgplot -lm
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cpgplot.h"
#include "usefunc.h"


# define ZERO_lf (0.000001) // round off zero error.
# define WORDSIZE (100)      // Size of filename.
# define LINELEN  (1000)      // Size of filename.

#define USAGE "\
   usage: sprkplt  inopen  trkfile  outname  ntime\n"


void sparkconfig (float *sprk_I, double *th_sprk_u, double *th_sprk_d, 
		     double *x_sprk, double *y_sprk, double *a_sprk, 
		     double *theta_sp, double h_sprk, double h_drft, 
		     double a_cap, double b_cap, double th_cap, double co_angl,
		     double x_cent, double y_cent, int N_trk, int *N_up, 
		     int *N_dn, int trk_max);

void plotcap (char *filename, float *sprk_I, float *x_los, float *y_los, 
		 double x_cent, double y_cent, double x_dist, double y_dist, 
		 double a_cap, double b_cap, double th_cap, int ncell_x, 
		 int ncell_y, int nlos, int nseq);

void unitvect (double *invect, double *transvect);

void rotvect (double *invect, double *rotaxis, double *transvect, 
		 double th_rot);

void crossprod (double *vect1, double *vect2, double *crosvect);

double dotprod (double *vect1, double *vect2);

int main(int argc, char *argv[])
{
   FILE        *fpt;
   char        *openfile, *trkfile, *outname;
   char        *readline, *readval[20];
   double       a_cap, b_cap, th_cap, co_angl, RS;
   double       x_cent, y_cent, z_cent;
   double       h_sprk, h_drft, a_sprk, b_sprk;
   double      *x_sprk, *y_sprk, *rad_sprk;
   double      *theta_sp, *th_sprk_u, *th_sprk_d;
   double      *outfield, *los_tr;
   double       inarr[3], rotaxis[3], trnsvect[3];
   double       max_phi, min_phi, phi_c;
   double       max_th, min_th, th_c;
   double       del_theta;
   double       a_in, a_out, a_trk, b_in, b_out;
   double       x_in, y_in, x_trns, y_trns;
   int         *N_up, *N_dn;
   int          Max_sprk, N_sprk, N_trk, N_cell, nlos;
   int          ncel, ncnt, ncap, nfield, ntime;
   int          ii, jj, indx_u, indx_d, trk_max;

   char        *label;
   float       *x_los, *y_los, *sprk_I;
   double       x_dist, y_dist;
   float        xmin, xmax, ymin, ymax;
   int          ncell_x, ncell_y;


   openfile = (char *)calloc(WORDSIZE, sizeof(char));
   trkfile = (char *)calloc(WORDSIZE, sizeof(char));
   outname = (char *)calloc(WORDSIZE, sizeof(char));
   readline = (char *)calloc(LINELEN, sizeof(char));
   label = (char *)calloc(WORDSIZE, sizeof(char));

   if (argc < 5)
      cmdlineError (USAGE);

   sprintf (openfile, "%s", argv[1]);
   sprintf (trkfile, "%s", argv[2]);
   sprintf (outname, "%s", argv[3]);
   ntime = atoi(argv[4]);


   //Reading polar cap configuration
   if ((fpt = fopen(openfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", openfile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);

   nfield = atoi(readval[1]);

   fprintf(stderr,"nfield %d\n", nfield);

   fgets(readline, LINELEN, fpt);

   splitstr (readline, readval, " ", 1);

   x_cent = atof(readval[1]);
   y_cent = atof(readval[3]);
   a_cap = atof(readval[5]);
   b_cap = atof(readval[7]);
   th_cap = M_PI - atof(readval[11])/180*M_PI;

   outfield = (double *)calloc(6*nfield+1, sizeof(double));

   for (ii=0;ii<nfield;ii++) {
      fgets(readline, LINELEN, fpt);
      splitstr (readline, readval, " ", 1);

      outfield[6*ii] = atof(readval[0]);
      outfield[6*ii+1] = atof(readval[1]);
      outfield[6*ii+2] = atof(readval[2]);
      outfield[6*ii+3] = atof(readval[6]);
      outfield[6*ii+4] = atof(readval[7]);
      outfield[6*ii+5] = atof(readval[8]);
   }
   fclose(fpt);

   //Finding minimum and maximum theta and phi
   max_th = outfield[4];
   min_th = outfield[4];

   max_phi = outfield[5];
   min_phi = outfield[5];

   for (ii=1;ii<nfield;ii++) {

      if (outfield[6*ii+4] > max_th)
         max_th = outfield[6*ii+4];

      if (outfield[6*ii+4] < min_th)
         min_th = outfield[6*ii+4];

      if (outfield[6*ii+5] > max_phi)
         max_phi = outfield[6*ii+5];

      if (outfield[6*ii+5] < min_phi)
         min_phi = outfield[6*ii+5];
   }

   th_c = 0.5*(max_th+min_th);
   phi_c = 0.5*(max_phi+min_phi);


   //Reading the LOS track on the surface
   if ((fpt = fopen(trkfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", trkfile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);
   nlos = atoi(readval[1]);

   los_tr = (double *)calloc(6*nlos+1, sizeof(double));
   x_los = (float *)calloc(nlos+1, sizeof(float));
   y_los = (float *)calloc(nlos+1, sizeof(float));

   fprintf(stderr,"nlos %d\n", nlos);

   for (ii=0;ii<nlos;ii++) {
      fgets(readline, LINELEN, fpt);
      splitstr (readline, readval, " ", 1);

      los_tr[6*ii] = atof(readval[0]);
      los_tr[6*ii+1] = atof(readval[1]);
      los_tr[6*ii+2] = atof(readval[2]);
      los_tr[6*ii+3] = atof(readval[3]);
      los_tr[6*ii+4] = atof(readval[4]);
      los_tr[6*ii+5] = atof(readval[5]);
   }
   fclose (fpt);

   RS = 10000.0; //Radius of neutron star

   z_cent = RS*cos(th_c);

   inarr[0] = sin(th_c)*cos(max_phi) - sin(th_c)*cos(min_phi);
   inarr[1] = sin(th_c)*sin(max_phi) - sin(th_c)*sin(min_phi);
   inarr[2] = 0.0;

   unitvect (inarr, rotaxis);

   for (ii=0;ii<nlos;ii++) {
      inarr[0] = los_tr[6*ii+3] - x_cent;
      inarr[1] = los_tr[6*ii+4] - y_cent;
      inarr[2] = los_tr[6*ii+5] - z_cent;       

      rotvect (inarr, rotaxis, trnsvect, -th_c);

      x_los[ii] = (float)(trnsvect[0]);
      y_los[ii] = (float)(trnsvect[1]);
   }


   //Initializing parameters

   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.
   h_sprk = a_cap/2.7;
   h_drft = h_sprk/25;

   a_sprk = h_sprk/2;
   b_sprk = a_sprk*b_cap/a_cap;

   co_angl = phi_c-M_PI/2;


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

   x_dist = sqrt(pow(a_cap*cos(th_cap),2) + pow(b_cap*sin(th_cap),2));
   y_dist = sqrt(pow(a_cap*sin(th_cap),2) + pow(b_cap*cos(th_cap),2));

   ncell_x = 2*x_dist/h_drft + 1;
   ncell_y = 2*y_dist/h_drft + 1;

   fprintf (stderr, "%d  %d  %d  %d\n", N_trk, Max_sprk, ncell_x, ncell_y); 


   //Initializing the spark configuration
   x_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   y_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   rad_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   th_sprk_u = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));
   th_sprk_d = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));
   
   sprk_I = (float *)calloc(ncell_x*ncell_y+1, sizeof(float));
  
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


   // Plotting the series of images of spark evolution
   for (ii=0;ii<nlos;ii++) {

      x_los[ii] = x_los[ii]*(float)(ncell_x)/(float)(2*x_dist) + (float)(ncell_x/2);
      y_los[ii] = y_los[ii]*(float)(ncell_y)/(float)(2*y_dist) + (float)(ncell_y/2);
   }

   for (ii=0;ii<ntime;ii++) {

      fprintf (stderr,"time %d/%d\r", ii, ntime);

      //Estimating the Spark configuration 
      sparkconfig (sprk_I, th_sprk_u, th_sprk_d, x_sprk, y_sprk, rad_sprk, 
                      theta_sp, h_sprk, h_drft, a_cap, b_cap, th_cap, co_angl, 
                      x_cent, y_cent, N_trk, N_up, N_dn, trk_max);

/*
   fpt = fopen ("tst_sprk.dat","w");

   y_in = y_cent + h_drft/2;
   for (ii=0;ii<ncell_y;ii++) {
      x_in = x_cent + h_drft/2;

      for (jj=0;jj<ncell_x;jj++) {
         fprintf (fpt, "%d  %d %f\n", jj, ii, sprk_I[ii*ncell_x+jj]);              
	 x_in+=h_drft;
      }
      y_in+=h_drft;
   }
*/

      //Plotting the sparks on the polar cap
      plotcap (outname, sprk_I, x_los, y_los, x_cent, y_cent, x_dist, y_dist, 
                 a_cap, b_cap, th_cap, ncell_x, ncell_y, nlos, ii+1);

 
      //Initializing spark location for next time
      indx_u = 0;
      indx_d = 0;

      a_out = a_cap;
      a_in = a_out - 2*a_sprk;

      a_trk = 0.5*(a_out+a_in);

      for (jj=0;jj<N_trk;jj++) {

         del_theta = h_drft/(2*a_trk);

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

         a_trk-=a_sprk;
      }
   }


   return 0;
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


void crossprod (double *vect1, double *vect2, double *crosvect)
{
   crosvect[0] = vect1[1]*vect2[2] - vect1[2]*vect2[1];

   crosvect[1] = vect1[2]*vect2[0] - vect1[0]*vect2[2];

   crosvect[2] = vect1[0]*vect2[1] - vect1[1]*vect2[0];

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


void sparkconfig (float *sprk_I, double *th_sprk_u, double *th_sprk_d, 
		     double *x_sprk, double *y_sprk, double *a_sprk, 
		     double *theta_sp, double h_sprk, double h_drft, 
		     double a_cap, double b_cap, double th_cap, double co_angl,
		     double x_cent, double y_cent, int N_trk, int *N_up, 
		     int *N_dn, int trk_max)
{
   double      trk_a_ul[N_trk], trk_b_ul[N_trk];
   double      trk_a_ut[N_trk], trk_b_ut[N_trk];
   double      trk_a_dl[N_trk], trk_b_dl[N_trk];
   double      trk_a_dt[N_trk], trk_b_dt[N_trk];
   double      a_in, a_out, b_in, b_out, e_trk, b_sprk;
   double      x_val, y_val, x_in, y_in, x_trns, y_trns, el_val;
   double      theta_c, theta_fl, theta_fh;
   double      a_trk, b_trk, rad_fh, rad_fl;
   double      x_dist, y_dist;
   float       sig_x, sig_y;
   int         indx_u, indx_d, indx_us, ncell_x, ncell_y;
   int         ii, jj, kk, ncnt, ncell, nsprk_l, nsprk_h, if_sprk;


   a_out = a_cap;
   a_in = a_out - h_sprk;

   b_out = b_cap;
   b_in = b_out - b_cap/a_cap*h_sprk; 

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
         a_sprk[ncnt] = h_sprk/2;

         //Spark at leading edge
         if (M_PI - co_angl - th_sprk_u[indx_u+jj] <= theta_sp[ii]/2) {
		 
            theta_fh = 0.5*(M_PI-co_angl-th_sprk_u[indx_u+jj]) + theta_sp[ii]/4;

            a_sprk[ncnt] = a_trk*sin(theta_fh);

            trk_a_ul[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

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

	    trk_a_ut[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

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
         a_sprk[ncnt] = h_sprk/2;
      
	 //Spark at leading edge
         if (th_sprk_d[indx_d+jj] - M_PI + co_angl <= theta_sp[ii]/2) {

            theta_fl = 0.5*(th_sprk_d[indx_d+jj]-M_PI+co_angl)+theta_sp[ii]/4;

            a_sprk[ncnt] = 0.5*(2*a_trk*sin(theta_fl));

            trk_a_dl[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

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

            trk_a_dt[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

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

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(th_sprk_d[indx_d]-th_sprk_u[indx_u])) - h_sprk);

	 e_trk = a_trk + h_sprk/2 - a_sprk[ncnt];

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

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(theta_fl+theta_fh)) - h_sprk);
    
	 e_trk = a_trk + h_sprk/2 - a_sprk[ncnt]; 

	 x_in = e_trk*cos(2*M_PI - th_cap - co_angl);
         y_in = e_trk*b_trk/a_trk*sin(2*M_PI - th_cap - co_angl);
 
	 coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+= x_cent;
         y_sprk[ncnt]+=y_cent;

         ncnt++;
      }

      indx_u+=trk_max;
      indx_d+=trk_max;

      a_out-=h_sprk;
      a_in-=h_sprk;

      b_out-=b_cap/a_cap*h_sprk;
      b_in-=b_cap/a_cap*h_sprk;
   }

   //Core Spark
   x_sprk[ncnt] = x_cent;
   y_sprk[ncnt] = y_cent;
   a_sprk[ncnt] = a_cap - N_trk*h_sprk - 1.5*h_drft;
   ncnt++;


   //Section of the polar cap within each spark
   x_dist = sqrt(pow(a_cap*cos(th_cap),2) + pow(b_cap*sin(th_cap),2));
   y_dist = sqrt(pow(a_cap*sin(th_cap),2) + pow(b_cap*cos(th_cap),2));

   ncell_x = 2*x_dist/h_drft + 1;
   ncell_y = 2*y_dist/h_drft + 1;

   y_val = y_cent - y_dist + h_drft/2;

   for (ii=0;ii<ncell_y;ii++) {

      x_val = x_cent - x_dist + h_drft/2;

      for (jj=0;jj<ncell_x;jj++) {
	      
         sprk_I[ii*ncell_x+jj] = 0.0;

	 kk = 0;

	 if_sprk = 0;

	 while (if_sprk < 1) {

	    if (kk == ncnt-1) 
	       if_sprk = 1;

            x_in = x_val - x_sprk[kk];
	    y_in = y_val - y_sprk[kk];

	    coortrans(x_in, y_in, th_cap, &x_trns, &y_trns);

	    b_sprk = a_sprk[kk]*b_cap/a_cap;

	    el_val = sqrt(x_trns*x_trns/a_sprk[kk]/a_sprk[kk] +
                            y_trns*y_trns/b_sprk/b_sprk);

            if (el_val < 1) {
               sig_x = (float)(a_sprk[kk]/1.3);
	       sig_y = (float)(b_sprk/1.3);

               sprk_I[ii*ncell_x+jj] = (float)(exp(-x_trns*x_trns/(2.0*sig_x*sig_x) - y_trns*y_trns/(2.0*sig_y*sig_y)));

	       if_sprk = 1;
            }
	    kk++;
	 }
         x_val+=h_drft;
      }
      y_val+=h_drft;
   }
}


void plotcap (char *filename, float *sprk_I, float *x_los, float *y_los, 
		 double x_cent, double y_cent, double x_dist, double y_dist, 
		 double a_cap, double b_cap, double th_cap, int ncell_x, 
		 int ncell_y, int nlos, int nseq)
{
   FILE              *fpt;
   char               psfile[10000];
   float             *plotstk, *x_elips, *y_elips;
   double             x_in, y_in, x_trns, y_trns;
   double             theta_c, del_th, a_out, b_out;
   int                ii, jj, indx, N_elips;

   // pgplot inputs
   float              tr[] = { 0, 1, 0, 0, 0, 1 };  // Identity Mapping
   float              HL[] = { 0.0, 0.2, 0.4, 0.6, 1.0 };
   float              HR[] = { 0.0, 0.5, 1.0, 1.0, 1.0 };
   float              HG[] = { 0.0, 0.0, 0.5, 1.0, 1.0 };
   float              HB[] = { 0.0, 0.0, 0.0, 0.3, 1.0 };
   float   RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
   float   RR[] = {0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
   float   RG[] = {0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
   float   RB[] = {0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};
   float   WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0 };
   float   WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
   float   WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
   float   WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};

   float              bri, cont, nc;
   float              xleft, xright, ybot, ytop;
   float              xmin, xmax, ymin, ymax, hi, lo;


   sprintf(psfile, "%s_sprkcap_%03d.ps/CPS", filename, nseq);

   plotstk  = (float *)calloc(ncell_x*ncell_y, sizeof(float));

   hi = sprk_I[0];
   lo = sprk_I[0];

   for (ii=0;ii<ncell_y;ii++) {
      for (jj=0;jj<ncell_x;jj++) {
         indx = ii*ncell_x+jj;

	 plotstk[indx] = sprk_I[indx];

         if (plotstk[indx] > hi)
            hi = plotstk[indx];

         if (plotstk[indx] < lo)
            lo = plotstk[indx];
      }
   }
   hi*=1.05;

   //Plotting the polar cap
   cpgbeg(0, psfile, 1, 1);

   xleft =  0.2;
   xright = 0.8;
   ybot = 0.2;
   ytop = 0.8;

   cpgsvp(xleft, xright, ybot, ytop);


   xmin = 0.0;
   xmax = ncell_x;

   ymin = 0.0;
   ymax = ncell_y;


//   cpgenv(xmin, xmax, ymin, ymax, 0, 0);
   cpgswin(xmin, xmax, ymin, ymax);

   cont = 0.95;
   bri = 0.5;
   nc = 5;

   cpgctab (HL, HR, HG, HB, nc, cont, bri);
//   cpgctab (RL, RR, RG, RB, nc, cont, bri);
//   cpgctab (WL, WR, WG, WB, nc, cont, bri);

   cpgimag (sprk_I, ncell_x, ncell_y, 1, ncell_x, 1, ncell_y, hi,
             lo, tr);

   cpgbox ("BC", 0.0, 0, "BC", 0.0, 0);

   cpgsch(1.4);

   cpglab("X (meters)","Y (meters)", "Sparking Evolution in Polar Cap");

   cpgaxis("N", xmin, ymin, xmin, ymax, -y_dist, y_dist, 0.0, 0, 
             0.5, 0.5, 0.5, -1.0, 0.0);

   cpgaxis("N", xmin, ymin, xmax, ymin, -x_dist, x_dist, 0.0, 0, 
              0.5, 0.5, 0.5, 1.0, 0.0);

   //The outline of the elliptical polar cap
   N_elips = 100;

   x_elips = (float *)calloc(N_elips+1, sizeof(float));
   y_elips = (float *)calloc(N_elips+1, sizeof(float));

   theta_c = 0;

   del_th = 2*M_PI/N_elips;

   ii = 0;

   a_out = a_cap*(double)(ncell_x)/(2*x_dist);
   b_out = b_cap/a_cap*a_out;

   while (theta_c < 2*M_PI) {
      
      x_in = a_out*cos(theta_c);
      y_in = b_out*sin(theta_c);

      coortrans (x_in, y_in, -th_cap, &x_trns, &y_trns);

      x_elips[ii] = (float)(x_trns) + xmax/2;
      y_elips[ii] = (float)(y_trns) + ymax/2;
      ii++;

      theta_c+=del_th;
   }
   x_elips[ii] = x_elips[0];
   y_elips[ii] = y_elips[0];

   cpgsci(4);
   cpgslw(3);

   cpgline(N_elips+1, x_elips, y_elips);


   //The LOS track on the polar cap

   cpgline(nlos, x_los, y_los);

   cpgclos();
}
