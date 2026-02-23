//  To read the open field boundary and simulate single pulses with drifting
//  from lagging behind corotation + boundary effect
//
//  New addition : also carry out LRFS analysis on the single pulses.
//
//  To compile: gcc PCdrift.c -o ~/bin/PCdrift -lfftw3 -lcpgplot -lm
//
   

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <fftw3.h>
#include "cpgplot.h"
#include "PCdrift.h"
#include "plotfunc.h"
#include "usefunc.h"
#include "sparkfunc.h"


# define WORDSIZE (100) // The size of file names.
# define LINELEN (1000) // The word size of line in files.
# define ZERO_lf (0.000001) // The floating point infinitesimal.

#define USAGE "\
 usage: PCdrift inopen trkfile outfile alpha beta P3 npulse\n"


int main(int argc, char *argv[])
{ 
   FILE        *fpt;
   char        *openfile, *trkfile, *outname, *snglfile;
   char        *readline, *readval[20];
   double      *th, *pulstk, *fftout;
   double      *outfield, *los_tr;
   double       alpha, beta, period, P3, int_time;
   double       RS, a_cap, b_cap, th_cap, co_angl;
   double       max_th, min_th, max_phi, min_phi;
   double       theta, th_c, phi_c;
   double       x_cent, y_cent;
   double       phi_i, phi_f, phi, del_phi, Amp;
   double       phase, samp_phs, bin_phs, sigma;
   int         *addnum;
   int          nfield, ntrain, npulse, nprd, num, nlos; 
   int          nbin, sbin, bwin, ewin;
   int          ii, jj, tt, pp, nok, nbad, indx, st_indx, if_bound;
   long         seed;

   //Sparking variables
   double       h_sprk, h_drft, a_sprk, b_sprk;
   double      *x_sprk, *y_sprk, *rad_sprk;
   double      *theta_sp, *th_sprk_u, *th_sprk_d;
   double       theta_c, del_theta;
   double       a_in, a_out, b_in, b_out;
   double       x_out, y_out, sprk_amp;
   int         *N_up, *N_dn;
   int          Max_sprk, N_sprk, N_trk, nsprk;
   int          indx_u, indx_d, trk_max;


   openfile = (char *)calloc(WORDSIZE, sizeof(char));
   trkfile = (char *)calloc(WORDSIZE, sizeof(char));
   outname = (char *)calloc(WORDSIZE, sizeof(char));
   readline = (char *)calloc(LINELEN, sizeof(char));

   if (argc < 8)
      cmdlineError (USAGE);

   sprintf (openfile, "%s", argv[1]);
   sprintf (trkfile, "%s", argv[2]);
   sprintf (outname, "%s", argv[3]);
   alpha = atof(argv[4])/180.0*M_PI;
   beta = atof(argv[5])/180.0*M_PI;
   P3 = atof(argv[6]);
   npulse = atoi(argv[7]);

   period = 1.0;
   int_time = 0.000488;

   fprintf(stderr, "alpha %lf beta %lf period %lf int_time %lf npulse %d\n", alpha, beta, period, int_time, npulse);


//*** Reading the polar cap and line of sight from files *** //

   //Reading the open field boundary
   if ((fpt = fopen(openfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", openfile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);

   nfield = atoi(readval[1]);

   fprintf(stderr,"nfield %d\n", nfield);

   //Reading the elliptical fit to polar cap
   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);

   x_cent = atof(readval[1]);
   y_cent = atof(readval[3]);
   a_cap = atof(readval[5]);
   b_cap = atof(readval[7]);
   th_cap = M_PI - atof(readval[11])/180*M_PI;

   //Reading the polar cap boundary values
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

   fprintf (stderr, "x_cent %lf  y_cent %lf  a_cap %lf  b_cap %lf  th_cap %lf  th_c %lf  phi_c %lf\n", x_cent, y_cent, a_cap, b_cap, th_cap, th_c, phi_c);


   // Reading the line of sight track on the surface.
   if ((fpt = fopen(trkfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", trkfile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);
   nlos = atoi(readval[1]);
  
   los_tr = (double *)calloc(3*nlos+1, sizeof(double)); 

   fprintf(stderr,"nlos %d\n", nlos);
 
   for (ii=0;ii<nlos;ii++) {
      fgets(readline, LINELEN, fpt);
      splitstr (readline, readval, " ", 1);

      los_tr[3*ii] = atof(readval[0]);
      los_tr[3*ii+1] = atof(readval[1]);
      los_tr[3*ii+2] = atof(readval[2]);
   }
   fclose (fpt);

   phi_i = los_tr[0];
   phi_f = los_tr[3*(nlos-1)];

   if (fabs(phi_i) < ZERO_lf && fabs(phi_f) < ZERO_lf) {
      phi_i = -M_PI;
      phi_f = +M_PI;
   }

   fprintf(stderr,"%lf  %lf  %lf\n", phi_i*180/M_PI, phi_f*180/M_PI, (alpha+beta)*180/M_PI);


// *** Intitalizing sparking configuration in the 2D - polar cap *** //

   //Initializing parameters
   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.
   h_sprk = a_cap/2.7;
   h_drft = h_sprk/25;

   a_sprk = h_sprk/2;
   b_sprk = a_sprk*b_cap/a_cap;

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

   trk_max = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk)/2 + 1;

   for (ii=0;ii<N_trk;ii++) {

      N_sprk = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk);

      theta_sp[ii] = 2*M_PI/N_sprk;

      Max_sprk+=N_sprk;

      fprintf (stderr, "%d  %lf  %d\n", ii+1, theta_sp[ii], N_sprk);

      a_out-=2*a_sprk;
      a_in-=2*a_sprk;

      b_out-=2*b_sprk;
      b_in-=2*b_sprk;
   }
   fprintf (stderr, "%d  %d\n", N_trk, Max_sprk);

   //Initializing the spark configuration
   x_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   y_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   rad_sprk = (double *)calloc(Max_sprk+10, sizeof(double));
   th_sprk_u = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));
   th_sprk_d = (double *)calloc(2*trk_max*N_trk+1, sizeof(double));

   indx_u = 0;
   indx_d = 0;

   co_angl = phi_c-M_PI/2; 

   fpt = fopen("tst.out", "w");

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


// *** Generating single pulses with transition to polar cap *** //

   seed = -(long)time(NULL);
   sigma = 0.05;             //20-sigma single pulse peak.

   nbin = period/int_time;

   pulstk = (double *)calloc(nbin*npulse, sizeof(double));
   fftout = (double *)calloc(nbin*npulse, sizeof(double));
   addnum = (int *)calloc(nbin*npulse, sizeof(int));
   snglfile = (char *)calloc(WORDSIZE, sizeof(char));


   for (ii=0;ii<nbin*npulse;ii++) {
      pulstk[ii] = 0.0;
      addnum[ii] = 0;
   }

   samp_phs = int_time/period;
   bin_phs  = 1.0/(double)(nbin);

   sprintf (snglfile, "%s_ascii.dat", outname);

   fpt = fopen("tst.out", "w");

   nprd = 0;
   num = 0;

   //Setting up the initial sparking system
   sparkconfig (th_sprk_u, th_sprk_d, x_sprk, y_sprk, rad_sprk, theta_sp, 
		  h_sprk, h_drft, a_cap, b_cap, th_cap, co_angl, x_cent, 
		  y_cent, N_trk, N_up, N_dn, trk_max, &nsprk);

   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   while (nprd < npulse) {
      phi = (double)(num)*int_time/period*M_PI*2.0 - M_PI;

      while (phi > M_PI)
         phi-=2.0*M_PI;

      Amp = noise(sigma, &seed);

      //Transiting to the polar cap along open field lines
      if (phi > phi_i && phi < phi_f) {

         pospolcap (los_tr, phi, th_c, phi_c, max_phi, min_phi, &x_out, &y_out,
                      nlos);

         sprk_amp = sparkamp (x_sprk, y_sprk, rad_sprk, a_cap, b_cap, th_cap, 
			        x_out, y_out, nsprk);

	 Amp+=sprk_amp;

//         fprintf (fpt, "%d  %lf  %lf\n", nprd, x_out, y_out);
      }


      //Writing out the pulse stack
      phase = ((double)(num)*int_time/period) - (double)(nprd);

      sbin = (int)(phase/bin_phs);

      pulstk[nprd*nbin+sbin]+=Amp;
      addnum[nprd*nbin+sbin]++;

//      fprintf (fpt, "%d %lf %lf\n", nprd, phase, Amp);

      num++;


      //Initilazing for next period
      if ((double)(num)*int_time > (double)(nprd+1)*period) {
         nprd++;
         fprintf (stderr,"period %d/%d\r", nprd, npulse);

         //Evolving the sparking configuration
         for (jj=0;jj<nsprk;jj++) 
            fprintf (fpt, "%d %lf %lf %lf\n", nprd, x_sprk[jj], y_sprk[jj], rad_sprk[jj]);

         indx_u = 0;
         indx_d = 0;

         for (jj=0;jj<N_trk;jj++) {

            del_theta = theta_sp[jj]/P3;

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

            while (th_sprk_d[indx_d+N_dn[jj]-1]<=2*M_PI-theta_sp[jj]-co_angl) {

               th_sprk_d[indx_d+N_dn[jj]] = th_sprk_d[indx_d+N_dn[jj]-1]+theta_sp[jj];
               N_dn[jj]++;
            }
            indx_d+=trk_max;
         }

         sparkconfig (th_sprk_u, th_sprk_d, x_sprk, y_sprk, rad_sprk, theta_sp,
                        h_sprk, h_drft, a_cap, b_cap, th_cap, co_angl, x_cent, 
			y_cent, N_trk, N_up, N_dn, trk_max, &nsprk);

      }
   }
   fclose (fpt);


   fpt = fopen(snglfile, "w");

   for (ii=0;ii<npulse;ii++) {
      for (jj=0;jj<nbin;jj++) {
         if (addnum[ii*nbin+jj] > 0) {
            pulstk[ii*nbin+jj]/=addnum[ii*nbin+jj];
	    fprintf (fpt, "%d  %d %lf\n", ii, jj, pulstk[ii*nbin+jj]);
	 }
      }
   }
   fclose (fpt);

   // Branching out to estimating LRFS
   bwin = (phi_i+M_PI)/2.0/M_PI*nbin;
   ewin = (phi_f+M_PI)/2.0/M_PI*nbin;

   pulsfft (pulstk, fftout, npulse, nbin, bwin, ewin);


   // Plotting the LRFS
   plotlrfs (outname, pulstk, fftout, npulse, nbin, bwin, ewin);


   // Plotting the single pulses
   plotsngl (outname, pulstk, alpha, beta, npulse, nbin, bwin, ewin);


   return 0;
}
