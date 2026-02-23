// Estimating the spark motion distribution in two dimension 
//
// To compile: gcc twoDgrid.c -o ~/bin/twoDgrid -lm

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


# define C_LIGHT (299792458) // The speed of light.
# define Hz2MHz (1000000.0)  // frequency conversion.
# define ZERO_lf (0.000001) // round off zero error.
# define LINELEN (10000)     // Size of line.
# define WORDSIZE (100)      // Size of filename.


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       D_cap, h_sprk, h_drft;
   double       x_val, y_val, rad_val, wt, tot_wt;
   double       x_cent, y_cent;
   double      *x_sprk, *y_sprk, x_cel, y_cel;
   double       theta_sp, theta_c, rad_in, rad_out;
   int         *N_sprk, tot_sprk, N_trk;
   int          ii, jj, ncel, ncnt;



   //Initializing parameters

   D_cap = 15.0;  // polar cap diameter in meters.
   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.

   fpt = fopen ("polarcap2d.dat", "w");

   x_val = h_drft/2;

   x_cent = D_cap/2;
   y_cent = D_cap/2;

   tot_wt = 0.0;
   ncel = 0;

   while (x_val < D_cap) {
   
      y_val = h_drft/2;

      while (y_val < D_cap) {

         rad_val = sqrt(pow((x_val-x_cent),2.0) + pow((y_val-y_cent),2.0));

	 if (rad_val < D_cap/2) { 

            ncel++;
            wt = 1.0;

            if (rad_val > D_cap/2-h_sprk/2)
               wt = 2.0;

	       tot_wt+=wt;

            fprintf (fpt, "%lf  %lf  %lf\n", x_val, y_val, wt);
	 }
	
         y_val+=h_drft;
      }

      x_val+=h_drft;
   }

   fclose(fpt);

   fprintf(stderr, "%d  %lf\n", ncel, tot_wt);


   // The packing of sparks
   N_trk = D_cap/(2.0*h_sprk);

   N_sprk = (int *)calloc(N_trk+1, sizeof(int));

   rad_out = D_cap/2;	   
   rad_in = rad_out - h_sprk;

   tot_sprk = 0;

   for (ii=0;ii<N_trk;ii++) {
      
      N_sprk[ii] = 3.3*(rad_out*rad_out - rad_in*rad_in)/(h_sprk*h_sprk);

      tot_sprk+=N_sprk[ii];

      rad_out-=h_sprk;
      rad_in-=h_sprk;
   }

   fprintf (stderr, "N_trk %d  tot_sprk %d\n", N_trk, tot_sprk);

   x_sprk = (double *)calloc(tot_sprk+1, sizeof(double));
   y_sprk = (double *)calloc(tot_sprk+1, sizeof(double));
   
   //Finding the central locations of sparks
   rad_out = D_cap/2;	   
   rad_in = rad_out - h_sprk;

   ncnt = 0;
   for (ii=0;ii<N_trk;ii++) {
     
      theta_sp = 2*M_PI/N_sprk[ii];

      theta_c = theta_sp/2;

      for (jj=0;jj<N_sprk[ii];jj++) {
         x_sprk[ncnt] = x_cent + 0.5*(rad_out+rad_in)*cos(theta_c);
         y_sprk[ncnt] = y_cent + 0.5*(rad_out+rad_in)*sin(theta_c);

         theta_c+=theta_sp;
	 ncnt++;
      }
      rad_out-=h_sprk;
      rad_in-=h_sprk;
   }

   //location of the polar cap within each spark
   fpt = fopen ("spark2d.dat", "w");

   x_val = h_drft/2;

   rad_in = D_cap/2-N_trk*h_sprk-1.5*h_drft;
   
   while (x_val < D_cap) {

      y_val = h_drft/2;

      while (y_val < D_cap) {

         for (ii=0;ii<tot_sprk;ii++) {
            x_cel = x_val - x_sprk[ii];
	    y_cel = y_val - y_sprk[ii];

            rad_val = sqrt(x_cel*x_cel + y_cel*y_cel);

            if (rad_val < 0.5*h_sprk) 
               fprintf (fpt, "%lf  %lf  %d\n", x_val, y_val, ii);
         }

	 //Core spark
	 x_cel = x_val - D_cap/2;
	 y_cel = y_val - D_cap/2;

         rad_val = sqrt(x_cel*x_cel + y_cel*y_cel);

	 if (rad_val < rad_in) 
            fprintf (fpt, "%lf  %lf  -1\n", x_val, y_val);


         y_val+=h_drft;
      }

      x_val+=h_drft;
   }

   fclose(fpt);


   return 0;
}
