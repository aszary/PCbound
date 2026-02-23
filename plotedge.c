// Estimating the separation of boundary from surrounding sparks for plotting
//
// To compile: gcc plotedge.c -o ~/bin/plotedge -lm
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       D_cap, h_sprk, h_drft;
   double       rad_in, rad_out;
   double       rad_lo, rad_hi, th_min, th_max;
   double      *rad_dist, *rad_sft, *th_sprk, *th_sft1, *th_sft2;
   double       theta_sp, theta_c;
   double      *x_sprk, *y_sprk, *x_sft, *y_sft;
   double       x_val, y_val, x_cent, y_cent, x_cel, y_cel;
   double       x_max, y_max, x_min, y_min;
   double       th_l, th_h, th_val, del_th;
   int          N_sprk, N_sprk_u, N_sprk_d;
   int          ii, jj, ncnt, nplt;


   //Initializing parameters

   D_cap = 30.0;  // polar cap diameter in meters.
   h_sprk = 5.2;  // spark diameter in meters.
   h_drft = 0.75;  // drift movement in meters.

   x_cent = D_cap/2;
   y_cent = D_cap/2;

   rad_out = D_cap/2;
   rad_in = rad_out - h_sprk;

   N_sprk = 3.3*(rad_out*rad_out - rad_in*rad_in)/(h_sprk*h_sprk);

   x_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   y_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   x_sft = (double *)calloc(N_sprk+1, sizeof(double));
   y_sft = (double *)calloc(N_sprk+1, sizeof(double));
   rad_dist = (double *)calloc(N_sprk+1, sizeof(double));
   th_sprk = (double *)calloc(N_sprk+1, sizeof(double));
   rad_sft = (double *)calloc(N_sprk+1, sizeof(double));
   th_sft1 = (double *)calloc(N_sprk+1, sizeof(double));
   th_sft2 = (double *)calloc(N_sprk+1, sizeof(double));


   theta_sp = 2*M_PI/N_sprk;

   theta_c = M_PI;

   N_sprk_u = 0;

   while (theta_c >= 0) {
      x_sprk[N_sprk_u] = x_cent + 0.5*(rad_out+rad_in)*cos(theta_c);
      y_sprk[N_sprk_u] = y_cent + 0.5*(rad_out+rad_in)*sin(theta_c);
      th_sprk[N_sprk_u] = theta_c;

      x_sft[N_sprk_u] = x_cent + 0.5*(rad_out+rad_in)*cos(theta_c) + h_drft;
      y_sft[N_sprk_u] = y_cent + 0.5*(rad_out+rad_in)*sin(theta_c);

      N_sprk_u++;
      theta_c-=theta_sp;
   }

   theta_c = M_PI;

   N_sprk_d = 0;

   while (theta_c <= 2*M_PI) {
      x_sprk[N_sprk_u+N_sprk_d] = x_cent + 0.5*(rad_out+rad_in)*cos(theta_c);
      y_sprk[N_sprk_u+N_sprk_d] = y_cent + 0.5*(rad_out+rad_in)*sin(theta_c);
      th_sprk[N_sprk_u+N_sprk_d] = theta_c;

      x_sft[N_sprk_u+N_sprk_d] = x_cent + 0.5*(rad_out+rad_in)*cos(theta_c) 
	                         + h_drft;
      y_sft[N_sprk_u+N_sprk_d] = y_cent + 0.5*(rad_out+rad_in)*sin(theta_c);

      N_sprk_d++;
      theta_c+=theta_sp;
   }


   //Writing the polar cap outline
   fpt = fopen ("pol_cap.dat", "w");

   th_val = 0.0;

   del_th = 2*M_PI/1000;

   while (th_val <= 2*M_PI) {

      x_val = x_cent + rad_out*cos(th_val); 
      y_val = y_cent + rad_out*sin(th_val);

      fprintf (fpt, "%lf  %lf\n", x_val, y_val);

      th_val+=del_th;
   }
   fclose(fpt);

   //Writing the Spark outline for upper track
   fpt = fopen ("spark_up.dat", "w");

   for (ii=0;ii<N_sprk_u-1;ii++) {

      th_val = 0;
      
      del_th = 2*M_PI/1000;

      while (th_val <= 2*M_PI) {

         x_val = x_sprk[ii] + (h_sprk/2)*cos(th_val);
         y_val = y_sprk[ii] + (h_sprk/2)*sin(th_val);

         fprintf (fpt, "%d  %lf  %lf  %lf\n", ii, x_val, y_val, x_val+h_drft);

         th_val+=del_th;
      }
      fprintf (fpt, "\n");
   }
   fclose(fpt);
   
 

   //Finding the distance of polar cap edge from spark centers
   fpt = fopen ("dist_pol.dat", "w");

   for (ii=0;ii<N_sprk_u-1;ii++) {

      th_val = th_sprk[ii];
      th_h = th_sprk[ii+1]; 

      rad_lo = 2*rad_out;
      rad_hi = 2*rad_out;
      th_min = th_val;
      th_max = th_val;

      while (th_val >= th_h) {

         x_val = x_cent + rad_out*cos(th_val);
         y_val = y_cent + rad_out*sin(th_val);

         for (jj=ii;jj<=ii+1;jj++) {

	    //Distance for edge spark
            x_cel = x_val - x_sprk[jj];
            y_cel = y_val - y_sprk[jj];

	    rad_dist[jj-ii] = sqrt(x_cel*x_cel + y_cel*y_cel);

	    //Distance for shifted spark
	    x_cel = x_val - x_sft[jj];
            y_cel = y_val - y_sft[jj];

	    rad_sft[jj-ii] = sqrt(x_cel*x_cel + y_cel*y_cel);
         }

         if (rad_dist[0]+rad_dist[1] < rad_lo) {
            rad_lo = rad_dist[0]+rad_dist[1];
	    th_min = th_val;
         }

         if (rad_sft[0]+rad_sft[1] < rad_hi) {
            rad_hi = rad_sft[0]+rad_sft[1];
            th_max = th_val;
         }

         th_val-=0.1*h_drft/(rad_out+rad_in);
      }
      x_min = x_cent + rad_out*cos(th_min);
      y_min = y_cent + rad_out*sin(th_min);

      x_max = x_cent + rad_out*cos(th_max);
      y_max = y_cent + rad_out*sin(th_max);

      fprintf(fpt, "%d  %lf %lf %lf %lf  %lf %lf %lf %lf %lf %lf %lf\n", ii, 0.5*(th_sprk[ii]+th_sprk[ii+1]), th_max-th_min, th_min, th_max, x_min, y_min, x_max, y_max, x_sprk[ii], y_sprk[ii], x_sprk[ii]+h_drft);

      th_sft1[ii] = th_max + theta_sp/2;
      th_sft2[ii] = th_max - theta_sp/2;
   }


   for (ii=N_sprk_u;ii<N_sprk_u+N_sprk_d-1;ii++) {

      th_val = th_sprk[ii];
      th_h = th_sprk[ii+1];

      rad_lo = 2*rad_out;
      rad_hi = 2*rad_out;
      th_min = th_val;
      th_max = th_val;

      while (th_val <= th_h) {

         x_val = x_cent + rad_out*cos(th_val);
         y_val = y_cent + rad_out*sin(th_val);

         for (jj=ii;jj<=ii+1;jj++) {

            //Distance for edge spark
            x_cel = x_val - x_sprk[jj];
            y_cel = y_val - y_sprk[jj];

            rad_dist[jj-ii] = sqrt(x_cel*x_cel + y_cel*y_cel);

            //Distance for shifted spark
            x_cel = x_val - x_sft[jj];
            y_cel = y_val - y_sft[jj];

            rad_sft[jj-ii] = sqrt(x_cel*x_cel + y_cel*y_cel);
         }

         if (rad_dist[0]+rad_dist[1] < rad_lo) {
            rad_lo = rad_dist[0]+rad_dist[1];
            th_min = th_val;
         }

         if (rad_sft[0]+rad_sft[1] < rad_hi) {
            rad_hi = rad_sft[0]+rad_sft[1];
            th_max = th_val;
         }

         th_val+=0.1*h_drft/(rad_out+rad_in);
      }
      

      fprintf(fpt, "%d  %lf %lf %lf %lf\n", ii, 0.5*(th_sprk[ii]+th_sprk[ii+1]), th_max-th_min, th_min, th_max);

      th_sft1[ii] = th_max - theta_sp/2;
      th_sft2[ii] = th_max + theta_sp/2;
   }
   fclose(fpt);


   fpt = fopen ("edge_sprk_sep.dat", "w");

   for (ii=1;ii<N_sprk_u-1;ii++) 
      fprintf (fpt, "%d  %lf %lf\n", ii, th_sprk[ii], th_sft1[ii]-th_sft2[ii-1]);

   
   for (ii=N_sprk_u+1;ii<N_sprk_u+N_sprk_d-1;ii++) 
      fprintf (fpt, "%d  %lf %lf\n", ii, th_sprk[ii], th_sft1[ii]-th_sft2[ii-1]);

   fclose(fpt);


   return 0;
}
