// Estimating the separation of boundary from surrounding sparks
//
// To compile: gcc edgEllipse.c -o ~/bin/edgEllipse -lcpgplot -lm
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cpgplot.h"
#include <time.h>


#define USAGE "\
   usage: edgEllipse  ntime  incl_angle (deg)  maj_axis (m) ellip (a/b) lag_angl (deg)\n"


void edgedist (float *x_plot, float *y_plot, double *th_sprk_u, 
		  double *th_sprk_d, double *edg_sft_u, double *edg_sft_d, 
		  double a_cap, double b_cap, double th_cap, double th_lag, 
		  double h_sprk, double h_drft, double x_cent, double y_cent, 
		  int N_sprk_u, int N_sprk_d, int *nplt);


void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out);


void delay(int number_of_seconds);


void cmdlineError (char *message);


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       a_cap, b_cap, th_cap;
   double       a_sprk, b_sprk, h_sprk, h_drft;
   double       a_val, b_val, a_out, a_in, b_out, b_in;
   double      *th_sprk_u, *th_sprk_d, *edg_sft_u, *edg_sft_d, *th_buf;
   double       theta_sp, theta_c, th_lag;
   double       x_val, y_val, x_cent, y_cent;
   double       x_cel, y_cel, x_in, y_in;
   int          N_sprk, N_sprk_u, N_sprk_d, sp_cnt;
   int          ii, jj, ncnt, ntime;

   char        *label;
   float       *x_plot, *y_plot;
   float        xmin, xmax, ymin, ymax, x_diff, y_diff;
   float        xleft, xright, ybot, ytop;
   int          del_sec, nplt;


   if (argc < 6)
      cmdlineError (USAGE);

   ntime = atoi(argv[1]);
   th_cap = atof(argv[2])/180*M_PI;
   a_cap = atof(argv[3]);
   b_cap = a_cap*atof(argv[4]);
   th_lag = atof(argv[5])/180*M_PI;

   label = (char *)calloc(100, sizeof(char));
   

   //Initializing parameters
   h_sprk = 5.2;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.

   a_sprk = h_sprk/2;
   b_sprk = a_sprk*b_cap/a_cap;

   //Finding centre of rotated ellipse
   x_cent = sqrt(pow(a_cap*cos(th_cap),2) + pow(b_cap*sin(th_cap),2));
   y_cent = sqrt(pow(a_cap*sin(th_cap),2) + pow(b_cap*cos(th_cap),2));

   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   b_out = b_cap;
   b_in = b_out - 2*b_sprk;

   N_sprk = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk);

   th_sprk_u = (double *)calloc(N_sprk+1, sizeof(double));
   th_sprk_d = (double *)calloc(N_sprk+1, sizeof(double));
   edg_sft_u = (double *)calloc(N_sprk+1, sizeof(double));
   edg_sft_d = (double *)calloc(N_sprk+1, sizeof(double));
   th_buf = (double *)calloc(N_sprk+1, sizeof(double));
   x_plot = (float *)calloc(N_sprk+1, sizeof(float));
   y_plot = (float *)calloc(N_sprk+1, sizeof(float));


   //Initializing the 
   theta_sp = 2*M_PI/N_sprk;

   theta_c = M_PI - theta_sp/4;

   N_sprk_u = 0;

   while (theta_c >= 0) {

      th_sprk_u[N_sprk_u++] = theta_c;

      theta_c-=theta_sp;
   }

   theta_c = M_PI + theta_c/4;

   N_sprk_d = 0;

   while (theta_c <= 2*M_PI) {

      th_sprk_d[N_sprk_d++] = theta_c;

      theta_c+=theta_sp;
   }


   //Finding the distance of polar cap edge from spark centers
   cpgbeg(0, "/xs", 1, 1);

   xleft =  0.1;
   xright = 0.9;
   ybot = 0.2;
   ytop = 0.8;

   del_sec = 100;

   fpt = fopen ("surf_shift.dat", "w");

   for (ii=0;ii<ntime;ii++) {
/*   
      for (jj=0;jj<N_sprk_u;jj++)
	 fprintf (fpt, " %d 0 %lf %d\n", ii, th_sprk_u[jj], N_sprk_u);

      for (jj=0;jj<N_sprk_d;jj++)
         fprintf (fpt, " %d 1 %lf %d\n", ii, th_sprk_d[jj], N_sprk_d);
        
      fprintf (fpt, "\n");
*/

      edgedist (x_plot, y_plot, th_sprk_u, th_sprk_d, edg_sft_u, edg_sft_d, 
		   a_cap, b_cap, th_cap, th_lag, h_sprk, h_drft, x_cent, 
		   y_cent, N_sprk_u, N_sprk_d, &nplt);


      //Plotting the angular shift of minimum distance
      for (jj=0;jj<nplt;jj++) {
         x_plot[jj]*=180/M_PI;
	 y_plot[jj]*=180/M_PI;

         fprintf (fpt, "%lf  %lf\n", x_plot[jj], y_plot[jj]);
      }

      xmin = 0.0;
      xmax = 360.0;

      ymin = y_plot[0];
      ymax = y_plot[0];

      for (jj=0;jj<nplt;jj++) {

         if (y_plot[jj] < ymin) 
            ymin = y_plot[jj];

         if (y_plot[jj] > ymax)
            ymax = y_plot[jj];
      }
      y_diff = 0.2*(ymax-ymin);

      ymin-=y_diff;
      ymax+=y_diff;

      cpgsvp(xleft, xright, ybot, ytop);

      cpgbox("BC", 0.0, 0, "BC", 0.0, 0);

      cpgswin(xmin, xmax, ymin, ymax);
   
      cpgaxis("N", xmin, ymin, xmin, ymax, ymin, ymax, 0.0, 0, 0.5, 0.5, 0.5, 
                   -1.0, 0.0);

      cpgaxis("N", xmin, ymin, xmax, ymin, xmin, xmax, 60.0, 6, 0.5, 0.5, 0.5, 
                   1.0, 0.0);

      cpgsch(1.4);

      cpglab ("angular dist (deg)", "angular shift (deg)", "");

      cpgsch(1.2);

      cpgslw(12);

      cpgpt (nplt, x_plot, y_plot, 17);

      cpgslw(1);

      delay (del_sec);

      cpgask (1);
     

      //Updating spark location for next time

      //Upper half

      th_buf[0] = th_sprk_u[0] + edg_sft_u[0];
      sp_cnt = 1;

      if (th_buf[0] < M_PI - theta_sp) {
	 th_buf[0]+=theta_sp;
         th_buf[1] = th_sprk_u[0] + edg_sft_u[0];
	 sp_cnt = 2;
      }

      for (jj=1;jj<N_sprk_u-1;jj++)  
         th_buf[sp_cnt++] = th_sprk_u[jj]+0.5*(edg_sft_u[jj-1]+edg_sft_u[jj]);
      

      if (th_sprk_u[N_sprk_u-1] + edg_sft_u[N_sprk_u-2] > 0)
	 th_buf[sp_cnt++] = th_sprk_u[N_sprk_u-1] + edg_sft_u[N_sprk_u-2];

      th_sprk_u[0] = th_buf[0];

      for (jj=1;jj<sp_cnt-1;jj++) {

         th_sprk_u[jj] = th_buf[jj];

         if (th_sprk_u[jj] < th_buf[jj+1]+theta_sp)
            th_sprk_u[jj] = th_buf[jj+1]+theta_sp;
      }
      th_sprk_u[sp_cnt-1] = th_sprk_u[sp_cnt-2]-theta_sp;

      N_sprk_u = sp_cnt;


      //Lower half
      sp_cnt = 0;

      th_buf[sp_cnt++] = th_sprk_d[0] + edg_sft_d[0];

      if (th_buf[sp_cnt-1] > M_PI + theta_sp) {
         th_buf[sp_cnt] = th_sprk_u[0] + edg_sft_u[0];
         th_buf[sp_cnt-1]-=theta_sp;
         sp_cnt++;
      }

      for (jj=1;jj<N_sprk_d-1;jj++) 
         th_buf[sp_cnt++] = th_sprk_d[jj]+0.5*(edg_sft_d[jj-1]+edg_sft_d[jj]);

      if (th_sprk_d[N_sprk_d-1] + edg_sft_d[N_sprk_d-2] < 2*M_PI)
         th_buf[sp_cnt++] = th_sprk_d[N_sprk_d-1] + edg_sft_d[N_sprk_d-2];
      
      th_sprk_d[0] = th_buf[0];

      for (jj=1;jj<sp_cnt-1;jj++) {

         th_sprk_d[jj] = th_buf[jj];

         if (th_sprk_d[jj] < th_buf[jj+1]-theta_sp)
            th_sprk_d[jj] = th_buf[jj+1]-theta_sp;
      }
      th_sprk_d[sp_cnt-1] = th_sprk_d[sp_cnt-2]+theta_sp;

      N_sprk_d = sp_cnt;
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


void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out)
{
   *x_out =  x_in*cos(theta) + y_in*sin(theta);

   *y_out = -x_in*sin(theta) + y_in*cos(theta);

}


void edgedist (float *x_plot, float *y_plot, double *th_sprk_u, 
		  double *th_sprk_d, double *edg_sft_u, double *edg_sft_d, 
		  double a_cap, double b_cap, double th_cap, double th_lag, 
		  double h_sprk, double h_drft, double x_cent, double y_cent, 
		  int N_sprk_u, int N_sprk_d, int *nplt)
{
   double       x_sprk[2], y_sprk[2], x_sft[2], y_sft[2];
   double       rad_dist[2], rad_sft[2], d_min[2], th_min[2];
   double       a_out, a_in, b_out, b_in;
   double       th_val, th_h;
   double       x_in, y_in, x_val, y_val, x_cel, y_cel;
   int          ii, jj, ncnt;


   a_out = a_cap;
   a_in = a_out - h_sprk;

   b_out = b_cap;
   b_in = b_out - b_cap/a_cap*h_sprk;

   ncnt = 0;

   //Upper half
   for (ii=0;ii<N_sprk_u-1;ii++) {

      //Cartesian position of sparks
      for (jj=0;jj<=1;jj++) {

	 x_in = 0.5*(a_out+a_in)*cos(th_sprk_u[ii+jj] - th_cap);
         y_in = 0.5*(b_out+b_in)*sin(th_sprk_u[ii+jj] - th_cap);

	 coortrans (x_in, y_in, -th_cap, &x_sprk[jj], &y_sprk[jj]);

	 x_sprk[jj]+=x_cent;
	 y_sprk[jj]+=y_cent;

         x_sft[jj] = x_sprk[jj] + h_drft*cos(th_lag);
         y_sft[jj] = y_sprk[jj] + h_drft*sin(th_lag);
      }

      //Finding the shift of minimum point
      th_val = th_sprk_u[ii];
      th_h = th_sprk_u[ii+1]; 

      d_min[0] = 2*a_cap;
      d_min[1] = 2*a_cap;

      while (th_val >= th_h) {

         x_in = a_out*cos(th_val-th_cap);
         y_in = b_out*sin(th_val-th_cap);

	 coortrans (x_in, y_in, -th_cap, &x_val, &y_val);

	 x_val+=x_cent;
	 y_val+=y_cent;

         for (jj=0;jj<=1;jj++) {

	    //Distance for edge spark
            x_cel = x_val - x_sprk[jj];
            y_cel = y_val - y_sprk[jj];

	    rad_dist[jj] = sqrt(x_cel*x_cel + y_cel*y_cel);

	    //Distance for shifted spark
	    x_cel = x_val - x_sft[jj];
            y_cel = y_val - y_sft[jj];

	    rad_sft[jj] = sqrt(x_cel*x_cel + y_cel*y_cel);
         }

         if (rad_dist[0]+rad_dist[1] < d_min[0]) {
            d_min[0] = rad_dist[0]+rad_dist[1];
	    th_min[0] = th_val;
         }

         if (rad_sft[0]+rad_sft[1] < d_min[1]) {
            d_min[1] = rad_sft[0]+rad_sft[1];
            th_min[1] = th_val;
         }

         th_val-=0.1*h_drft/(a_out+a_in);
      }
      edg_sft_u[ii] = th_min[1]-th_min[0];

      x_plot[ncnt] = (float)(0.5*(th_sprk_u[ii]+th_sprk_u[ii+1]));
      y_plot[ncnt] = (float)(edg_sft_u[ii]);
      ncnt++;
   }


   //lower half
   for (ii=0;ii<N_sprk_d-1;ii++) {

      //Cartesian position of sparks
      for (jj=0;jj<=1;jj++) {

	 x_in = 0.5*(a_out+a_in)*cos(th_sprk_d[ii+jj] - th_cap);
         y_in = 0.5*(b_out+b_in)*sin(th_sprk_d[ii+jj] - th_cap);

	 coortrans (x_in, y_in, -th_cap, &x_sprk[jj], &y_sprk[jj]);

	 x_sprk[jj]+=x_cent;
	 y_sprk[jj]+=y_cent;

         x_sft[jj] = x_sprk[jj] + h_drft*cos(th_lag);
         y_sft[jj] = y_sprk[jj] + h_drft*sin(th_lag);
      }

      //Estimating shift in minimum point
      th_val = th_sprk_d[ii];
      th_h = th_sprk_d[ii+1];

      d_min[0] = 2*a_out;
      d_min[1] = 2*a_out;

      while (th_val <= th_h) {

         x_in = a_out*cos(th_val-th_cap);
         y_in = b_out*sin(th_val-th_cap);

         coortrans (x_in, y_in, -th_cap, &x_val, &y_val);

         x_val+=x_cent;
         y_val+=y_cent;

         for (jj=0;jj<=1;jj++) {

            //Distance for edge spark
            x_cel = x_val - x_sprk[jj];
            y_cel = y_val - y_sprk[jj];

            rad_dist[jj] = sqrt(x_cel*x_cel + y_cel*y_cel);

            //Distance for shifted spark
            x_cel = x_val - x_sft[jj];
            y_cel = y_val - y_sft[jj];

            rad_sft[jj] = sqrt(x_cel*x_cel + y_cel*y_cel);
         }

         if (rad_dist[0]+rad_dist[1] < d_min[0]) {
            d_min[0] = rad_dist[0]+rad_dist[1];
            th_min[0] = th_val;
         }

         if (rad_sft[0]+rad_sft[1] < d_min[1]) {
            d_min[1] = rad_sft[0]+rad_sft[1];
            th_min[1] = th_val;
         }

         th_val+=0.1*h_drft/(a_out+a_in);
      }
      edg_sft_d[ii] = th_min[1]-th_min[0];

      x_plot[ncnt] = (float)(0.5*(th_sprk_d[ii]+th_sprk_d[ii+1]));
      y_plot[ncnt] = (float)(edg_sft_d[ii]);
      ncnt++;
   }
   *nplt = ncnt;
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
