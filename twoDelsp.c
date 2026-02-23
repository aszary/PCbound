// Estimating the spark distribution in elliptical polar cap
//
// To compile: gcc twoDelsp.c -o ~/bin/twoDelsp -lm
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


# define C_LIGHT (299792458) // The speed of light.
# define Hz2MHz (1000000.0)  // frequency conversion.
# define ZERO_lf (0.000001) // round off zero error.
# define LINELEN (10000)     // Size of line.
# define WORDSIZE (100)      // Size of filename.


#define USAGE "\
   usage: twoDelsp  incl_angle (deg)  maj_axis (m) ellip (a/b)\n"

void coortrans (double x_in, double y_in, double theta, double *x_out, 
		   double *y_out);

void cmdlineError (char *message);


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       a_cap, b_cap, h_sprk, h_drft;
   double       x_val, y_val, x_cent, y_cent;
   double       a_val, b_val, el_val, th_cap;
   double       a_sprk, b_sprk, a_in, a_out, b_in, b_out;
   double      *x_sprk, *y_sprk, *th_sprk;
   double       x_cel, y_cel, x_trfm, y_trfm;
   double       theta_sp, theta_c, theta_t;
   int         *N_sprk, tot_sprk, N_trk;
   int          ii, jj, ncnt;


   if (argc < 3)
      cmdlineError (USAGE);

   th_cap = atof(argv[1])/180*M_PI;
   a_cap = atof(argv[2]);
   b_cap = a_cap*atof(argv[3]);
   

   //Initializing parameters
   h_sprk = 2.6;     // spark diameter in meters.
   h_drft = 0.1;     // drift movement in meters.

   a_sprk = h_sprk/2;
   b_sprk = a_sprk*b_cap/a_cap;

   
   //Finding centre of rotated ellipse
   x_cent = sqrt(pow(a_cap*cos(th_cap),2) + pow(b_cap*sin(th_cap),2));
   y_cent = sqrt(pow(a_cap*sin(th_cap),2) + pow(b_cap*cos(th_cap),2));


   // The packing of sparks
   N_trk = b_cap/b_sprk/2;

   N_sprk = (int *)calloc(N_trk+1, sizeof(int));

   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   b_out = b_cap;
   b_in = b_out - 2*b_sprk;

   tot_sprk = 0;

   for (ii=0;ii<N_trk;ii++) {
      
      N_sprk[ii] = 0.75*(a_out*b_out - a_in*b_in)/(a_sprk*b_sprk);

      tot_sprk+=N_sprk[ii];

      a_out-=2*a_sprk;
      a_in-=2*a_sprk;
      
      b_out-=2*b_sprk;
      b_in-=2*b_sprk;
   }

   fprintf (stderr, "N_trk %d  tot_sprk %d\n", N_trk, tot_sprk);

   x_sprk = (double *)calloc(tot_sprk+1, sizeof(double));
   y_sprk = (double *)calloc(tot_sprk+1, sizeof(double));
   
   //Finding the central locations of sparks
   a_out = a_cap;
   a_in = a_out - 2*a_sprk;

   b_out = b_cap;  
   b_in = b_out - 2*b_sprk;

   ncnt = 0;

   for (ii=0;ii<N_trk;ii++) {
     
      theta_sp = 2*M_PI/N_sprk[ii];

      theta_c = theta_sp/2 - th_cap;

      for (jj=0;jj<N_sprk[ii];jj++) {

	 x_val = 0.5*(a_out+a_in)*cos(theta_c);

	 y_val = 0.5*(b_out+b_in)*sin(theta_c);

	 coortrans (x_val, y_val, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

	 th_sprk[ncnt] = theta_c;

	 x_sprk[ncnt]+=x_cent;
	 y_sprk[ncnt]+=y_cent;

         theta_c+=theta_sp;
	 ncnt++;
      }

      a_out-=2*a_sprk;
      a_in-=2*a_sprk;

      b_out-=2*b_sprk;
      b_in-=2*b_sprk;
   }


   //location of the polar cap within each spark
   fpt = fopen ("spark2Delisp.dat", "w");

   x_val = h_drft/2;

   a_in = a_cap - N_trk*2*a_sprk - 1.5*h_drft;
   b_in = b_cap - N_trk*2*b_sprk - 1.5*h_drft;
   
   while (x_val < 2*x_cent) {

      y_val = h_drft/2;

      while (y_val < 2*y_cent) {

         for (ii=0;ii<tot_sprk;ii++) {

            x_cel = x_val - x_sprk[ii];
	    y_cel = y_val - y_sprk[ii];

            coortrans(x_cel, y_cel, th_cap, &x_trfm, &y_trfm); 

            el_val = sqrt(x_trfm*x_trfm/a_sprk/a_sprk + 
			    y_trfm*y_trfm/b_sprk/b_sprk);

            if (el_val < 1) 
               fprintf (fpt, "%lf  %lf  %d\n", x_val, y_val, ii);
         }

	 //Core spark
	 x_cel = x_val - x_cent;
	 y_cel = y_val - y_cent;

	 coortrans(x_cel, y_cel, th_cap, &x_trfm, &y_trfm);

	 el_val = (x_trfm*x_trfm/a_in/a_in + y_trfm*y_trfm/b_in/b_in);

	 if (el_val < 1) 
            fprintf (fpt, "%lf  %lf  -1\n", x_val, y_val);


         y_val+=h_drft;
      }

      x_val+=h_drft;
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

