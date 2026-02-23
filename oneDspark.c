// Estimating the spar motion in one dimension 
//
// To compile: gcc oneDspark.c -o ~/bin/oneDspark -lm

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
   usage: oneDspark ntime\n"


void cmdlineError (char *message);


int main(int argc, char *argv[])
{
   FILE        *fpt;
   double       D_cap, h_sprk, h_drft;
   double       x_val, x_strt;
   double       sprk_pos, sprk_wid, sprk_cnt, I_val;
   int          ii, jj, ntime;


   if (argc < 2)
      cmdlineError (USAGE);

   ntime = atoi(argv[1]);


   //Initializing parameters

   D_cap = 15.0;  // polar cap diameter in meters.
   h_sprk = 2.6;  // spark diameter in meters.
   h_drft = 0.1;  // drift movement in meters.

   fpt = fopen ("sparktime.dat", "w");

   x_strt = 0;

   for (ii=0;ii<ntime;ii++) {
      
      x_val = x_strt+h_sprk/2;

      //First small spark
      sprk_pos = 0;
      sprk_cnt = (x_val-h_sprk/2)/2;
      sprk_wid = sprk_cnt/3;

      while (sprk_pos <= x_val-h_sprk/2-2*h_drft) {

         I_val = exp(-0.5*(pow(sprk_pos-sprk_cnt,2.0))/(sprk_wid*sprk_wid));

         fprintf (fpt, " %lf  %lf %d\n", sprk_pos, I_val, ii);

         sprk_pos+=h_drft;
      }

      while (x_val <= D_cap-h_sprk/2) {
        
	 //Writing Gaussian Sparks
         sprk_pos = x_val-h_sprk/2;
	 sprk_cnt = x_val;
	 sprk_wid = h_sprk/3;

	 while (sprk_pos < x_val+h_sprk/2) {

            I_val = exp(-0.5*(pow(sprk_pos-sprk_cnt,2.0))/(sprk_wid*sprk_wid));
	    
            fprintf (fpt, " %lf  %lf %d\n", sprk_pos, I_val, ii);

	    sprk_pos+=h_drft;
	 }

	 x_val+=h_sprk;
      }
      fprintf (fpt, "\n");

      x_strt+=h_drft;

      if (x_strt > h_sprk)
         x_strt-=h_sprk;
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

