# define LINELEN (1000) // The word size of line in files.


double pulsfft (double *pulstk, double *fftout, int npulse, int nbin, int bwin,                  int ewin)
{
   int                   ii, jj, ncnt, indx;

   // FFT variables
   double *in;
   fftw_complex  *out;
   fftw_plan p1;


   in = (double*) fftw_malloc(sizeof(double)*(npulse+2));
   out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex)*(npulse/2+1));

   for (jj=0;jj<nbin;jj++) {
      ncnt = 0;

      for (ii=0;ii<npulse;ii++)
         in[ncnt++] = pulstk[ii*nbin+jj];

      p1 = fftw_plan_dft_r2c_1d(ncnt, in, out, FFTW_ESTIMATE);

      fftw_execute(p1);

      for (ii=0;ii<=npulse/2;ii++) {
         indx = jj*npulse+ii*2;
         fftout[indx] = out[ii][0]*2.0/(double)(npulse);
         fftout[indx+1] = out[ii][1]*2.0/(double)(npulse);
      }
   }

   fftw_destroy_plan(p1);
   fftw_free(out);
   fftw_free(in);
}


//Plotting function
void plotlrfs (char *infile, double *pulstk, double *fftout, int fftlen,
                int nbin, int bwin, int ewin)
{
   FILE              *fpt;
   char               lrfsfile[1000], foldfile[1000];
   char               psfile[1000], label[1000];
   float             *plotstk, *plotavg, *yavg;
   float             *xfold, *yfold;
   double             rad, theta, fftx, ffty;
   double            *radmax, *thetamax, statmax[2];
   float             *plotphs, *xphs;
   int                maxpos, spos, epos, ncnt;
   int                ii, jj, indx;

   // pgplot inputs
   float              tr[] = { 0, 1, 0, 0, 0, 1 };  // Identity Mapping
   float              HL[] = { 0.0, 0.2, 0.4, 0.6, 1.0 };
   float              HR[] = { 0.0, 0.5, 1.0, 1.0, 1.0 };
   float              HG[] = { 0.0, 0.0, 0.5, 1.0, 1.0 };
   float              HB[] = { 0.0, 0.0, 0.0, 0.3, 1.0 };
   float              bri, cont, nc;
   float              xleft, xright, ybot, ytop;
   float              xmin, xmax, ymin, ymax, hi, lo;


   plotstk  = (float *)calloc(nbin*fftlen/2, sizeof(float));
   plotavg  = (float *)calloc(fftlen/2, sizeof(float));
   yavg     = (float *)calloc(fftlen/2, sizeof(float));
   xfold    = (float *)calloc(nbin, sizeof(float));
   yfold    = (float *)calloc(nbin, sizeof(float));
   radmax   = (double *)calloc(nbin, sizeof(double));
   thetamax = (double *)calloc(nbin, sizeof(double));
   plotphs  = (float *)calloc(nbin, sizeof(float));
   xphs  = (float *)calloc(nbin, sizeof(float));


   sprintf(psfile,"%s_lrfsout.ps/VCPS", infile);

   // Plotting FFT
   fftx = fftout[2];
   ffty = fftout[3];
   rec2pol (fftx, ffty, &rad, &theta);
   hi = (float)rad;
   lo = (float)rad;

   for (ii=1;ii<=fftlen/2;ii++) {
      for (jj=bwin;jj<ewin;jj++) {
         indx = (ii-1)*(ewin-bwin)+jj-bwin;
         fftx = fftout[jj*fftlen+ii*2];
         ffty = fftout[jj*fftlen+ii*2+1];

         rec2pol (fftx, ffty, &rad, &theta);
         plotstk[indx] = (float)rad;

         if (plotstk[indx] > hi)
            hi = plotstk[indx];

         if (plotstk[indx] < lo)
            lo = plotstk[indx];
      }
   }
   cpgbeg(0, psfile, 1, 1);

   cpgslw(2);
   cpgsch(1.2);
   cpgscf(2);

   // Plotting the FFT stack  
   xleft =  0.4;
   xright = 0.92;
   ybot = 0.25;
   ytop = 0.78;
   cpgsvp(xleft, xright, ybot, ytop);

   xmin = 0;
   xmax = ewin-bwin-1;

   fprintf(stderr,"xmin %f xmax %f  %f  %f\n", xmin, xmax, hi, lo);
   ymin = 0;
   ymax = fftlen/2;
   cpgswin(xmin, xmax, ymin, ymax);

   cont = -1.0;
   bri = 0.5;
   nc = 5;
   cpgctab (HL, HR, HG, HB, nc, cont, bri);

   cpgimag (plotstk, ewin-bwin, fftlen/2-1, 1, ewin-bwin, 1, fftlen/2-1, hi,
             lo, tr);

   cpgbox ("BC", 0.0, 0, "BC", 0.0, 0);


   //estimating full profile
   for (ii=0;ii<nbin;ii++)
      yfold[ii] = 0.0;

   for (ii=0;ii<nbin;ii++) {
      for (jj=0;jj<fftlen;jj++)
         yfold[ii]+=(float)(pulstk[jj*nbin+ii]/fftlen);
   }

   sprintf (foldfile, "%s_foldout.dat", infile);

   fpt = fopen(foldfile, "w");

   for (jj=0;jj<nbin;jj++)
      fprintf (fpt,"%d %f\n", jj, yfold[jj]);

   fclose(fpt);

   // Plotting folded profile
   for (ii=0;ii<(ewin-bwin);ii++)
      yfold[ii] = 0.0;

   for (ii=0;ii<(ewin-bwin);ii++) {
      xfold[ii] = (float)(ii+bwin)/(float)(nbin)*360-180;
      for (jj=0;jj<fftlen;jj++)
         yfold[ii]+=(float)(pulstk[jj*nbin+ii+bwin]/fftlen);
   }

   ymin = yfold[0];
   ymax = yfold[0];

   for (ii=0;ii<(ewin-bwin);ii++) {

      if (yfold[ii] > ymax)
         ymax = yfold[ii];

      if (yfold[ii] < ymin)
         ymin = yfold[ii];
   }

   xmin = (float)(bwin)/(float)(nbin)*360-180;
   xmax = (float)(ewin)/(float)(nbin)*360-180;

   xleft  = 0.4;
   xright = 0.92;
   ybot = 0.075;
   ytop = 0.25;

   cpgsvp (xleft, xright, ybot, ytop);

   cpgswin (xmin, xmax, ymin-0.1*(ymax-ymin), ymax+0.2*(ymax-ymin));

   cpglab("Pulsar Longitude [degrees]"," ", " ");

   cpgbox("BCNTS", 0.0, 0, "BCNTSV", 0.0, 0);

   cpgsci(2);
   cpgline(ewin-bwin, xfold, yfold);


   //Plotting average fourier spectra
   sprintf (lrfsfile, "%s_lrfs.dat", infile);

   fpt = fopen(lrfsfile, "w");

   for (ii=0;ii<fftlen/2-1;ii++)
      plotavg[ii] = 0.0;

   for (ii=1;ii<fftlen/2;ii++) {
      yavg[ii-1] = (float)ii/(float)(fftlen);

      for (jj=0;jj<nbin;jj++) {
         fftx = fftout[jj*fftlen+ii*2];
         ffty = fftout[jj*fftlen+ii*2+1];

         rec2pol (fftx, ffty, &rad, &theta);
         plotavg[ii-1]+=(float)rad/(float)(ewin-bwin);
      }
      fprintf(fpt,"%d %lf\n", ii, plotavg[ii-1]);
   }
   fclose(fpt);

   xmin = plotavg[0];
   xmax = plotavg[0];
   maxpos = 1;

   for (ii=0;ii<fftlen/2-1;ii++) {

      if (plotavg[ii] > xmax) {
         xmax = plotavg[ii];
         maxpos = ii+1;
      }

      if (plotavg[ii] < xmin)
         xmin = plotavg[ii];
   }

   //Normalising by max Peak
   for (ii=0;ii<fftlen/2;ii++)
      plotavg[ii]/=xmax;

   xmin/=xmax;
   xmax = 1;

   ymin = 0.0;
   ymax = 0.5;

   fprintf(stderr,"xmin %f xmax %f  maxpos %d\n", xmin, xmax, maxpos);

   xleft  = 0.1;
   xright = 0.4;
   ybot = 0.25;
   ytop = 0.78;

   cpgsvp (xleft, xright, ybot, ytop);
   cpgsci(1);
   cpgswin ( xmax+0.1*(xmax-xmin), xmin-0.2*(xmax-xmin), ymin, ymax);

   cpgslw(2);
   cpglab(" ","frequency (1/P)", " ");

   cpgbox ("BCMTS", 0.0, 0, "BCNTS", 0.0, 0);

   cpgmtxt("T", 2, -0.07, 0.0, "");

   cpgslw(2);
   cpgsci(4);
   cpgline(fftlen/2-1, plotavg, yavg);


   //Plotting the phase at maximum FFT
   for (jj=0;jj<nbin;jj++) {
      fftx =fftout[jj*fftlen+maxpos*2];
      ffty =fftout[jj*fftlen+maxpos*2+1];

      rec2pol (fftx, ffty, &rad, &theta);

      radmax[jj] = rad;
      thetamax[jj] = atan2(ffty,fftx)*180.0/M_PI;
   }

   findbase(radmax, nbin, 3, &spos, &epos, statmax);

   ncnt = 0;
   for (jj=bwin;jj<ewin;jj++) {
      if (radmax[jj] > statmax[0]+3.0*statmax[1]) {
         xphs[ncnt] = (float)(jj)/(float)(nbin)*360-180;
         plotphs[ncnt] = (float)thetamax[jj];

          
//         plotphs[ncnt] = (float)thetamax[jj]-60.0;

         if (plotphs[ncnt] < -180)
            plotphs[ncnt]+=360;
         ncnt++;
      }
   }

   xmin = (float)(bwin)/(float)(nbin)*360-180;
   xmax = (float)(ewin)/(float)(nbin)*360-180;

   ymin = -179.999;
//   ymin = 0.0;
   ymax = 179.999;

   xleft  = 0.4;
   xright = 0.92;
   ybot = 0.78;
   ytop = 0.88;

   cpgsvp (xleft, xright, ybot, ytop);
   cpgsci(1);
   cpgswin (xmin, xmax, ymin, ymax);

   cpgslw(2);
   cpgsch(0.9);

   cpgmtxt("T", 2, 0.4, 0.0, "Phase [deg]");

   cpgbox ("BCMTS", 0.0, 0, "BCMTV", 90.0, 0);

   cpgslw(8);
   cpgsci(1);
   cpgpt(ncnt, xphs, plotphs, -1);

   cpgclos();
}


// Plotting the single pulses
void plotsngl (char *filename, double *pulstk, double alpha, double beta,
                int npulse, int nbin, int bwin, int ewin)
{
   char               psfile[LINELEN];
   float             *plotstk, *xfold, *yfold, *plotavg, *yavg, *arrppa;
   double             xppa, yppa, phi;
   int                sbin, ebin;
   int                ii, jj, kk, indx;

   // pgplot inputs
   float              tr[] = { 0, 1, 0, 0, 0, 1 };  // Identity Mapping
   float              HL[] = { 0.0, 0.2, 0.4, 0.6, 1.0 };
   float              HR[] = { 0.0, 0.5, 1.0, 1.0, 1.0 };
   float              HG[] = { 0.0, 0.0, 0.5, 1.0, 1.0 };
   float              HB[] = { 0.0, 0.0, 0.0, 0.3, 1.0 };

   float   WL[] = {0.0, 0.5, 0.5, 0.7, 0.7, 0.85, 0.85, 0.95, 0.95, 1.0 };
   float   WR[] = {0.0, 1.0, 0.0, 0.0, 0.3,  0.8,  0.3,  1.0,  1.0, 1.0};
   float   WG[] = {0.0, 0.5, 0.4, 1.0, 0.0,  0.0,  0.2,  0.7,  1.0, 1.0};
   float   WB[] = {0.0, 0.0, 0.0, 0.0, 0.4,  1.0,  0.0,  0.0, 0.95, 1.0};

   float   RL[] = {-0.5, 0.0, 0.17, 0.33, 0.50, 0.67, 0.83, 1.0, 1.7};
   float   RR[] = {0.0, 0.0,  0.0,  0.0,  0.6,  1.0,  1.0, 1.0, 1.0};
   float   RG[] = {0.0, 0.0,  0.0,  1.0,  1.0,  1.0,  0.6, 0.0, 1.0};
   float   RB[] = {0.0, 0.3,  0.8,  1.0,  0.3,  0.0,  0.0, 0.0, 1.0};

   float   AL[] = {0.0, 0.1, 0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5,
                   0.5, 0.6, 0.6, 0.7, 0.7, 0.8, 0.8, 0.9, 0.9, 1.0};
   float   AR[] = {0.0, 0.0, 0.3, 0.3, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0,
                   0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
   float   AG[] = {0.0, 0.0, 0.3, 0.3, 0.0, 0.0, 0.0, 0.0, 0.8, 0.8,
                   0.6, 0.6, 1.0, 1.0, 1.0, 1.0, 0.8, 0.8, 0.0, 0.0};
   float   AB[] = {0.0, 0.0, 0.3, 0.3, 0.7, 0.7, 0.7, 0.7, 0.9, 0.9,
                   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

   float              bri, cont, nc;
   float              xleft, xright, ybot, ytop;
   float              xmin, xmax, ymin, ymax, hi, lo;


   sbin = bwin-(ewin-bwin)/4;
   if (sbin < 0)
      sbin = 0;

   ebin = ewin+(ewin-bwin)/4;
   if (ebin > nbin-1)
      ebin = nbin-1;

   fprintf (stderr," sbin %d  ebin %d \n", sbin, ebin);

   sprintf(psfile,"%s_pulstack.ps/VCPS", filename);

   // Plotting pulse stack within window 
   plotstk  = (float *)calloc(nbin*npulse, sizeof(float));

   hi = (float)pulstk[sbin];
   lo = (float)pulstk[sbin];

   for (ii=0;ii<npulse;ii++) {
      for (jj=sbin;jj<ebin;jj++) {
         indx = ii*(ebin-sbin)+jj-sbin;
         plotstk[indx] = (float)pulstk[ii*nbin+jj];

         if (plotstk[indx] > hi)
            hi = plotstk[indx];

         if (plotstk[indx] < lo)
            lo = plotstk[indx];
      }
   }

   cpgbeg(0, psfile, 1, 1);

   cpgslw(2);
   cpgsch(1.2);
   cpgscf(2);

   // Plotting the pulse stack  
   xleft =  0.4;
   xright = 0.92;
   ybot = 0.25;
   ytop = 0.78;
   cpgsvp(xleft, xright, ybot, ytop);

   xmin = 0;
   xmax = ebin-sbin-1;

   fprintf(stderr,"xmin %f xmax %f  %f  %f\n", xmin, xmax, hi, lo);
   ymin = 0;
   ymax = npulse-1;
   cpgswin(xmin, xmax, ymin, ymax);

   cont = -1.0;
   bri = 0.5;
   nc = 5;
   cpgctab (HL, HR, HG, HB, nc, cont, bri);

   cpgbox ("BC", 0.0, 0, "BC", 0.0, 0);

   cpgimag (plotstk, ebin-sbin, npulse-1, 1, ebin-sbin, 1, npulse-1, hi,
             lo, tr);


   // Plotting folded profile
   xfold = (float *)calloc(nbin, sizeof(float));
   yfold = (float *)calloc(nbin, sizeof(float));

   for (ii=0;ii<(ebin-sbin);ii++)
      yfold[ii] = 0.0;

   for (ii=0;ii<(ebin-sbin);ii++) {
      xfold[ii] = (float)(ii+sbin)/(float)(nbin)*360-180;
      for (jj=0;jj<npulse;jj++)
         yfold[ii]+=(float)(pulstk[jj*nbin+ii+sbin]/npulse);
   }

   ymin = yfold[0];
   ymax = yfold[0];

   for (ii=0;ii<(ebin-sbin);ii++) {

      if (yfold[ii] > ymax)
         ymax = yfold[ii];

      if (yfold[ii] < ymin)
         ymin = yfold[ii];
   }

   xmin = (float)(sbin)/(float)(nbin)*360-180;
   xmax = (float)(ebin)/(float)(nbin)*360-180;

   xleft  = 0.4;
   xright = 0.92;
   ybot = 0.075;
   ytop = 0.25;

   cpgsvp (xleft, xright, ybot, ytop);

   cpgswin (xmin, xmax, ymin-0.1*(ymax-ymin), ymax+0.2*(ymax-ymin));

   cpglab("Pulsar Longitude [degrees]"," ", " ");

   cpgbox("BCNTS", 0.0, 0, "BCNTSV", 0.0, 0);

   cpgsci(2);
   cpgline(ebin-sbin, xfold, yfold);


   //Plotting average energy
   plotavg = (float *)calloc(npulse, sizeof(float));
   yavg = (float *)calloc(npulse, sizeof(float));

   for (ii=0;ii<npulse;ii++)
      plotavg[ii] = 0.0;

   for (ii=0;ii<npulse;ii++) {
      yavg[ii] = (float)ii;

      for (jj=bwin;jj<ewin;jj++)
         plotavg[ii]+=(float)pulstk[ii*nbin+jj]/(float)(ewin-bwin);
   }

   xmin = plotavg[0];
   xmax = plotavg[0];

   for (ii=0;ii<npulse;ii++) {

      if (plotavg[ii] > xmax)
         xmax = plotavg[ii];

      if (plotavg[ii] < xmin)
         xmin = plotavg[ii];
   }

   xleft  = 0.1;
   xright = 0.4;
   ybot = 0.25;
   ytop = 0.78;

   ymin = 0;
   ymax = (float)(npulse-1);

   cpgsvp (xleft, xright, ybot, ytop);
   cpgsci(1);
   cpgswin (xmax+0.2*(xmax-xmin), xmin-0.2*(xmax-xmin), ymin, ymax);

   cpgslw(2);
   cpgbox ("BCMT", 0.0, 0, "BCNTS", 0.0, 0);
   cpglab("","Pulse Number", "energy[arb. units]");

   cpgmtxt("T", 2, -0.07, 0.0, "");

   cpgslw(2);
   cpgsci(4);
   cpgline(npulse-1, plotavg, yavg);
   cpgsci(1);


   // Writing Polarization Position Angle (PPA)
   arrppa = (float *)calloc(nbin, sizeof(float));

   for (ii=0;ii<(ebin-sbin);ii++) {
      phi = (double)xfold[ii]/180*M_PI;

      xppa = sin(alpha)*sin(phi);
      yppa = sin(alpha+beta)*cos(alpha)-sin(alpha)*cos(alpha+beta)*cos(phi);

      arrppa[ii] = atan2(xppa,yppa)*180/M_PI;
   }

   xmin = (float)(sbin)/(float)(nbin)*360-180;
   xmax = (float)(ebin)/(float)(nbin)*360-180;

   xleft  = 0.4;
   xright = 0.92;
   ybot = 0.78;
   ytop = 0.88;

   ymin = -179.999;
   ymax = 179.999;

   cpgsvp (xleft, xright, ybot, ytop);

   cpgswin (xmin, xmax, ymin, ymax);

//   cpglab("Pulsar Longitude (in degrees)"," ", " ");

   cpgslw(2);
   cpgsch(0.9);

   cpgmtxt("T", 2, 0.4, 0.0, "PPA [deg]");

   cpgbox("BCMTS", 0.0, 0, "BCMTV", 90.0, 0);

   cpgslw(4);
   cpgsci(1);
//   cpgline(ebin-sbin, xfold, arrppa);
   cpgpt (ebin-sbin, xfold, arrppa, -1);

   cpgsci(1);


   cpgclos();
}

