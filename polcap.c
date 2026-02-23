//  To generate the open field line region on the surface from a given 
//  magnetic field configuration.

//  To compile: gcc polcap.c -o ~/bin/polcap -lm
   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ellipse_fit.h"


# define WORDSIZE (100) // The size of file names.
# define LINELEN (1000) // The size of line in files.

#define USAGE "\
   usage: polcap infield outopen alpha\n"


void cmdlineError (char *message);

void nrerror(char error_text[]);

void coortrans (double *theta, double *phi, double alpha, double theta_in,
                 double phi_in);

double *vector(int nl, int nh);

int splitstr (char *iline, char *strarray[], char *sep, int nsep);

double dotprod (double *vect1, double *vect2);

void crossprod (double *vect1, double *vect2, double *crosvect);

void rotvect (double *invect, double *rotaxis, double *transvect,
                double th_rot);

void unitvect (double *invect, double *transvect);


int kmax=0,kount=0;
double *xp=0,**yp=0,dxsav=0;

/* The main program. */
int main(int argc, char *argv[])
{ 
   FILE        *fpt;
   char        *infile, *outfile;
   char        *readline, *readval[10];
   void         odeint(), equations(), free_vector();;
   double      *vector();
   double      *th, *dthdr;
   double      *outfield, *surf_m;
   double       RS, RLC, r1, r2, alpha, beta;
   double       delr, r, th_b, del_th, del_ph;
   double       th_in, ph_in;
   double       th_min, th_max, ph_min, ph_max, th_c, ph_c;
   double      *x_cap, *y_cap;
   double       inarr[3], rotaxis[3], trnsvect[3];
   double       x1, x2, y1, y2, z1, z2;
   double       x_cent, y_cent, z_cent, R_maj, R_min, R_avg;
   double       el_fact, el_angl, el_sig;
   int          nsurf, ncap;
   int          FixedCenter, FixedAngle, FixedEll;
   int          ii, jj, tt, pp, nok, nbad, nphi, indx, if_fit;


   infile = (char *)calloc(WORDSIZE, sizeof(char));
   outfile = (char *)calloc(WORDSIZE, sizeof(char));
   readline = (char *)calloc(LINELEN, sizeof(char));

   if (argc < 4)
      cmdlineError (USAGE);

   sprintf (infile, "%s", argv[1]);
   sprintf (outfile, "%s", argv[2]);
   alpha = atof(argv[3])/180*M_PI;

   //Reading the surface magnetic config file 
   if ((fpt = fopen(infile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", infile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);
   splitstr (readline, readval, " ", 1);
   nsurf = atoi(readval[1]);

   surf_m = (double *)calloc(6*nsurf, sizeof(double));

   for (ii=0;ii<nsurf;ii++) {
      fgets(readline, LINELEN, fpt);
      splitstr (readline, readval, " ", 1);

      for (jj=0;jj<6;jj++)
         surf_m[6*ii+jj] = atof(readval[jj]);
   }


   th = vector(1,2); 
   dthdr = vector(1,2); 

   RS = 10000;  // neutron star radius in meters.
   RLC = 299792458/(2*M_PI);  // Light cylinder radius.


   //Initial Conditions non aligned rotator
   r1 = 30*RS; 
   r2 = RS;
   delr = 0.01*RS;


   //Estimating the boundary of the outer field lines   
   th_b = asin(sqrt((r1+delr)/RLC));
   nphi = 100;

   outfield = (double *)calloc(12*(nphi+1), sizeof(double)); 
   x_cap = (double *)calloc(2*nphi+1, sizeof(double)); 
   y_cap = (double *)calloc(2*nphi+1, sizeof(double)); 

   th_in = th_b;
   del_ph = M_PI/nphi;
  
   pp = 0;
   ph_in = -M_PI;

   while (ph_in <= M_PI) {

      fprintf (stderr, "phi %lf \r", ph_in*180/M_PI);
      indx = 6*pp;
      r = r1;
      coortrans (&th[1], &th[2], alpha, th_in, ph_in);

      outfield[indx] = r;
      outfield[indx+1] = th[1];
      outfield[indx+2] = th[2];
 
 
      while (r > r2) {
         odeint(th, 2, r+delr, r, 1.0e-3, 0.05*delr, 1.0e-6, &nok, &nbad, 
                 equations, surf_m, nsurf, alpha);

         r-=delr;
      }
      outfield[indx+3] = r;
      if (th[1] < 0)
         th[1]+=2*M_PI;

      outfield[indx+4] = th[1];

      if (th[2] < -M_PI)
         th[2]+=2.0*M_PI;

      outfield[indx+5] = th[2];

      ph_in+=del_ph;
      pp++;
   }
   ncap = pp;


   //Rotating the Elliptical polar cap in the x-y plane
   th_min = outfield[4];
   th_max = outfield[4];
   ph_min = outfield[5];
   ph_max = outfield[5];

   for (ii=1;ii<ncap;ii++) {
      indx = 6*ii;

      if (th_min > outfield[indx+4])
         th_min = outfield[indx+4];
      if (th_max < outfield[indx+4])
         th_max = outfield[indx+4];

      if (ph_min > outfield[indx+5])
         ph_min = outfield[indx+5];
      if (ph_max < outfield[indx+5])
         ph_max = outfield[indx+5];
   }

   th_c = (th_min+th_max)/2;
   ph_c = (ph_min+ph_max)/2;

   fprintf (stderr, "%lf  %lf  %lf  %lf\n", th_min, th_max, th_c, ph_c);

   //Rotating the polar cap and fitting ellipse
   z_cent = RS*cos(th_c);
   x_cent = RS*sin(th_c)*cos(ph_c);
   y_cent = RS*sin(th_c)*sin(ph_c);

   inarr[0] = sin(th_c)*cos(ph_max) - sin(th_c)*cos(ph_min);
   inarr[1] = sin(th_c)*sin(ph_max) - sin(th_c)*sin(ph_min);
   inarr[2] = 0.0;

   unitvect (inarr, rotaxis);

   for (ii=0;ii<ncap;ii++) {
      indx = 6*ii;

      inarr[0] = outfield[indx+3]*sin(outfield[indx+4])*cos(outfield[indx+5]) -
	          x_cent;
      inarr[1] = outfield[indx+3]*sin(outfield[indx+4])*sin(outfield[indx+5]) -
	          y_cent;
      inarr[2] = outfield[indx+3]*cos(outfield[indx+4]) - z_cent;

      rotvect (inarr, rotaxis, trnsvect, -th_c);

      x_cap[ii] = trnsvect[0] + x_cent;
      y_cap[ii] = trnsvect[1] + y_cent;
   }

   fprintf (stderr, "%lf  %lf  %lf\n", x_cent, y_cent, z_cent);

   FixedCenter = 0;
   FixedAngle = 0;
   FixedEll = 0;

   if_fit = FitEllipse (x_cap, y_cap, ncap, &x_cent, &y_cent, FixedCenter,
                         FixedAngle, FixedEll, &R_maj, &R_min, &R_avg,
                         &el_fact, &el_angl, &el_sig);


   // Writing field lines in output file
   fpt = fopen(outfile,"w");
  
   fprintf(fpt, "#nfield %d  theta = %lf (rad)/ %lf (deg)\n", pp, th_b, th_b*180/M_PI);

   fprintf (fpt, "#x_cent %lf  y_cent %lf  R_maj %lf  R_min %lf  el_fact %lf  el_angl %lf\n", x_cent, y_cent, R_maj, R_min, el_fact, el_angl);


   for (jj=0;jj<ncap;jj++ ) { 
      indx = 6*jj;

      x1 = outfield[indx]*sin(outfield[indx+1])*cos(outfield[indx+2]);
      y1 = outfield[indx]*sin(outfield[indx+1])*sin(outfield[indx+2]);
      z1 = outfield[indx]*cos(outfield[indx+1]);

      x2 = outfield[indx+3]*sin(outfield[indx+4])*cos(outfield[indx+5]);
      y2 = outfield[indx+3]*sin(outfield[indx+4])*sin(outfield[indx+5]);
      z2 = outfield[indx+3]*cos(outfield[indx+4]);

      fprintf (fpt, " %lf %lf %lf  %lf %lf %lf    %lf %lf %lf  %lf %lf %lf\n", outfield[indx], outfield[indx+1], outfield[indx+2], x1, y1, z1, outfield[indx+3], outfield[indx+4], outfield[indx+5], x2, y2, z2);
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


void nrerror(char error_text[])
{ 
   void exit();
   fprintf(stderr,"Numerical Recipes run-time error...\n");
   fprintf(stderr,"%s\n",error_text);
   fprintf(stderr,"...now exiting to the system...\n");
}

/* A routine to allocate a subscripted array. */
double *vector(int nl, int nh)
{
   double *v;
   v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
   if (!v) nrerror("allocation failure in vector()");
   return v-nl;
}

/* This frees the array when no longer needed. */
void free_vector(double *v, int nl, int nh)
{ 
   free((char*) (v+nl)); 
}


int splitstr(char *iline, char *strarray[], char *sep, int nsep)
{
   int          ii, jj, kk, nword;
   int          llen;
   int          inword_f, if_sep;
   char         *newline;

   newline = (char *)calloc(LINELEN, sizeof(char));

   llen = strlen(iline);
   if (llen > LINELEN) {
      fprintf (stderr, "Input line longer than 10000 characters. ");
      fprintf (stderr, "Exiting\n");
      exit (1);
   }

   ii = 0;
   jj = 0;
   nword = 0;
   inword_f = 0;
   while (ii < llen) {
      if_sep = 0;
      for (kk=0;kk<nsep;kk++) {
         if (*(iline+ii) == *(sep+kk))
            if_sep = 1;
      }
      if (if_sep == 0) {
         *(newline+jj) = *(iline+ii);
         if (inword_f == 0) {
            inword_f = 1;
            strarray[nword++] = newline+jj;
         }
         jj++;
      }
      else {
         if (inword_f == 1) {
            inword_f = 0;
            *(newline+jj) = '\0';
            jj++;
         }
      }
      ii++;
   }

   return (nword);
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

double dotprod (double *vect1, double *vect2)
{
   double       scal;
   int          ii;

   scal = 0;

   for (ii=0;ii<3;ii++)
      scal+=vect1[ii]*vect2[ii];

   return scal;
}

void crossprod (double *vect1, double *vect2, double *crosvect)
{
   crosvect[0] = vect1[1]*vect2[2] - vect1[2]*vect2[1];

   crosvect[1] = vect1[2]*vect2[0] - vect1[0]*vect2[2];

   crosvect[2] = vect1[0]*vect2[1] - vect1[1]*vect2[0];

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


/* Here are the differential equations to be integrated. */
/* Non-aligned dipolar magnetic field (in 3D) : Global + Crust */

void equations(r, th, dthdr, surf_m, nsurf, alpha)
double r, alpha, th[], dthdr[], surf_m[];
int nsurf;
{
   double       Bd_r, Bd_th, Bd_ph;
   double       Bm_r, Bm_th, Bm_ph;
   double       md_r, md_th, md_ph;
   double       ms_r, ms_th, ms_ph;
   double       rs_r, rs_th, rs_ph;
   double       T, D;
   double       d, th_d;
   double      *m, *th_m, *ph_m; 
   double      *rs,*th_r, *ph_r;
   double       RS;
   int          ii;


   // Initializing Global dipole field
   RS = 10000.0;
   d = -1.0;
   th_d = alpha;

   //Reading the surface magnetic config file 
   rs = (double *)calloc(nsurf, sizeof(double));
   th_r = (double *)calloc(nsurf, sizeof(double));
   ph_r = (double *)calloc(nsurf, sizeof(double));
   m = (double *)calloc(nsurf, sizeof(double));
   th_m = (double *)calloc(nsurf, sizeof(double));
   ph_m = (double *)calloc(nsurf, sizeof(double));

   for (ii=0;ii<nsurf;ii++) {
      rs[ii] = surf_m[ii*6]*RS;
      th_r[ii] = surf_m[ii*6+1]/180.0*M_PI;
      ph_r[ii] = surf_m[ii*6+2]/180.0*M_PI;
      m[ii] = surf_m[ii*6+3]*d;
      th_m[ii] = surf_m[ii*6+4]/180.0*M_PI;
      ph_m[ii] = surf_m[ii*6+5]/180.0*M_PI;
   }

   // Global dipole  
   md_r  = d*(sin(th_d)*sin(th[1])*cos(th[2]) + cos(th_d)*cos(th[1]));
   md_th = d*(sin(th_d)*cos(th[1])*cos(th[2]) - cos(th_d)*sin(th[1]));
   md_ph = -d*sin(th_d)*sin(th[2]);
 
   Bd_r  = 2.0*md_r/pow(r,3.0);
   Bd_th = -md_th/pow(r,3.0);
   Bd_ph = -md_ph/pow(r,3.0);


   // Crust dipole
   Bm_r = 0.0;
   Bm_th = 0.0;
   Bm_ph = 0.0;

   for (ii=0;ii<nsurf;ii++) {
      ms_r  = m[ii]*(sin(th_m[ii])*sin(th[1])*cos(th[2]-ph_m[ii]) +
                         cos(th_m[ii])*cos(th[1]));

      ms_th = m[ii]*(sin(th_m[ii])*cos(th[1])*cos(th[2]-ph_m[ii]) -
                         cos(th_m[ii])*sin(th[1]));

      ms_ph = -m[ii]*sin(th_m[ii])*sin(th[2]-ph_m[ii]);

      rs_r  = rs[ii]*(sin(th_r[ii])*sin(th[1])*cos(th[2]-ph_r[ii]) +
                         cos(th_r[ii])*cos(th[1]));

      rs_th = rs[ii]*(sin(th_r[ii])*cos(th[1])*cos(th[2]-ph_r[ii]) -
                         cos(th_r[ii])*sin(th[1]));

      rs_ph = -rs[ii]*sin(th_r[ii])*sin(th[2]-ph_r[ii]);

      T = ms_r*r - (ms_r*rs_r + ms_th*rs_th + ms_ph*rs_ph);
      D = rs[ii]*rs[ii] + r*r - 2.0*r*rs_r;

      Bm_r+= -(3.0*T*rs_r - 3.0*T*r + D*ms_r)/pow(D,2.5);
      Bm_th+= -(3.0*T*rs_th + D*ms_th)/pow(D,2.5);
      Bm_ph+= -(3.0*T*rs_ph + D*ms_ph)/pow(D,2.5);
   }

   // Equation of field lines
   dthdr[1] = (Bd_th + Bm_th)/(r*(Bd_r + Bm_r));
   dthdr[2] = (Bd_ph + Bm_ph)/(r*(Bd_r + Bm_r)*sin(th[1]));

}

void coortrans (double *theta, double *phi, double alpha, double theta_in, 
                 double phi_in) 
{

  *theta = acos(cos(theta_in)*cos(alpha) - 
             sin(theta_in)*sin(alpha)*cos(phi_in));

  *phi = atan2(sin(theta_in)*sin(phi_in), 
                sin(theta_in)*cos(phi_in)*cos(alpha)+cos(theta_in)*sin(alpha));
}


void rk4(y,dydx,n,x,h,yout,derivs,surf_m,nsurf,alpha)
double y[],dydx[],x,h,yout[], surf_m[], alpha;
int nsurf;
void (*derivs)();
int n;
{ int i;
  double xh,hh,h6,*dym,*dyt,*yt,*vector();
  void free_vector();

  dym=vector(1,n);
  dyt=vector(1,n);
  yt=vector(1,n);
  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
  (*derivs)(xh,yt,dyt, surf_m, nsurf, alpha);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(xh,yt,dym, surf_m, nsurf, alpha);
  for (i=1;i<=n;i++) {
      yt[i]=y[i]+h*dym[i];
      dym[i] += dyt[i];    }
  (*derivs)(x+h,yt,dyt, surf_m, nsurf, alpha);
  for (i=1;i<=n;i++)
      yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
  free_vector(yt,1,n);
  free_vector(dyt,1,n);
  free_vector(dym,1,n);
}

#define PGROW -0.20
#define PSHRINK -0.25
#define FCOR 0.06666666
#define SAFETY 0.9
#define ERRCON 6.0e-4

void rkqc(y,dydx,n,x,htry,eps,yscal,hdid,hnext,derivs, surf_m, nsurf, alpha)
double y[],dydx[],*x,htry,eps,yscal[],*hdid,*hnext, surf_m[], alpha;
int nsurf;
void (*derivs)();
int n;
{
  int i;
  double xsav,hh,h,temp,errmax;
  double *dysav,*ysav,*ytemp,*vector();
  void rk4(),nerror(),free_vector();

  ysav=vector(1,n);
  dysav=vector(1,n);
  ytemp=vector(1,n);
  xsav=(*x);
  for (i=1;i<=n;i++) {
      ysav[i]=y[i];
      dysav[i]=dydx[i];  }
  h=htry;
  for (;;) {
      hh=0.5*h;
      rk4(ysav,dysav,n,xsav,hh,ytemp,derivs, surf_m, nsurf, alpha);
      *x=xsav+hh;
      (*derivs)(*x,ytemp,dydx, surf_m, nsurf, alpha);
      rk4(ytemp,dydx,n,*x,hh,y,derivs, surf_m, nsurf, alpha);
      *x=xsav+h;
      if (*x == xsav) nrerror("Step size too small in routine RKQC");
      rk4(ysav,dysav,n,xsav,h,ytemp,derivs, surf_m, nsurf, alpha);
      errmax=0.0;
      for (i=1;i<=n;i++) {
          ytemp[i]=y[i]-ytemp[i];
          temp=fabs(ytemp[i]/yscal[i]);
          if (errmax < temp) errmax=temp;  }
      errmax /= eps;
      if (errmax <= 1.0)  {
         *hdid=h;
         *hnext=(errmax > ERRCON ?
             SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
         break;  }
      h=SAFETY*h*exp(PSHRINK*log(errmax)); }
  for (i=1;i<=n;i++) y[i] += ytemp[i]*FCOR;
  free_vector(ytemp,1,n);
  free_vector(dysav,1,n);
  free_vector(ysav,1,n);
}

#define MAXSTP 10000
#define TINY 1.0e-30

void odeint(ystart, nvar, x1, x2, eps, h1, hmin, nok, nbad, derivs, surf_m, 
              nsurf, alpha)
double ystart[], surf_m[], x1, x2, eps, h1, hmin, alpha;
int nvar,*nok,*nbad, nsurf;
void (*derivs)();
{ int nstp,i;
  double xsav,x,hnext,hdid,h;
  double *yscal,*y,*dydx,*vector();
  void rkqc(), nrerror(), free_vector();

  yscal=vector(1,nvar);
  y=vector(1,nvar);
  dydx=vector(1,nvar);
  x=x1;
  h=(x2 > x1) ? fabs(h1) : -fabs(h1);
  *nok = (*nbad) = kount = 0;
  for (i=1;i<=nvar;i++) y[i]=ystart[i];
  if (kmax > 0) xsav=x-dxsav*2.0;
  for (nstp=1;nstp<=MAXSTP;nstp++)  {
     (*derivs)(x, y, dydx, surf_m, nsurf, alpha);
     for (i=1;i<=nvar;i++)
         yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
     if (kmax > 0)  {
        if (fabs(x-xsav) > fabs(dxsav))  {
           if (kount < kmax-1)  {
              xp[++kount]=x;
              for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
              xsav=x;
     }  }  }
     if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
     rkqc(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs,surf_m,nsurf,alpha);
     if (hdid == h) ++(*nok); else ++(*nbad);
     if ((x-x2)*(x2-x1) >= 0.0)  {
        for (i=1;i<=nvar;i++) ystart[i]=y[i];
        if (kmax)  {
           xp[++kount]=x;
           for (i=1;i<=nvar;i++) yp[i][kount]=y[i]; }
        free_vector(dydx,1,nvar);
        free_vector(y,1,nvar);
        free_vector(yscal,1,nvar);
        return; }
     if (fabs(hnext) <= hmin) nrerror("Step size too small in ODEINT");
     h=hnext;  }
  nrerror("Too many steps in the routine ODEINT");  
}
