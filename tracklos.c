//  To read the open field line boundary on the emission region and polar cap
//   and for a given geometry produce the line of sight tracks.

//  To compile: gcc tracklos.c -o ~/bin/tracklos -lm
   
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>


# define WORDSIZE (100) // The size of file names.
# define LINELEN (1000) // The word size of line in files.
# define ZERO_lf (0.0000001) // round off zero error.


#define USAGE "\
   usage: tracklos infield inopen trkfile alpha beta nlos\n"


void cmdlineError (char *message);

void nrerror(char error_text[]);

double *vector(int nl, int nh);

int splitstr (char *iline, char *strarray[], char *sep, int nsep);

double findedge (double x1, double y1, double x2, double y2, double val);


int kmax=0,kount=0;
double *xp=0,**yp=0,dxsav=0;

/* The main program. */
int main(int argc, char *argv[])
{ 
   FILE        *fpt;
   char        *magfile, *openfile, *trackfile;
   char        *readline, *readval[20];
   void         odeint(), equations();
   double      *th, *dthdr;
   double      *outfield, *surf_m, *los_tr;
   double       alpha, beta, cen_thet;
   double       RS, r, delr, r1, r2;
   double       phi_i, phi_f, phi, del_phi;
   int          nsurf, nfield, nlos;
   int          ii, jj, nok, nbad;


   magfile = (char *)calloc(WORDSIZE, sizeof(char));
   openfile = (char *)calloc(WORDSIZE, sizeof(char));
   trackfile = (char *)calloc(WORDSIZE, sizeof(char));
   readline = (char *)calloc(LINELEN, sizeof(char));

   if (argc < 7)
      cmdlineError (USAGE);

   sprintf (magfile, "%s", argv[1]);
   sprintf (openfile, "%s", argv[2]);
   sprintf (trackfile, "%s", argv[3]);
   alpha = atof(argv[4])/180.0*M_PI;
   beta = atof(argv[5])/180.0*M_PI;
   nlos = atoi(argv[6]);

   //Reading the surface magnetic config file 
   if ((fpt = fopen(magfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", magfile);
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
   fclose(fpt);

   //Reading the open field boundary
   if ((fpt = fopen(openfile,"r")) == NULL){
      fprintf(stderr, "Error in opening %s file\n", openfile);
      exit(1);
   }

   fgets(readline, LINELEN, fpt);

   nfield = 0;
   while (fgets(readline, LINELEN, fpt)) 
      nfield++;

   fclose(fpt);

   outfield = (double *)calloc(6*nfield+1, sizeof(double)); 

   fpt = fopen(openfile,"r");
   fgets(readline, LINELEN, fpt);

   nfield = 0;
   while (fgets(readline, LINELEN, fpt)) {
      splitstr (readline, readval, " ", 1);

      outfield[6*nfield] = atof(readval[0]);
      outfield[6*nfield+1] = atof(readval[1]);
      outfield[6*nfield+2] = atof(readval[2]);
      outfield[6*nfield+3] = atof(readval[6]);
      outfield[6*nfield+4] = atof(readval[7]);
      outfield[6*nfield+5] = atof(readval[8]);
      nfield++;
   }
   fclose(fpt);

   //The boundaries of line of sight along open field lines.
   phi_i = 0.0;
   for (ii=0;ii<nfield/2-1;ii++) {
      if (alpha+beta > outfield[6*ii+1] && alpha+beta < outfield[6*(ii+1)+1]) 
         phi_i = findedge (outfield[6*ii+2], outfield[6*ii+1], 
                  outfield[6*(ii+1)+2], outfield[6*(ii+1)+1], alpha+beta);
   }

   phi_f = 0.0;
   for (ii=nfield/2;ii<nfield-1;ii++) {
      if (alpha+beta <= outfield[6*ii+1] && alpha+beta > outfield[6*(ii+1)+1])
         phi_f = findedge (outfield[6*ii+2], outfield[6*ii+1], 
                  outfield[6*(ii+1)+2], outfield[6*(ii+1)+1], alpha+beta);
   }

   fprintf(stderr,"%lf  %lf  %lf\n", phi_i, phi_f, alpha+beta);

   if (fabs(phi_i) < ZERO_lf && fabs(phi_f) < ZERO_lf) {
      phi_i = -M_PI;
      phi_f = +M_PI;
   }

   // Generating line of sight tracks on the surface.
   los_tr = (double *)calloc(3*(nlos+1), sizeof(double));

   th = vector(1,2); 
   dthdr = vector(1,2);
   RS = 10000;  // neutron star radius in meters.
   r1 = 30*RS; 
   r2 = RS;
   delr = 0.01*RS;
   del_phi = (phi_f-phi_i)/nlos;

   phi = phi_i;
   ii = 0;
   while (phi <= phi_f) {
      fprintf(stderr, "%lf\r", phi*180/M_PI);
      los_tr[3*ii] = phi;
      r = r1;
      th[1] = alpha+beta;
      th[2] = phi;

      while (r > r2) {
         odeint(th, 2, r+delr, r, 1.0e-3, 0.05*delr, 1.0e-6, &nok, &nbad,
                 equations, surf_m, nsurf, alpha);

         r-=delr;
      }
      los_tr[3*ii+1] = th[1];
      los_tr[3*ii+2] = th[2];
      phi+=del_phi;
      ii++;
   }

   nlos = ii;

   //Estimating the axis inclination
   r = r1;
   th[1] = alpha;
   th[2] = 0.0;

   while (r > r2) {
      odeint(th, 2, r+delr, r, 1.0e-3, 0.05*delr, 1.0e-6, &nok, &nbad,
              equations, surf_m, nsurf, alpha);

      r-=delr;
   }
   cen_thet = th[1];

   //Writing out track in file
   fpt = fopen(trackfile, "w");

   fprintf(fpt, "#nlos %d  axis %lf\n", nlos, cen_thet);

   for (ii=0;ii<nlos;ii++)
      fprintf (fpt, "%lf  %lf  %lf  %lf  %lf  %lf\n", los_tr[3*ii], los_tr[3*ii+1], los_tr[3*ii+2], RS*sin(los_tr[3*ii+1])*cos(los_tr[3*ii+2]), RS*sin(los_tr[3*ii+1])*sin(los_tr[3*ii+2]), RS*cos(los_tr[3*ii+1]));

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


double findedge (double x1, double y1, double x2, double y2, double val)
{
   double       m_ln, b_ln;

   m_ln = (y2-y1)/(x2-x1);
   b_ln = y2-m_ln*x2;

   return (val-b_ln)/m_ln;
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
   d = 1.0;
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
      m[ii] = surf_m[ii*6+3];
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
