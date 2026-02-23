#ifndef PCdrift_H
#define PCdrift_H

void rec2pol (double x, double y, double *r, double *theta);

void findbase(double *inarray, int nsize, int basefact, int *spos, int *epos,
               double *minstat);

void meanrms(double *inarray, int arrsize, double *arrstat);

double urandom(long *idum);

void sparkconfig (double *th_sprk_u, double *th_sprk_d, double *x_sprk,
                    double *y_sprk, double *a_sprk, double *theta_sp,
                    double h_sprk, double h_drft, double a_cap, double b_cap,
                    double th_cap, double co_angl, double x_cent,
                    double y_cent, int N_trk, int *N_up, int *N_dn, 
		    int trk_max, int *ncap);

void unitvect (double *invect, double *transvect);

void pospolcap (double *los_tr, double phi, double th_c, double ph_c,
                  double max_phi, double min_phi, double *x_out, double *y_out,
                  int nlos);

void coortrans (double x_in, double y_in, double theta, double *x_out,
                   double *y_out);
#endif

