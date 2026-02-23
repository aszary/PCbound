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


void crossprod (double *vect1, double *vect2, double *crosvect)
{
   crosvect[0] = vect1[1]*vect2[2] - vect1[2]*vect2[1];

   crosvect[1] = vect1[2]*vect2[0] - vect1[0]*vect2[2];

   crosvect[2] = vect1[0]*vect2[1] - vect1[1]*vect2[0];

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


void sparkconfig (double *th_sprk_u, double *th_sprk_d, double *x_sprk,
                    double *y_sprk, double *a_sprk, double *theta_sp,
                    double h_sprk, double h_drft, double a_cap, double b_cap,
                    double th_cap, double co_angl, double x_cent,
                    double y_cent, int N_trk, int *N_up, int *N_dn, 
		    int trk_max, int *ncap)
{
   double      trk_a_ul[N_trk], trk_b_ul[N_trk];
   double      trk_a_ut[N_trk], trk_b_ut[N_trk];
   double      trk_a_dl[N_trk], trk_b_dl[N_trk];
   double      trk_a_dt[N_trk], trk_b_dt[N_trk];
   double      a_in, a_out, b_in, b_out, e_trk, b_sprk;
   double      x_val, y_val, x_in, y_in, x_trns, y_trns, el_val;
   double      theta_c, theta_fl, theta_fh;
   double      a_trk, b_trk, rad_fh, rad_fl;
   int         indx_u, indx_d, indx_us;
   int         ii, jj, ncnt, ncell, nsprk_l, nsprk_h;


   a_out = a_cap;
   a_in = a_out - h_sprk;

   b_out = b_cap;
   b_in = b_out - b_cap/a_cap*h_sprk;

   ncnt = 0;

   indx_u = 0;
   indx_d = 0;

   indx_us = 0;

   for (ii=0;ii<N_trk;ii++) {

      a_trk = 0.5*(a_out+a_in);
      b_trk = 0.5*(b_out+b_in);

      //Upper half, clockwise track
      for (jj=0;jj<N_up[ii];jj++) {

         x_in = a_trk*cos(th_sprk_u[indx_u+jj] - th_cap);
         y_in = b_trk*sin(th_sprk_u[indx_u+jj] - th_cap);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;
         a_sprk[ncnt] = h_sprk/2;

         //Spark at leading edge
         if (M_PI - co_angl - th_sprk_u[indx_u+jj] <= theta_sp[ii]/2) {

            theta_fh = 0.5*(M_PI-co_angl-th_sprk_u[indx_u+jj]) + theta_sp[ii]/4;

            a_sprk[ncnt] = a_trk*sin(theta_fh);

            trk_a_ul[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

            x_in = trk_a_ul[ii]*cos(M_PI - co_angl - theta_fh - th_cap);
            y_in = trk_a_ul[ii]*b_trk/a_trk*sin(M_PI-co_angl-theta_fh-th_cap);

            coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

            x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

         //Spark at trailing edge
         if (th_sprk_u[indx_u+jj] <= theta_sp[ii]/2 - co_angl) {

            theta_fh = 0.5*(th_sprk_u[indx_u+jj] + theta_sp[ii]/2 + co_angl);

            a_sprk[ncnt] = a_trk*sin(theta_fh);

            trk_a_ut[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

            x_in = trk_a_ut[ii]*cos(theta_fh - th_cap - co_angl);
            y_in = trk_a_ut[ii]*b_trk/a_trk*sin(theta_fh - th_cap - co_angl);

            coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

            x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

         ncnt++;
      }


      //lower half, anti-clockwise track
      for (jj=0;jj<N_dn[ii];jj++) {

         x_in = a_trk*cos(th_sprk_d[indx_d+jj] - th_cap);
         y_in = b_trk*sin(th_sprk_d[indx_d+jj] - th_cap);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;
         a_sprk[ncnt] = h_sprk/2;

         //Spark at leading edge
         if (th_sprk_d[indx_d+jj] - M_PI + co_angl <= theta_sp[ii]/2) {

            theta_fl = 0.5*(th_sprk_d[indx_d+jj]-M_PI+co_angl)+theta_sp[ii]/4;

            a_sprk[ncnt] = 0.5*(2*a_trk*sin(theta_fl));

            trk_a_dl[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

            x_in = trk_a_dl[ii]*cos(M_PI - co_angl + theta_fl - th_cap);
            y_in = trk_a_dl[ii]*b_trk/a_trk*sin(M_PI-co_angl+theta_fl-th_cap);

            coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

            x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

         //Spark at trailing edge
         if (2*M_PI-th_sprk_d[indx_d+jj] - co_angl <= theta_sp[ii]/2) {

            theta_fl = 0.5*(2*M_PI-th_sprk_d[indx_d+jj]-co_angl)+theta_sp[ii]/4;

            a_sprk[ncnt] = 0.5*(2*a_trk*sin(theta_fl));

            trk_a_dt[ii] = a_trk + h_sprk/2 - a_sprk[ncnt];

            x_in = trk_a_dt[ii]*cos(2*M_PI - theta_fl - th_cap - co_angl);
            y_in = trk_a_dt[ii]*b_trk/a_trk*sin(2*M_PI-theta_fl-th_cap-co_angl);

            coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

            x_sprk[ncnt]+=x_cent;
            y_sprk[ncnt]+=y_cent;
         }

         ncnt++;
      }

/*
      //Plugging gap in the beginning
      if (th_sprk_d[indx_d] - th_sprk_u[indx_u] >= theta_sp[ii]) {

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(th_sprk_d[indx_d]-th_sprk_u[indx_u])) - h_sprk);

         e_trk = a_trk + h_sprk/2 - a_sprk[ncnt];

         x_in = e_trk*cos(M_PI - th_cap - co_angl);
         y_in = e_trk*b_trk/a_trk*sin(M_PI - th_cap - co_angl);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+=x_cent;
         y_sprk[ncnt]+=y_cent;


         ncnt++;
      }

      //Plugging gap in the end
      theta_fh = th_sprk_u[indx_u+N_up[ii]-1];
      theta_fl = 2*M_PI - th_sprk_d[indx_d+N_dn[ii]-1];

      if (theta_fl+theta_fh >= theta_sp[ii] &&
                                 theta_fl+theta_fh <= 2*theta_sp[ii]) {

         a_sprk[ncnt] = 0.5*(2*a_trk*sin(0.5*(theta_fl+theta_fh)) - h_sprk);

         e_trk = a_trk + h_sprk/2 - a_sprk[ncnt];

         x_in = e_trk*cos(2*M_PI - th_cap - co_angl);
         y_in = e_trk*b_trk/a_trk*sin(2*M_PI - th_cap - co_angl);

         coortrans (x_in, y_in, -th_cap, &x_sprk[ncnt], &y_sprk[ncnt]);

         x_sprk[ncnt]+= x_cent;
         y_sprk[ncnt]+=y_cent;

         ncnt++;
      }
*/
      indx_u+=trk_max;
      indx_d+=trk_max;

      a_out-=h_sprk;
      a_in-=h_sprk;

      b_out-=b_cap/a_cap*h_sprk;
      b_in-=b_cap/a_cap*h_sprk;
   }

   //Core Spark
   x_sprk[ncnt] = x_cent;
   y_sprk[ncnt] = y_cent;
   a_sprk[ncnt] = a_cap - N_trk*h_sprk - 1.5*h_drft;
   ncnt++;

   *ncap = ncnt;
}


void pospolcap (double *los_tr, double phi, double th_c, double phi_c, 
		  double max_phi, double min_phi, double *x_out, double *y_out,
		  int nlos)
{
   double       inarr[3], rotaxis[3], trnsvect[3];
   double       RS, th_cp, phi_cp;
   double       x_in, y_in, x_cent, y_cent, z_cent;
   int          ii;


   for (ii=0;ii<nlos-1;ii++) { 

      if (phi >= los_tr[3*ii] && phi < los_tr[3*(ii+1)]) {
	 
         phi_cp = findbet (los_tr[3*ii], los_tr[3*(ii+1)], phi, 
			          los_tr[3*ii+2], los_tr[3*(ii+1)+2]);

	 th_cp = findedge (los_tr[3*ii+1], los_tr[3*ii+2], los_tr[3*(ii+1)+1],
			     los_tr[3*(ii+1)+2], phi_cp);
      }
   }


   RS = 10000.0; //Radius of neutron star
 
   x_cent = RS*sin(th_c)*cos(phi_c);
   y_cent = RS*sin(th_c)*sin(phi_c);
   z_cent = RS*cos(th_c);

   inarr[0] = sin(th_c)*cos(max_phi) - sin(th_c)*cos(min_phi);
   inarr[1] = sin(th_c)*sin(max_phi) - sin(th_c)*sin(min_phi);
   inarr[2] = 0.0;

   unitvect (inarr, rotaxis);

   inarr[0] = RS*sin(th_cp)*cos(phi_cp) - x_cent;
   inarr[1] = RS*sin(th_cp)*sin(phi_cp) - y_cent;
   inarr[2] = RS*cos(th_cp) - z_cent;

   rotvect (inarr, rotaxis, trnsvect, -th_c);

  *x_out = trnsvect[0] + x_cent;
  *y_out = trnsvect[1] + y_cent;
}


double sparkamp (double *x_sprk, double *y_sprk, double *a_sprk, double a_cap, 
		   double b_cap, double th_cap, double x_in, double y_in, 
		   int nsprk)
{
   double       del_x, del_y, x_trns, y_trns;
   double       b_sprk, el_val, sig_x, sig_y, sprk_int;
   int          ii, if_chk;


   if_chk = 0;
   ii = 0;
 
   sprk_int = 0.0;

   while (if_chk < 1) {

      del_x = x_in - x_sprk[ii];
      del_y = y_in - y_sprk[ii];

      coortrans(del_x, del_y, th_cap, &x_trns, &y_trns);

      b_sprk = a_sprk[ii]*b_cap/a_cap;

      el_val = sqrt(x_trns*x_trns/a_sprk[ii]/a_sprk[ii] +
                            y_trns*y_trns/b_sprk/b_sprk);

      if (el_val < 1) {
	    sig_x = a_sprk[ii]/1.75;
            sig_y = b_sprk/1.75;

            sprk_int+=exp(-x_trns*x_trns/(2.0*sig_x*sig_x) - 
                             y_trns*y_trns/(2.0*sig_y*sig_y));

//   fprintf (stderr, "x_in %lf   y_in  %lf   ii %d   x_sprk[ii]  %lf   y_sprk[ii] %lf  x_trns %lf  y_trns %lf\n", x_in, y_in, ii, x_sprk[ii], y_sprk[ii], x_trns, y_trns);
//   getchar();
	 if_chk = 1;
      }
      
      if (ii == nsprk-1)
         if_chk = 1;

      ii++;      
   }

   return sprk_int;
}
