/****************************************************************************
 *
 * Regularized logistic regression
 * Copyright (C) 2011, The University of Chicago
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301  USA
 *
 * Questions? Contact Robert B. Gramacy (rbgramacy@chicagobooth.edu)
 *
 ****************************************************************************/


#include <stdlib.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>
#ifdef _OPENMP
  #include <omp.h>
#endif 
#include "rand_draws.h"

#define TOL sqrt(DOUBLE_EPS)


/*
 * draw_lambda_prior:
 *
 * draw Lambda from its (infinite mixture) prior
 */

double draw_lambda_prior(double *psii, int kmax, rk_state *state)
{
  double lambda;
  int k;

  lambda = 0.0;
  for(k=0; k<=kmax; k++) {
    lambda += psii[k] * expo_rand(state);
  }
  
  return lambda;
}


/*
 * draw_lambda:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its 
 * full conditional -- uses proposals from the prior
 */

double draw_lambda(double lambda_old, double xbeta, 
		   double kappa, int kmax, int thin, rk_state *state)
{
  int t, k;
  double lambda, lp, lpold, m, s;
  double *psii;

  /* calculate the probability og the previous lambda */
  s = sqrt(lambda_old);
  m = xbeta + 0.5*(1.0-kappa)*lambda_old;
  lpold = pnorm(0.0, m, s, 0, 1);

  /* allocate psii */
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] =  2.0/((1.0+k)*(kappa+k));
    
  /* thinning is essential when kappa is large */
  for(t=0; t<thin; t++) {

    /* propose a new lambda from the prior */
    lambda = draw_lambda_prior(psii, kmax, state);
    
    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta + 0.5*(1.0-kappa)*lambda;
    lp = pnorm(0.0, m, s, 0, 1);
    
    /* MH accept or reject */
    if(runi(state) < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }

  /* possibly clean up psii */
  free(psii);

  return lambda_old;
}



/*
 * draw_lambda_zz:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its 
 * full conditional under the z=0 model -- uses proposals 
 * from the prior
 */

double draw_lambda_zz(double lambda_old, double xbeta, 
			       double kappa, int kmax, int thin, rk_state *state)
{
  int t, k;
  double lambda, lp, lpold, m, s;
  double *psii;

  /* calculate the probability og the previous lambda */
  s = sqrt(lambda_old);
  m = xbeta + 0.5*(1.0-kappa)*lambda_old;
  lpold = dnorm(0.0, m, s, 1);

  /* allocate psii */
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] = 2.0/((0.5+k)*(kappa-0.5+k));
    
  /* thinning is essential when kappa is large */
  for(t=0; t<thin; t++) {

    /* propose a new lambda from the prior */
    lambda = draw_lambda_prior(psii, kmax, state);
    
    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta + 0.5*(1.0-kappa)*lambda;
    lp = dnorm(0.0, m, s, 1);

    /* MH accept or reject */
    if(runi(state) < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }

  /* possibly clean up psii */
  free(psii);

  return lambda_old;
}


/* 
 * draw_lambda_R:
 *
 * R interface to lambda draw from the posterior conditional
 * by rejection sampling and proposals from the prior
 */

void draw_lambda_R(int *n_in, double *xbeta_in, double *kappa_in, 
		   int *kl_in, int *kmax_in, int *zzero_in, 
		   int *thin_in, int *tl_in, double *lambda_inout)
{
  int thin;
  double kappa;

  /* initial values */
  kappa = *kappa_in;
  thin = *thin_in;

  #ifdef _OPENMP
  #pragma omp parallel
  {
  int i, start, step;
  double aux[2];
  double lambda_sqrt, xbeta;
  rk_state *state;
  start = omp_get_thread_num();
  step = omp_get_max_threads();
#else 
  int i, start, step;
  double aux[2];
  double lambda_sqrt, xbeta;
  rk_state *state;
  start = 0;
  step = 1;
#endif

  state = states[start];

  for(i=start; i<*n_in; i += step) {

    /* depends on whether we are vectorizing or not */
    if(*kl_in > 1) kappa = kappa_in[i];
    if(*tl_in > 1) thin = thin_in[i];

    /* the actual draw */
    if(!*zzero_in)
      lambda_inout[i] = draw_lambda(lambda_inout[i], xbeta_in[i], 
				    kappa, *kmax_in, thin, state);
    else {
      lambda_inout[i] = draw_lambda_zz(lambda_inout[i], xbeta_in[i], kappa, 
			    *kmax_in, thin, state);
    }
  }
#ifdef _OPENMP
  }
#endif

}


/* 
 * draw_z_R:
 *
 * C-side for R function draw.z, which draws from the full
 * conditional of the latent z-values given y, beta, and 
 * the latent lambdas
 */

void draw_z_R(int *n_in, double *xbeta_in, double *beta_in, 
	      double *lambda_in, double *kappa_in, int *kl_in, 
	      double *z_out)
{
  int n;
  double madj;

  /* copy scalars */
  n = *n_in;
 
  /* calculate mean multiplicative adjustment due to kappa */
  madj = 0.5*(1.0 - *kappa_in);

#ifdef _OPENMP
  #pragma omp parallel
  {
  int i, start, step;
  double aux[2];
  double lambda_sqrt, xbeta;
  rk_state *state;
  start = omp_get_thread_num();
  step = omp_get_max_threads();
#else 
  int i, start, step;
  double aux[2];
  double lambda_sqrt, xbeta;
  rk_state *state;
  start = 0;
  step = 1;
#endif

  state = states[start];

  /* loop over rows of X */
  for(i=start; i<n; i += step) {

    /* depends on whether we are vectorizing or not */
    if(*kl_in > 1) madj = 0.5*(1.0 - kappa_in[i]);
    
    /* calculate the mean and variance of the normal */
    xbeta = xbeta_in[i] + madj * lambda_in[i];
    lambda_sqrt = sqrt(lambda_in[i]);

    /* draw until we get one in the right half-plane */
    if(xbeta >= 0) {
      do { 
        rnor(aux, state);
        z_out[i] = aux[0]*lambda_sqrt + xbeta;
        if(z_out[i] >= 0) break;
        z_out[i] = aux[1]*lambda_sqrt + xbeta;

      } while (z_out[i] < 0.0);
    } else { 
      z_out[i] = rtnorm_reject(xbeta, 0.0, lambda_sqrt, state);
      assert(z_out[i] > 0.0);
    }

    // #ifdef _OPENMP
    //     #pragma omp master
    // #endif
    // fprintf(stdout, "i=%d\n", i);
  }
#ifdef _OPENMP
  }
#endif

}


/* 
 * draw_omega_R:
 *
 * C-side for R function draw.omega, which draws from the full
 * conditional of the latent omega-values given beta, 
 * sigma, and nu
 */

void draw_omega_R(int *m_in, double *beta_in, double *sigma_in, 
		  double *nu_in, double *omega_out)
{
  int j;
  double mu;

  GetRNGstate();

  /* for each latent omega */
  for(j=0; j<*m_in; j++) {
    
    /* calculate the mean of the inverse gaussian */
    mu = *nu_in * sigma_in[j] / fabs(beta_in[j]);

    /* sample from the inverse gaussian and invert */
    omega_out[j] = 1.0/rinvgauss(mu, 1.0);
  }

  PutRNGstate();
}


/*
 * rinvgauss_R:
 *
 * R interface to R_invgauss
 */

void rinvgauss_R(int* n_in, double *mu_in, double *lambda_in, 
		 double *x_out)
{
  int i;

  GetRNGstate();

  for(i=0; i<*n_in; i++) 
    x_out[i] = rinvgauss(*mu_in, *lambda_in);

  PutRNGstate();
}
