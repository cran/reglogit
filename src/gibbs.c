// #include "rhelp.h"
#include <stdlib.h>
#include <assert.h>
#include <R.h>
#include <Rmath.h>

#define TOL sqrt(DOUBLE_EPS)

/*
 * sq:
 * 
 * calculate the square of x
 */

double sq(double x)
{
  return x*x;
}


/*
 * rinvgauss:
 *
 * Michael/Schucany/Haas Method for generating Inverse Gaussian
 * random variable, as given in Gentle's book on page 193
 */

double rinvgauss(const double mu, const double lambda)
{
  double u, y, x1, mu2, l2;

  y = sq(norm_rand());
  mu2 = sq(mu);
  l2 = 2*lambda;
  x1 = mu + mu2*y/l2 - (mu/l2)* sqrt(4*mu*lambda*y + mu2*sq(y));

  u = unif_rand();
  if(u <= mu/(mu + x1)) return x1;
  else return mu2/x1;
}


/*
 * rtnorm_reject:
 *
 * dummy function in place of the Robert (1995) algorithm 
 * based on proposals from the exponential, should have that 
 * mean < tau 
 */

double rtnorm_reject(double mean, double tau, double sd)
{
  double x, z, lambda;
  int cnt;

  /* Christian Robert's way */
  assert(mean < tau);
  tau = (tau - mean)/sd;

  /* optimal exponential rate parameter */
  lambda = 0.5*(tau + sqrt(sq(tau) + 4.0));

  /* do the rejection sampling */
  cnt = 0;
  do {
    z = rexp(1.0/lambda) + tau;
  } while (unif_rand() > exp(0.0-0.5*sq(z - lambda)));

  /* put x back on the right scale */
  x = z*sd + mean;

  assert(x > 0);
  return(x);

}


/*
 * draw_lambda_prior:
 *
 * draw Lambda from its (infinite mixture) prior
 */

double draw_lambda_prior(double *psii, int kmax)
{
  double lambda;
  int k;

  lambda = 0.0;
  for(k=0; k<=kmax; k++) {
    lambda += psii[k] * rexp(1.0);
  }
  
  return lambda;
}


/* 
 * draw_lambda_prior_R:
 *
 * R interface to lambda draw from the posterior conditional
 * by rejection sampling and proposals from the prior
 */

void draw_lambda_prior_R(int *n_in, double *psii_in, int *kmax_in, 
			 double *lambda_out)
{
  int i;

  GetRNGstate();

  for(i=0; i<*n_in; i++) 
    lambda_out[i] = draw_lambda_prior(psii_in, *kmax_in);

  PutRNGstate();  
}


/*
 * draw_lambda:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its 
 * full conditional -- uses proposals from the prior
 */

double draw_lambda(double lambda_old, double xbeta, 
		   double kappa, int kmax, int thin)
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
    lambda = draw_lambda_prior(psii, kmax);
    
    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta + 0.5*(1.0-kappa)*lambda;
    lp = pnorm(0.0, m, s, 0, 1);
    
    /* MH accept or reject */
    if(unif_rand() < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }

  /* possibly clean up psii */
  free(psii);

  return lambda_old;
}



/*
 * draw_lambda_zz_MH:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its 
 * full conditional under the z=0 model -- uses proposals 
 * from the prior
 */

double draw_lambda_zz_MH(double lambda_old, double xbeta, 
			       double kappa, int kmax, int thin)
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
    lambda = draw_lambda_prior(psii, kmax);
    
    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta + 0.5*(1.0-kappa)*lambda;
    lp = dnorm(0.0, m, s, 1);

    /* MH accept or reject */
    if(unif_rand() < exp(lp - lpold)) {
      lambda_old = lambda;
      lpold = lp;
    }
  }

  /* possibly clean up psii */
  free(psii);

  return lambda_old;
}


/*
 * draw_lambda_zz_slice:
 *
 * Metropolis-Hastings algorithm for drawing lambda from its 
 * full conditional under the z=0 model -- uses slice sampling
 */

double draw_lambda_zz_slice(double lambda_old, double xbeta, 
			       double kappa, int kmax)
{
  int k, cnt;
  double lambda, p, pold, m, s, u;
  double *psii;

  /* calculate the probability og the previous lambda */
  s = sqrt(lambda_old);
  m = xbeta + 0.5*(1.0-kappa)*lambda_old;
  pold = dnorm(0.0, m, s, 0);

  /* draw the uniform latent variable */
  u = runif(0, pold);

  /* allocate psii */
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] = 2.0/((0.5+k)*(kappa-0.5+k));
  
  cnt = 0;
  do {

    /* propose a new lambda from the prior */
    lambda = draw_lambda_prior(psii, kmax);
    
    /* calculate the probability of the propsed lambda */
    s = sqrt(lambda);
    m = xbeta + 0.5*(1.0-kappa)*lambda;
    p = dnorm(0.0, m, s, 0);

    /* sanity print for poor rejection rate */
    cnt++;
    if(cnt % 1000000 == 0) 
      warning("slice=%d, xbeta=%g, lambda=%g, u=%g, kappa=%g\n", 
	      cnt, xbeta, lambda_old, u, kappa);

  } while(p <= u);


  /* clean up */
  free(psii);

  return lambda;
}


/*
 * draw_lambda_zz_vaduva:
 *
 * Vaduva's generalization of the rejection method for
 * sampling lambda_i in the z=0 version of the model
 */

double draw_lambda_zz_vaduva(double xbeta, double kappa, int kmax)
{
  int k, cnt;
  double lambda, INv, INm, y;
  double *psii;

  /* lambda parmater for IN */
  INv = 0.25*(2.0-kappa)*(2.0-kappa) + kappa - 1.0;
  INm = sqrt(INv)/fabs(xbeta);

  /* allocate psii */
  psii = (double*) malloc(sizeof(double) * (kmax+1));
  for(k=0; k<=kmax; k++) psii[k] = 2.0/((1.0+k)*(kappa+1.0+k));
  
  cnt = 0;
  do {
    
    /* x draw in the notes */
    lambda = 1.0/rinvgauss(INm, INv);

    /* y draw in the notes */
    y = draw_lambda_prior(psii, kmax);

    /* sanity print for poor rejection rate */
    cnt++;
    if(cnt % 100000 == 0) 
      warning("vaduva=%d, x=%g, y=%g, xbeta=%g, kappa=%g, INv=%g\n",
	      cnt, lambda, y, xbeta, kappa, INv);

  } while(lambda < y);

  
  /* clean up */
  free(psii);

  return lambda;
}


/* 
 * draw_lambda_R:
 *
 * R interface to lambda draw from the posterior conditional
 * by rejection sampling and proposals from the prior
 */

void draw_lambda_R(int *n_in, double *xbeta_in, double *kappa_in, 
		   int *kl_in, int *kmax_in, int *zzero_in, 
		   int *thin_in, int *tl_in, int *method_in, 
		   double *lambda_inout)
{
  int i, thin;
  double kappa;

  GetRNGstate();

  /* initial values */
  kappa = *kappa_in;
  thin = *thin_in;

  for(i=0; i<*n_in; i++) {

    /* depends on whether we are vectorizing or not */
    if(*kl_in > 1) kappa = kappa_in[i];
    if(*tl_in > 1) thin = thin_in[i];

    /* the actual draw */
    if(!*zzero_in)
      lambda_inout[i] = draw_lambda(lambda_inout[i], xbeta_in[i], 
				    kappa, *kmax_in, thin);
    else {
      switch(*method_in) {
      case 1: 
	lambda_inout[i] = draw_lambda_zz_vaduva(xbeta_in[i], kappa, 
						*kmax_in); 
	break;
      case 2:
        lambda_inout[i] = 
	  draw_lambda_zz_slice(lambda_inout[i], xbeta_in[i], kappa, 
			       *kmax_in);
	break;
      case 3:
	lambda_inout[i] = 
	  draw_lambda_zz_MH(lambda_inout[i], xbeta_in[i], kappa, 
			    *kmax_in, thin);
	break;
      default: error("no default defined");
      }
    }
  }

  PutRNGstate();
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
  int n, i;
  double xbeta, madj, lambda_sqrt;

  /* copy scalars */
  n = *n_in;
 
  /* calculate mean multiplicative adjustment due to kappa */
  madj = 0.5*(1.0 - *kappa_in);

  /* get the RNG state */
  GetRNGstate();

  /* loop over rows of X */
  for(i=0; i<n; i++) {

    /* depends on whether we are vectorizing or not */
    if(*kl_in > 1) madj = 0.5*(1.0 - kappa_in[i]);
    
    /* calculate the mean and variance of the normal */
    xbeta = xbeta_in[i] + madj * lambda_in[i];
    lambda_sqrt = sqrt(lambda_in[i]);

    /* draw until we get one in the right half-plane */
    if(xbeta >= 0) {
      do { z_out[i] = rnorm(xbeta, lambda_sqrt);
      } while (z_out[i] < 0.0);
    } else { 
      z_out[i] = rtnorm_reject(xbeta, 0.0, lambda_sqrt);
      assert(z_out[i] > 0.0);
    }
  }

  /* put the RNG state back */
  PutRNGstate();
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
