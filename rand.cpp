#include "rand.h"
#include <stdlib.h>
#include <math.h>
				// should be in math.h but fuck it
# define M_PI               3.14159265358979323846  /* pi */

/* 
   exponential distribution, specified mean
 */
double exponl_distn(double mean)
{
				// 1 - unif allows log(1) but not log(0)
				// supposedly
  return -mean*log(1 - unif_distn());
}

/*
  poisson deviate, from numerical recipes in C pp. 294ff
*/
double poidev(double mean)
{
  double gammln(double xx);
  static double sq, alxm, g, oldm=(-1.0);
  double em, t, y;

  if (mean < 12.0)
  {
    if (mean != oldm)
    {
      oldm = mean;
      g=exp(-mean);
    }
    em = -1;
    t=1.0;
    do
    {
      ++em;
      t *= unif_distn();
    }  while (t > g);
  }
  else
  {
    if (mean != oldm)
    {
      oldm=mean;
      sq=sqrt(2.0*mean);
      alxm=log(mean);
      g=mean*alxm-gammln(mean+1.0);
    }
    do
    {
      do
      {
	y = tan(M_PI*unif_distn());
	em = sq*y+mean;
      }  while (em < 0.0);
      em = floor(em);
      t=0.9 * (1.0 + y*y) * exp(em*alxm-gammln(em+1.0)-g);
    }  while (unif_distn() > t);
  }
  return em;
}

/* gammln, used by poidev
   numerical recipes in C, p. 214
*/

double gammln(double xx)
{
  double x,y,tmp,ser;
  static double cof[6] = 
    { 76.18009172947146, -86.50532032941677,
      24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5 };
  int j;

  y = x = xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser = 1.000000000190015;
  for (j=0; j<=5; j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005 * ser/x);
}

/*
  normal distribution, specified mean + std. dev.
  based on numerical recipes in C, pp. 289ff
*/
double norml_distn(double mean, double stdev)
{
  static int iset=0;
  static double gset;
  double fac,r,v1,v2;

  if (iset == 0) {
    do {
      v1=2.0*unif_distn()-1.0;
      v2=2.0*unif_distn()-1.0;
      r=v1*v1+v2*v2;
    } while (r >= 1.0);
    fac=sqrt(-2.0*log(r)/r);

    gset=v1*fac;
    iset=1;
    return v2*fac*stdev + mean;
  } else {
    iset=0;
    return gset*stdev + mean;
  }
}

/*
  binomial dist'n
  numerical recipes pp. 295,296
*/

#include <math.h> 
#define PI 3.141592654 
double bnldev(double pp, int n/*, long *idum*/) 
{
  double gammln(double xx); 
  //  float ran1(long *idum); 
  int j; 
  static int nold=(-1); 
  double am,em,g,angle,p,bnl,sq,t,y; 
  static double pold=(-1.0),pc,plog,pclog,en,oldg; 
  p=(pp <= 0.5 ? pp : 1.0-pp); 
  //The binomial distribution is invariant under changing pp to 1-pp,
  // if we also change the answer to n minus itself; we'll remember to 
  // do this below.
  am=n*p; 
  // This is the mean of the deviate to be produced. 
  if (n < 25) { 
    // Use the direct method while n is not too large. This can require up
    // to 25 calls to ran1. 
    bnl=0.0; 
    for (j=1;j<=n;j++) 
      //if (ran1(idum) < p) 
      if (unif_distn() < p)
	++bnl; 
  } 
  else if (am < 1.0)
  { 
    //If fewer than one event is expected out of 25 or more trials, 
    // then the distribution is quite accurately Poisson. Use direct
    // Poisson method. 
    g=exp(-am); 
    t=1.0; 
    for (j=0;j<=n;j++) 
    { 
      //t *= ran1(idum); 
      t *= unif_distn();
      if (t < g) break; 
    } 
    bnl=(j <= n ? j : n); 
  } 
  else 
  { 
    // Use the rejection method.
    if (n != nold) 
    { 
      // If n has changed, then compute useful quantities. 
      en=n; 
      oldg=gammln(en+1.0); 
      nold=n; 
    } 
    if (p != pold) 
    { 
      // If p has changed, then compute useful quantities. 
      pc=1.0-p; 
      plog=log(p); 
      pclog=log(pc); 
      pold=p;
    } 
    sq=sqrt(2.0*am*pc); 
    // The following code should by now seem familiar: 
    // rejection method with a Lorentzian comparison function. 
    do { 
      do { 
	//angle=PI*ran1(idum); 
	angle = PI*unif_distn();
	y=tan(angle); 
	em=sq*y+am; 
      } 
      while (em < 0.0 || em >= (en+1.0)); 
      // Reject. 
      em=floor(em); 
      // Trick for integer-valued distribution. 
      t=1.2*sq*(1.0+y*y)*exp(oldg-gammln(em+1.0)
			     -gammln(en-em+1.0)+em*plog+(en-em)*pclog); 
    } 
    //while (ran1(idum) > t); 
    while (unif_distn() > t);
    // Reject. This happens about 1.5 times per deviate, on average. 
    bnl=em; 
  } 
  if (p != pp)   // Remember to undo the symmetry transformation. 
    bnl=n-bnl; 
  return bnl; 
}

/*
  chooses a random number between 0 and nmax-1 given nmax
  */
int rand_index(int nmax)
{
  // can it go wrong?
  int ri;
  do{
    ri = (int)(nmax*unif_distn());
  } while (ri >= nmax);
  return ri;
}
