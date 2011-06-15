#include "genrand2.h"
#ifndef RANDOM_H
#define RANDOM_H

/*
  uniform distribution between 0 and 1
 */
inline double unif_distn(void) {return genrand2d();}
#define uniform_distribution unif_distn

/* 
   exponential distribution, specified mean
 */
double exponl_distn(double mean);
#define exponential_distribution exponl_distn

/*
  poisson deviate, specified mean
*/
double poidev(double mean);

/*
  normal distribution, specified mean + std. dev.
*/
double norml_distn(double mean, double stdev);
#define normal_distribution norml_distn

/*
  binomial distribution - # of successes from n trials with prob p of success
*/
double bnldev(double p, int n);

/*
  chooses a random number between 0 and nmax-1
*/
unsigned rand_index(unsigned nmax); // not including nmax

#define generate_poisson_number(M) ((int)poidev(M))
#define uniform_int(a,b)  ((a) + rand_index((b)-(a))) // not including b
#define uniform_double(a,b) ((a) + ((b) - (a))*unif_distn()) 

#endif //RANDOM_H
