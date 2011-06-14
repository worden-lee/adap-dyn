/* -*- C++ -*- */
#ifndef INTEGRATOR_H
#define INTEGRATOR_H
//#include "number_classes.h"
#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::vector;

#include "Indexing.h"
#if VNL
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#endif
class Community;
class Site;

// cute little function to duplicate an array,
// make sure it's suitable to go to free() or realloc() later
void *memdup(void *ptr, size_t sz);

class Integrator
{
protected:
  Site *site;
  Indexing indexing;
  double t;
  int/*static*/ evoTime; // need one per integrator
  double noChangeSince;
  double lastEventWas;
  double speciationThreshold;
  double nextImmigration;
  bool timeToSpeciate;
  bool exploded;  
  bool equil;
  bool tiredOfWaiting;
  bool atESS;
public:
  Integrator();
  // copy constructor makes integrators with shared indexing
  // does not copy site or state
  Integrator(const Integrator &in);

  virtual ~Integrator() {};
  void siteIs(Site *);  

  int evolutionaryTime() {
    return evoTime;
  }
  void increaseEvolutionaryTime();

  // set the state, time before beginning integration
  virtual void setState(double t, VectorAccess<double> *v) = 0;

  // try to integrate to time t1, stop if an exceptional condition arises
  // such as equilibrium, extinction, etc.
  virtual void integratePartially(double t1) = 0;
  // do extinctions, speciations, keep integrating until t1 if possible
  virtual void integrateNonstop(double t1) = 0;
  
  virtual void integrate(double t1)
  { integrateNonstop(t1);
  } 

  // should be valid even during integration
  double time(void) { return t; }

  virtual VectorAccess<double> &state(void) = 0;
  double &state(Index i)
  { return state(indexing.index(i)); }
  virtual double &state(int i) //variable indexing
  {
    VectorAccess<double> &sa = state();
    return sa[i];
  }

  void set(Index i, double xi)
  { set(indexing.index(i),xi); }    
  virtual void set(int i, double xi) //variable indexing
  { VectorAccess<double> &sa = state();
    sa[i] = xi;
  } 
  
  Indexing &variableIndexing()
  { return indexing; }

  // notification of community changes
  virtual void newVariable(Index ni);
  // these 2 don't call newVariable, you should call it first
  virtual void speciation(Index parent_index,
			  Index daughter_index);
  // this notifies for either introduction in case of an uncorrelated
  //  assembly model, or diffusion of a species from a neighboring site
  //  in case of a spatial model
  virtual void immigration(Index new_index);
  //  void newResource(long new_j);
  virtual void extinction(Index the_deceased)
  { set(the_deceased,0); }

  virtual void reachedESS(void)
  { atESS = true; }
  bool currentlyAtESS(void)
  { return atESS; }

  // might want to override for efficiency or whatever
  virtual int nextIndex(void)
  { VectorAccess<double> &sa = state();
    int nx = sa.size();
    sa.resize(nx+1);
    return nx;
  }

  Index newIndex(void)
  {
    int ni = nextIndex();
    indexing.registerIndex(ni);
    return Index(ni,&indexing);
  }
  
  // reveal the state of the integrator
  virtual bool hadExplosion(void) { return exploded; }
  virtual bool atEquilibrium(void) { return equil; } 
  virtual bool readyToSpeciate(void) { return timeToSpeciate; }
  
  virtual double invasionFunction(Index&);
#if VNL
  virtual vnl_matrix<double> jacobian(const VectorAccess<double>*x);
#endif
  
protected:
  virtual void resetSpeciationThreshold(void);
  virtual void resetNextImmigration(void);

  friend class CVAccess;
};

#endif //INTEGRATOR_H
