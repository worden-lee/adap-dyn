/* -*- C++ -*- */
#ifndef ITERATOR_H
#define ITERATOR_H
#include "Integrator.h"
#include "Indexing.h"

class Iterator : public Integrator
{
public:
  Iterator();
  
  void integrateNonstop(double t1);
  // called by discrete implementation of integrate
  virtual void iterate_to(double t1, double eps, 
		  double extinctionthreshold/*, double dtsav*/);
  virtual double invasionFunction(Index&);
  void integratePartially(double);
  VectorAccess<double> &state(void);
  void setState(double t, VectorAccess<double> *s);
private:
  vector<double> statevector;
  stl_vectorAccess<double> stateaccess;
};

#endif //ITERATOR_H
