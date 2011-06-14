/* -*- C++ -*- */
#ifndef NRINTEGRATOR_H
#define NRINTEGRATOR_H
// #include <stdio.h>
// #include <iostream>
// #include <stdlib.h>
// #include <string.h>
#if VNL
#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#endif
#include "Integrator.h"

class NRIntegrator: public Integrator
{
public:
  NRIntegrator() : svx(statevector) {}
  virtual ~NRIntegrator() {}

  void setState(double t, VectorAccess<double> *v);
  VectorAccess<double> &state(void);

  virtual void integratePartially(double t1);
  virtual void integrateNonstop(double t1);

protected:
  vector<double> statevector;
  stl_dvectorAccess svx;

  void ierror(const char *);
  void odeint(double t2, double eps, 
	      double h1, double hmin, 
	      long *nok, long *nbad/*,
	      double dxsav*/);
  void rkqs(vector<double> &x, vector<double> &dxdt,
	    double htry, double eps,
	    vector<double> &xscal, double *hdid, double *hnext);
  void rkck(vector<double> &x, vector<double> &dxdt, double h,
	    vector<double> &xout, vector<double> &xerr);
};

#endif //NRINTEGRATOR_H
