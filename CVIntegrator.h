/* -*- C++ -*- */
#ifndef CVINTEGRATOR_H
#define CVINTEGRATOR_H
#include "Integrator.h"
#include "Matrix.h"

  // This project assumes that SUNDIALS is installed with
  //  realtype == double!  It's true currently, but beware.
#include "sundials/sundials_types.h"
#include "cvode/cvode.h"
#include "nvector/nvector_serial.h"
#include "cvode/cvode_dense.h"

#include <errno.h>
using namespace std;

// this relies on hoping that N_VNew_Serial uses malloc
static N_Vector L_N_VResize(N_Vector v, long int N)
{
  if (N < 0) N = 0;
  if (!v)
    v = N_VNew_Serial(N);
  else
  {
    if (!NV_OWN_DATA_S(v))
      cerr << "Warning in L_N_VResize: doesn't own data!" << endl;
    NV_DATA_S(v) =
      (typeof(NV_DATA_S(v))) realloc(NV_DATA_S(v), N * sizeof(NV_Ith_S(v,0)));
    if (errno == ENOMEM) {
      cerr << "Realloc failed in L_V_Resize()\n";
    }
    if (NV_DATA_S(v) == NULL) {
      free(v);
      return 0;
    }
    NV_LENGTH_S(v) = N;
  }
  return v;
}

// if cv_real ever isn't double, bad luck will ensue
class N_VectorAccess : public VectorAccess<realtype>
{
private:
  N_Vector *v;
  int real_length; // how much allocated
public:
  inline N_VectorAccess(N_Vector *nv, int rl=0) :
    v(nv), real_length(rl)
  { if (rl==0 && nv!=0 && (*nv)!=0) real_length = NV_LENGTH_S((*nv)); }
  virtual ~N_VectorAccess() {}
  //  inline double &operator[] (int i)
  double &operator[] (int i)
  {
    if (i < 0)
      cerr << "!!! index " << i << " too small " << endl;
    else if (i >= (int)size())
      cerr << "!!! index " << i << " too big -- size = " << size() << endl;
    return NV_Ith_S((*v),i);
  }
  inline const double &operator[] (int i) const
  {
    if (i < 0)
      cerr << "index " << i << " too small " << endl;
    else if (i >= (int)size())
      cerr << "index " << i << " too big -- size = " << size() << endl;
    return NV_Ith_S((*v),i);
  }
  inline unsigned int size(void) const
  { return *v ? NV_LENGTH_S((*v)) : 0; }
  void resize(unsigned int r)
  { if (r == size())
      return;
    if (!*v)
    { (*v) = N_VNew_Serial(r);
      real_length = 0;
    }
    else if ((int)r<real_length && real_length < 2*(int)r) // optimize
      NV_LENGTH_S((*v)) = r;
    else
    { L_N_VResize(*v,r);
      real_length = r;
    }
  }
};

class CVIntegrator;
class CVAccess
{
public:
  CVAccess(CVIntegrator&ci) : in(ci) {}
  virtual ~CVAccess() {}
  
  virtual void calcDerivatives(const VectorAccess<double> *x,
			       VectorAccess<double> *dxdt);
  virtual bool doJacobian(void)
  { return false; }  
//   virtual void jacobian(VectorAccess<double> *x, RhsFn f,
// 			MatrixAccess<double> *J) {}
  virtual int throwIfNeeded(void);
  virtual int check(const VectorAccess<double> *x, realtype t);
  //virtual void considerExtinction(VectorAccess<double> &x, const Index &i);
private:
  CVIntegrator &in;
};

class CVIntegrator : public Integrator
{
public:
  CVIntegrator() :
    ax(*this), nv(0), nvx(&nv), cvode_mem(0), cvode_mem_N(0), reltol(0),
    integrating(0), discovery(false), extnctn(false),
    restartCVODE(false)
  {
    //    integrateFlag = SUCCESS;
  }

  // using this copy constructor suffices for a deepCopy():
  //  newsite->integrator = new CVIntegrator(oldsite->integrator);
  CVIntegrator(CVIntegrator &other) :
    Integrator(other), ax(*this), nv(0), nvx(&nv),
    cvode_mem(0), cvode_mem_N(0), reltol(0),
    integrating(0), discovery(false), extnctn(false),
    restartCVODE(false)
  { //setState(other.time(), &other.state());
    //    integrateFlag = SUCCESS;
  }
  
  ~CVIntegrator()
  {
    if (cvode_mem!=0)
      CVodeFree(&cvode_mem);
    nvx.resize(0);
    //    integrateFlag = SUCCESS;
  }
  
  void integrateNonstop(double t1);
  // this (and check) should have no side effects to outside objects
  void integratePartially(double t1);
  
  virtual void setState(double t, VectorAccess<double> *v);
  virtual VectorAccess<double> &state(void);

  void extinction(Index the_deceased);

  int check(const VectorAccess<double> *x, realtype t);

  virtual void newVariable(Index ni)
  {
    Integrator::newVariable(ni);
    discovery = true;
  }

  // integrateFlag is the error code from integratePartially
  int integrateFlag;

  bool hadExtinction(void) { return extnctn; }
  
protected:
  CVAccess ax;
  N_Vector nv;
  N_VectorAccess nvx;
  void *cvode_mem;
  unsigned int cvode_mem_N;
  double reltol, abstol;
  int integrating;  // nonzero while within a call to integrate()
  bool discovery;   // true if a new variable has been introduced
		    // and the cvode_mem is out of date
  bool extnctn;     // true if something needs to be eliminated
  bool restartCVODE; // true if cvode_mem is out of date
  friend class CVAccess;
};

#endif //CVINTEGRATOR_H
