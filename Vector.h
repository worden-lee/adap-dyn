/* -*- C++ -*- 
 * This file includes Vector (and some VectorAccess classes, soon?)
*/
#ifndef VECTOR
#define VECTOR
#include <stdlib.h>
#include <vector>
//#include "Integrator.h"
#ifdef VNL
#include "vnl/vnl_vector.h"
#endif //VNL

template<typename T>
class Vector
{
protected:
  T *x;
  long n;
public:
  Vector<T>(): x(NULL), n(0) {}
  Vector<T>(long nv) : x(NULL), n(0) { resize(nv); }
  Vector<T>(const Vector<T>& v): x(NULL), n(0)
    { *this = v; }
  Vector<T> &operator=(const Vector<T> &y)
  { 
    if (this != &y)
    {
      if (n != y.n)
      { 
	if (n > 0) delete[] x;
	n=y.n;
	x = new T[n];
      }
      for (int i=0; i < n; i++) x[i] = y.x[i];
    }
    return *this;
  }
  ~Vector<T>() { if (n > 0) delete[] x; }
   
  bool operator==(Vector<T> &other)
  { if (n != other.n)
      return false;
    for (int i=0; i < n; i++)
      if (x[i] != other.x[i])
	return false;
    return true;
  }
  bool operator!=(Vector<T> &other)
  { return !(*this == other); }

  long size(void) const { return n; }
  void resize(long nv) // always reallocs, does not initialize new values
  { //cout << " resize " << nv << endl;
    if ( n != nv ) 
    { T *y=x;
      x = new T[nv];
      for (int i=0; i < n; i++) x[i] = y[i];
      delete[] y;
      n = nv;
    }
  }
  Vector<T> &set(const T &v)
  { for (int i = 0; i < n; i++)
      x[i] = v;
    return *this;
  }

  //T &operator[] (const VariableIndex &i) const { return x[i.index()]; }
  T &operator[] (const long &i) const { return x[i]; }

  //VariableIndex firstIndex(void) const { return 0; }
  //VariableIndex finalIndex(void) const { return n; } // this is going too far

#ifdef VNL
  operator vnl_vector<T>()
  {
    return vnl_vector<T>(x, n);
  }
  Vector<T>(const vnl_vector<T>& vn): x(NULL), n(0)
    {
      resize(vn.size());
      for (int i=0; i < n; i++) x[i] = vn[i];
    }
#endif //VNL
};

#endif //VECTOR
