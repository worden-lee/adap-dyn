#ifndef MATRIX
#define MATRIX
#include <stdlib.h>
#include "Vector.h"
#include <vector>
#ifdef VNL
#include "vnl/vnl_matrix.h"
#include "vnl/algo/vnl_matrix_inverse.h"
#include "vnl/vnl_diag_matrix.h"
#endif //VNL
#include <iostream>
using std::cerr;
using std::vector;

// this class is now dissed.  prefer to use matrix, below.
template<typename T>
class Matrix : public Vector< Vector<T> >
{
public:
  using Vector< Vector<T> >::n;
  using Vector< Vector<T> >::x;

  Matrix<T>() : Vector< Vector<T> >() { }
  Matrix<T>(long nr, long nc) : Vector< Vector<T> >(nr)
  { setNCols(nc); }
  void resize(long nr, long nc)
  {
    Vector< Vector<T> >::resize(nr);
    setNCols(nc);
  }
  void setNCols(long nc)
  {
    for (long i=0; i<n; i++)
      x[i].resize(nc);
  }
  Matrix<T> &set(const T &v)
  { for (int i = 0; i < n; i++)
      x[i].set(v);
    return *this;
  }

  long rows()
  { return Vector< Vector<T> >::size(); }
  long cols()
  { if (Vector< Vector<T> >::size() > 0)
      return x[0].size();
    else
      return 0;
  }

#ifdef VNL // cast different brands of matrix back and forth
  operator vnl_matrix<T>()
  {
    vnl_matrix<T> v(n, (n>0?x[0].size() : 0));
    for (int i = 0; i<n; i++)
      for (int j=0; j<x[i].size(); j++)
	v.put(i,j,x[i][j]);
    return v;
  }
  Matrix<T>(vnl_matrix<T>&v) : Vector< Vector<T> >(v.rows())
  { setNCols(v.cols());
    for (int i = 0; i<n; i++)
      for (int j=0; j<x[i].size(); j++)
	x[i][j] = v(i,j);  
  }
#endif //VNL
};

template<typename T>
class matrix : public vector< vector<T> >
{
public:
  matrix<T>() : vector< vector<T> >() { }
  matrix<T>(long nr, long nc) : vector< vector<T> >(nr)
  { setNCols(nc); }
  void resize(long nr, long nc, const T &_val = T())
  {
    vector< vector<T> >::resize(nr);
    setNCols(nc,_val);
  }
  void setNCols(long nc, const T &_val = T())
  { for (long i=0; i<rows(); i++)
      (*this)[i].resize(nc,_val);
  }
  matrix<T> &set(const T &v)
  { for (int i = 0; i < rows(); i++)
      (*this)[i].assign(cols(),v);
    return *this;
  }

  inline long rows() const
  { return vector< vector<T> >::size(); }
  inline long cols() const
  { if (rows() > 0)
      return this->begin()->size();
    else
      return 0;
  }

#ifdef VNL // cast different brands of matrix back and forth
  operator vnl_matrix<T>()
  {
    vnl_matrix<T> v(rows(), cols());
    for (int i=0; i<rows(); i++)
      for (int j=0; j<cols(); j++)
	v.put(i,j, (*this)[i][j]);
    return v;
  }
  matrix<T>(vnl_matrix<T>&v) : vector< vector<T> >(v.rows())
  { setNCols(v.cols());
    for (int i=0; i<rows(); i++)
      for (int j=0; j<cols(); j++)
	(*this)[i][j] = v(i,j);  
  }
#endif //VNL
};

template<typename T>
class SquareMatrix : public Matrix<T>
{
public:
  SquareMatrix<T>() : Matrix<T>() {}
  SquareMatrix<T>(long nn) : Matrix<T>(nn, nn) {}
  void resize(long nn)
  { Matrix<T>::resize(nn, nn); }

#ifdef VNL
  SquareMatrix<T>(vnl_matrix<T>&v) : Matrix<T>(v) {}
  // this can automatically be cast to SquareMatrix
  static /* SquareMatrix<T> */vnl_matrix<T> fromEigensystem(SquareMatrix<T> &v, Vector<T> &l)
  { // m = v diag(l) v^-1
    // lots of casting takes place here-- I don't care about efficiency
    vnl_matrix<T> s = vnl_matrix<T>(v) * vnl_diag_matrix<T>(l) * vnl_matrix_inverse<T>(v);
    return s;//SquareMatrix<T>(s);
  }
#endif //VNL
  
};

template<typename T>
class MatrixAccess
{
public:
  class MatrixRowAccess
  {
  public:
    T &operator[] (const long &c)
    { return mx.entry(row,c); } 
    MatrixRowAccess(MatrixAccess&m, long r) : mx(m), row(r) {}
  protected:
    MatrixAccess&mx;
    long row;
  };
  virtual ~MatrixAccess() {}
  //MatrixRowAccess<T> operator[] (const long &r)
  MatrixRowAccess operator[] (const long &r)
  { return MatrixRowAccess(*this,r); }
  // subclass?  implement these three
  virtual unsigned int rows(void) const=0;
  virtual unsigned int cols(void) const=0;
  virtual void resize(unsigned int nr, unsigned int nc) = 0;
  virtual T &entry(const long &r, const long &c)=0;
};

template <class T>
class vnl_matrixAccess : public MatrixAccess<T>
{
protected:
  vnl_matrix<T> *m;
  bool ownm;
public:
  vnl_matrixAccess(vnl_matrix<T> *_m) : m(_m), ownm(false)
  {}
  vnl_matrixAccess(unsigned int rr=0, unsigned int cc=0)
    : m(new vnl_matrix<T>(rr,cc,0)), ownm(true)
  {} 
  unsigned int rows(void) const
  { return m->rows(); } 
  unsigned int cols(void) const
  { return m->cols(); }
  T &entry(const long &r, const long &c)
  { return (*m)[r][c]; }
  // this preserves existing matrix entries
  void resize(unsigned int r, unsigned int c)
  {
    if (r >= rows() && c >= cols())
    { vnl_matrix<T> n(*m,vnl_tag_grab());
      m->set_size(r,c);
      m->fill(0);
      m->update(n);
    }
    else if (r < rows() && c < cols())
      *m = m->extract(r,c);
    else // too bad for you
      cerr << "resize: too much trouble\n";
  }
  ~vnl_matrixAccess()
  { if (ownm) delete m; }
};

template <typename T>
class matrixAccess : public MatrixAccess<T>
{
  template<typename I, typename T_>
  friend class IndexedMatrix;
protected:
  matrix<T> *m;
  bool ownm;
public:
  matrixAccess(matrix<T> &_m) : m(&_m), ownm(false) {}
  matrixAccess(unsigned int rr=0, unsigned int cc=0)
    : m(new matrix<T>(rr,cc,0)), ownm(true)
  {} 
  unsigned int rows(void) const
  { return m->rows(); } 
  unsigned int cols(void) const
  { return m->cols(); }
  T &entry(const long &r, const long &c)
  { return (*m)[r][c]; }
  // this preserves existing matrix entries
  void resize(unsigned int r, unsigned int c)
  { m->resize(r,c); }
  ~matrixAccess()
  { if (ownm) delete m; }
};

// IndexedMatrix uses this weird thing
template<typename T>
class MatrixAccessAccess : public MatrixAccess<T>
{
public:
  MatrixAccess<T> &va;

  MatrixAccessAccess(MatrixAccess<T>&inner) : va(inner) {}
  virtual ~MatrixAccessAccess() {}
  virtual typename MatrixAccess<T>::MatrixRowAccess operator[] (int i)
  { return va[i]; }
  virtual const typename MatrixAccess<T>::MatrixRowAccess operator[] (int i) const
  { return va[i]; }
  virtual unsigned int rows(void) const
  { return va.rows(); }
  virtual unsigned int cols(void) const
  { return va.cols(); }
  virtual void resize(unsigned int nr, unsigned int nc)
  { va.resize(nr,nc); }
  T &entry(const long &r, const long &c)
  { return va[r][c]; }
};

#include "Indexing.h"

template <typename T, typename I=index_type>
class IndexedMatrix : public MatrixAccessAccess<T>
{
public:
  using MatrixAccessAccess<T>::va;
  _Indexing<I> &indexing;
private:
  static MatrixAccess<T> &createV(void)
  { matrix<T>*vv = new matrix<T>;
    matrixAccess<T>*vx = new matrixAccess<T>(*vv);
    return *vx;
  }
  bool ownv;
public:
  class IndexedRowAccess : public MatrixAccess<T>::MatrixRowAccess
  {
  public:
    IndexedRowAccess(MatrixAccess<T>&m, long r)
      : MatrixAccess<T>::MatrixRowAccess(m,r) {}
    using MatrixAccess<T>::MatrixRowAccess::row;
    using MatrixAccess<T>::MatrixRowAccess::mx;
    T &operator[] (_Index<I> c)
    { return
	mx.entry(row, static_cast<IndexedMatrix&>(mx).indexing.index(c));
    }
  };

  IndexedMatrix(MatrixAccess<T> &x, _Indexing<I> &i) :
    MatrixAccessAccess<T>(x), indexing(i), ownv(false) {}
  IndexedMatrix(_Indexing<I> &i) :
    MatrixAccessAccess<T>(createV()), indexing(i), ownv(true) {} 

  virtual ~IndexedMatrix()
  { if (ownv)
    { matrixAccess<T>*vx = static_cast<matrixAccess<T>*>(&va);
      delete &(vx->m);
      delete vx;
    }
  }
//   unsigned int rows()
//   { return va->rows(); } 
//   unsigned int cols()    
//   { return va->cols(); } 
  virtual IndexedRowAccess operator[](_Index<I> r)
  { return IndexedRowAccess(*this, indexing.index(r)); }
  // this one should be IndexedRowAccess but the g++ is giving me troubles
  virtual typename MatrixAccess<T>::MatrixRowAccess operator[] (I i)
  { return IndexedRowAccess(*this, i); }
  // const versions??
};

#endif //MATRIX
