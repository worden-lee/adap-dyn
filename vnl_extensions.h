/*   -*- C++ -*-
 * extra stuff to go with the vnl matrix/vector library
 *  by lee w. 7/04
 */
#ifndef VNL_EXTENSIONS_H
#define VNL_EXTENSIONS_H

#include "vnl/vnl_vector.h"
#include "vnl/vnl_matrix.h"
#include "vnl/vnl_diag_matrix.h"

template <class T>
vnl_vector<T> tensor_product(vnl_vector<T> v, vnl_vector<T> w)
{
  vnl_vector<T> vxw(v.size()*w.size());
  const unsigned int n = w.size();
  for (unsigned int i = 0; i < v.size(); i++)
    for (unsigned int j = 0; j < n; j++)
      vxw(n*i + j) = v[i]*w[j];
  return vxw;
}

template <class T>
vnl_matrix<T> tensor_product(vnl_matrix<T> v, vnl_matrix<T> w)
{
  vnl_matrix<T> vxw(v.rows()*w.rows(),v.cols()*w.cols());
  const unsigned int m = w.rows();
  const unsigned int n = w.cols();
  for (unsigned int i = 0; i < v.rows(); i++)
    for (unsigned int j = 0; j < v.cols(); j++)
      for (unsigned int k = 0; k < m; k++)
	for (unsigned int l = 0; l < n; l++)
	  vxw(m*i + k, n*j + l) = v(i,j)*w(k,l);
  return vxw;
}

// tensor_product not defined as yet for matrix/vector and v.v.

// generalization of apply method, which is only given
// for the case T0 == T1
template <class T0, class T1>
vnl_vector<T1> apply(T1 (*f)(const T0 &), const vnl_vector<T0>&vec)
{ unsigned n = vec.size();
  vnl_vector<T1> dest(n);
  for (unsigned i = 0; i < n; ++i)
    dest[i] = f(vec[i]);
  return dest;
}

// specialize a tiny bit so apply(vcl_real...) works okay
template <class T0, class T1>
vnl_vector<T1> apply(const T1& (*f)(const T0 &), const vnl_vector<T0>&vec)
{ unsigned n = vec.size();
  vnl_vector<T1> dest(n);
  for (unsigned i = 0; i < n; ++i)
    dest[i] = f(vec[i]);
  return dest;
}

template <class T0, class T1>
vnl_matrix<T1> apply(T1 (*f)(const T0 &), const vnl_matrix<T0>&M)
{ unsigned r = M.rows();
  unsigned c = M.cols();
  vnl_matrix<T1> dest(r,c);
  for (unsigned i = 0; i < r; ++i)
    for (unsigned j = 0; j < c; ++j)
      dest[i][j] = f(M[i][j]);
  return dest;
}

// template <class T0, class T1>
// T1 cast_function(T0&x)
// {
//   return static_cast<T1>(x);
// }

// probably outmoded: in vxl-1.1.0 better to use
//  vnl_matrix(len,len,vnl_matrix_identity)
template <class T>
class identity_matrix : public vnl_diag_matrix<T> 
{
public:
  identity_matrix() {}
  identity_matrix(int len) : vnl_diag_matrix<T>(len,1) {}
  // notice .as_matrix() function may be very useful
};

// see vnl_svd : this is inspired by that.
// it decomposes your matrix M into form V D W,
// with D a diagonal matrix of eigenvalues,
// V a matrix whose columns are right eigenvectors,
// W (the inverse of V) with left eigenvectors as rows.
#include <vnl_real_eigensystem.h>
#include <vnl_matrix_inverse.h>
#include <vnl_complexify.h>
#include <vnl_real.h>
class vnl_real_diagonalization
{
public:
  vnl_real_diagonalization(vnl_matrix<double> const &M)
    : eig(M), V_(fixPhase(eig.V)), D_(eig.D), W_(0,0), W__(V_)
  { // sort eigenvalues?
  }

  vnl_matrix< vcl_complex<double> > &V() 
  { return V_; }
  vnl_diag_matrix< vcl_complex<double> > &D()
  { return D_; }
  vnl_matrix< vcl_complex<double> > &W()
  { ensureW();
    return W_;
  }

  //: Solve the matrix equation M X = B, returning X
  vnl_matrix< vcl_complex<double> >
     solve (vnl_matrix< vcl_complex<double> > const& B) 
  { ensureDinverse();
    return V_*Dinverse_*W__*B;
  }
  vnl_matrix< double >
     solve (vnl_matrix< double > const& B) 
  { return vnl_real(solve(vnl_complexify(B)));
  }
  
  //: Solve the matrix-vector system M x = y, returning x.
  vnl_vector< vcl_complex<double> >
     solve (vnl_vector< vcl_complex<double> > const& y) 
  { ensureDinverse();
    return V_*Dinverse_*W__*y;
  }
  vnl_vector< double >
     solve (vnl_vector< double > const& y) 
  { return vnl_real(solve(vnl_complexify(y)));
  }

private:
  vnl_real_eigensystem eig;
  vnl_matrix< vcl_complex<double> > V_;
  vnl_diag_matrix< vcl_complex<double> > D_;
  vnl_diag_matrix< vcl_complex<double> > Dinverse_;
  vnl_matrix< vcl_complex<double> > W_;
  vnl_matrix_inverse< vcl_complex<double> > W__;

  // line up phase of top entry of eigenvectors, in place
  vnl_matrix< vcl_complex<double> >
    &fixPhase(vnl_matrix< vcl_complex<double> > &V)
  {
    vnl_diag_matrix< vcl_complex<double> > sgn( V.cols() );
    for (unsigned int i = 0; i < V.cols(); ++i)
      sgn[i] = 1.0/(V[0][i]/abs(V[0][i]));
    V = V*sgn;
    V.normalize_columns();
    return V;
  }
  void ensureW(void)
  { if (W_.rows()==0) W_ = W__; // force computation of inverse
  }
  void ensureDinverse(void)
  { if (Dinverse_.size()!=0) return;
    Dinverse_ = vnl_diag_matrix< vcl_complex<double> >(D_.size(),0);
    for (unsigned int i = 0; i < D_.size(); ++i)
      if (D_[i] != 0.0)
	Dinverse_[i] = 1.0/D_[i];
  }
  
    // Disallow assignment.
  vnl_real_diagonalization& operator=(vnl_real_diagonalization const &)
  { return *this; }
};
    
#include "Indexing.h"
class vnl_dvectorAccess : public VectorAccess<double>
{
private:
  vnl_vector<double> &v;
public:
  vnl_dvectorAccess(vnl_vector<double>&vv) : v(vv) {}
  virtual ~vnl_dvectorAccess() {}
  virtual double &operator[] (int i)
  { return v[i]; }
  virtual const double &operator[] (int i) const
  { return v[i]; }
  virtual unsigned int size(void) const
  { return v.size(); }
  virtual void resize(unsigned int ns)
  { cerr << "error: vnl_dvectorAccess doesn't resize" << endl;
    return; 
  }
};

#if 0
// can't use this class any more because in vxl 1.5 too much is private
// after this code block is an alternate approach
#include <vnl_real_npolynomial.h>
// fixing some holes in the vnl_real_npolynomial class
class l_vnl_real_npolynomial : public vnl_real_npolynomial
{
public:
  l_vnl_real_npolynomial(const vnl_vector<double>&c,
			 const vnl_matrix<unsigned int>& p)
    : vnl_real_npolynomial(c, p)
  { simplify(); }  // unfortunately these call the parent
		   // class version as well
  l_vnl_real_npolynomial(const vnl_real_npolynomial &v)
    : vnl_real_npolynomial(v)
  { simplify(); }
  l_vnl_real_npolynomial() {}
  l_vnl_real_npolynomial &operator=(const vnl_real_npolynomial &v)
  {
    vnl_real_npolynomial::operator=(v);
    simplify();
    return *this;
  }

  void print()
  {
    cout << *this << "\n";
  }
  
  //: Combine terms with identical exponents (i.e., identical rows in polyn_).
  // Remove terms with zero coefficient.
  void simplify()
  {
    print();
    for (unsigned int row1=0; row1<nterms_; ++row1)
      for (unsigned int row2=row1+1; row2<nterms_; ++row2) {
	unsigned int col=0;
	while (col<nvar_ && polyn_(row1,col) == polyn_(row2,col)) ++col;
	if (col < nvar_) continue; // not all exponents are identical
	coeffs_(row1) += coeffs_(row2); coeffs_(row2) = 0;
      }
    print();
    while (nterms_>0 && coeffs_(nterms_-1)==0) --nterms_;
    for (unsigned int row=0; row<nterms_; ++row)
      if (coeffs_(row) == 0) {
	--nterms_; // decrement nterms, and move last element to vacant place:
	coeffs_(row) = coeffs_(nterms_);
	coeffs_(nterms_) = 0; // not really necessary; to keep coeffs_ consistent
	for (unsigned int i=0; i<nvar_; ++i)
	  polyn_(row,i) = polyn_(nterms_,i);
      }
    print();
  }
  
  // bug fixes: handles -1 coefficient correctly;
  // doesn't append a newline.
  // also places spaces between factors involving powers
  friend vcl_ostream& operator<<(vcl_ostream& os,
				 l_vnl_real_npolynomial const& P)
  {
    if (P.nvar_ <= 3)
      for (unsigned int i=0; i<P.nterms_; ++i)
      {
	if (i>0)
	{
	  os << ' ';
	  if (P.coeffs_(i) >= 0)
	    os << "+ ";
	}
	if (P.coeffs_(i) < 0)
	  os << "- ";
	if (fabs(P.coeffs_(i)) != 1)
	  os << fabs(P.coeffs_(i)) << ' ';
	unsigned int totaldeg = 0;
	if (P.nvar_ > 0 && P.polyn_(i,0) > 0)
	{ os << 'X'; totaldeg += P.polyn_(i,0); }
	if (P.nvar_ > 0 && P.polyn_(i,0) > 1)
	  os << '^' << P.polyn_(i,0) << ' ';
	if (P.nvar_ > 1 && P.polyn_(i,1) > 0)
	{ os << 'Y'; totaldeg += P.polyn_(i,1); }
	if (P.nvar_ > 1 && P.polyn_(i,1) > 1)
	  os << '^' << P.polyn_(i,1) << ' ';
	if (P.nvar_ > 2 && P.polyn_(i,2) > 0)
	{ os << 'Z'; totaldeg += P.polyn_(i,2); }
	if (P.nvar_ > 2 && P.polyn_(i,2) > 1)
	  os << '^' << P.polyn_(i,2) << ' ';
	if (totaldeg == 0 && vcl_fabs(P.coeffs_(i)) == 1)
	  os << fabs(P.coeffs_(i));
      }
    else
      for (unsigned int i=0; i<P.nterms_; ++i)
      {
	if (i>0)
	{
	  os << ' ';
	  if (P.coeffs_(i) >= 0)
	    os << "+ ";
	}
	if (P.coeffs_(i) < 0)
	  os << "- ";
	if (fabs(P.coeffs_(i)) != 1)
	  os << fabs(P.coeffs_(i)) << ' ';
	unsigned int totaldeg = 0;
	for (unsigned int j=0; j<P.nvar_; ++j) {
	  if (P.polyn_(i,j) > 0)  os << 'X' << j;
	  if (P.polyn_(i,j) > 1)  os << '^' << P.polyn_(i,j) << ' ';
	  totaldeg += P.polyn_(i,j);
	}
	if (totaldeg == 0 && vcl_fabs(P.coeffs_(i)) == 1)
	  os << fabs(P.coeffs_(i));
      }
    //os << endl;
    return os;
  }

};
#endif //0
#include <vnl_real_npolynomial.h>
// vnl_real_npolynomial::simplify() is buggy; use this after creating
//one
// Combine terms with identical exponents (i.e., identical rows
//in polyn_).  Remove terms with zero coefficient.
inline void simplify(vnl_real_npolynomial &poly)
{
  vnl_vector<double> coeffs_ = poly.coefficients();
  vnl_matrix<unsigned int> polyn_ = poly.polyn();
  unsigned int             nvar_ = polyn_.cols();
  unsigned int             nterms_ = polyn_.rows();

  for (unsigned int row1=0; row1<nterms_; ++row1)
    for (unsigned int row2=row1+1; row2<nterms_; ++row2) {
      unsigned int col=0;
      while (col<nvar_ && polyn_(row1,col) == polyn_(row2,col)) ++col;
      if (col < nvar_) continue; // not all exponents are identical
      coeffs_(row1) += coeffs_(row2); coeffs_(row2) = 0;
    }
  //print();
  while (nterms_>0 && coeffs_(nterms_-1)==0) --nterms_;
  for (unsigned int row=0; row<nterms_; ++row)
    if (coeffs_(row) == 0) {
      --nterms_; // decrement nterms, and move last element to
		 // vacant place:
      coeffs_(row) = coeffs_(nterms_);
      coeffs_(nterms_) = 0; // not really necessary; to keep
			    // coeffs_ consistent
      for (unsigned int i=0; i<nvar_; ++i)
	polyn_(row,i) = polyn_(nterms_,i);
    }
  if (nterms_ < polyn_.rows())
    poly.set(coeffs_.extract(nterms_),polyn_.extract(nterms_,nvar_));
  else
    poly.set(coeffs_,polyn_);
}

// use stream << asString(P) instead of stream << P
// it looks better
inline string asString(vnl_real_npolynomial &P)
{
  ostringstream os;
  vnl_vector<double> coeffs_ = P.coefficients();
  vnl_matrix<unsigned int> polyn_ = P.polyn();
  unsigned int             nvar_ = polyn_.cols();
  unsigned int             nterms_ = polyn_.rows();

  if (nvar_ <= 3)
    for (unsigned int i=0; i<nterms_; ++i)
    {
      if (i>0)
      {
	os << ' ';
	if (coeffs_(i) >= 0)
	  os << "+ ";
      }
      if (coeffs_(i) < 0)
	os << "- ";
      if (fabs(coeffs_(i)) != 1)
	os << fabs(coeffs_(i)) << ' ';
      unsigned int totaldeg = 0;
      if (nvar_ > 0 && polyn_(i,0) > 0)
      { os << 'X'; totaldeg += polyn_(i,0); }
      if (nvar_ > 0 && polyn_(i,0) > 1)
	os << '^' << polyn_(i,0) << ' ';
      if (nvar_ > 1 && polyn_(i,1) > 0)
      { os << 'Y'; totaldeg += polyn_(i,1); }
      if (nvar_ > 1 && polyn_(i,1) > 1)
	os << '^' << polyn_(i,1) << ' ';
      if (nvar_ > 2 && polyn_(i,2) > 0)
      { os << 'Z'; totaldeg += polyn_(i,2); }
      if (nvar_ > 2 && polyn_(i,2) > 1)
	os << '^' << polyn_(i,2) << ' ';
      if (totaldeg == 0 && vcl_fabs(coeffs_(i)) == 1)
	os << fabs(coeffs_(i));
    }
  else
    for (unsigned int i=0; i<nterms_; ++i)
    {
      if (i>0)
      {
	os << ' ';
	if (coeffs_(i) >= 0)
	  os << "+ ";
      }
      if (coeffs_(i) < 0)
	os << "- ";
      if (fabs(coeffs_(i)) != 1)
	os << fabs(coeffs_(i)) << ' ';
      unsigned int totaldeg = 0;
      for (unsigned int j=0; j<nvar_; ++j) {
	if (polyn_(i,j) > 0)  os << 'X' << j;
	if (polyn_(i,j) > 1)  os << '^' << polyn_(i,j) << ' ';
	totaldeg += polyn_(i,j);
      }
      if (totaldeg == 0 && vcl_fabs(coeffs_(i)) == 1)
	os << fabs(coeffs_(i));
    }
  //os << endl;
  return os.str();
}
  
#include <vnl_real_polynomial.h>
inline void simplify(vnl_real_polynomial &p)
{
  vnl_vector<double> &coeffs_ = p.coefficients();
  // vnl_rpoly_roots insists on no leading zeros
  unsigned first = 0;
  const double threshold = 1e-8;
  while (first<coeffs_.size() && coeffs_[first]<threshold)
    ++first;
  if (first > 0)
    coeffs_ = coeffs_.extract(coeffs_.size()-first,first);

  p.set_coefficients(coeffs_);
}

inline string asString(vnl_real_polynomial &p)
{
  ostringstream os;
  int d = p.degree();
  vnl_vector<double> &coeffs_ = p.coefficients();
  int i = 0;
  while (i <= d && coeffs_[i] == 0) ++i;
  if (i > d) { os << "0 "; return os.str(); }
  //  bool b = (coeffs_[i+1] > 0); // to avoid '+' in front of equation
  bool b = (coeffs_[i] > 0); // to avoid '+' in front of equation

  for (; i <= d; ++i) {
    if (coeffs_[i] == 0) continue;
    if (coeffs_[i] < 0)        os << "- ";
    else if (!b)               os << "+ ";
    b = false;
    if (d == 0 || fabs(coeffs_[i]) != 1)
                               os << fabs(coeffs_[i]) << ' ';

    if (d-i > 1)               os << "X^" << d-i << ' ';
    else if (d-i == 1)         os << "X ";
  }
  return os.str();
}

#endif //VNL_EXTENSIONS_H
