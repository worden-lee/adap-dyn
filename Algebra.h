#ifndef ALGEBRA_H
#define ALGEBRA_H
/* -*- C++ -*- */

// this is fun

#include "util.h"
#include <map>
#include <list>
#include <iostream>
#include <sstream>
#include "vnl_math.h"
#include "vnl_real_polynomial.h"
#include "vnl_real_npolynomial.h"
#include "vnl_numeric_traits.h"
#include "vnl_complex_traits.h"
using namespace std;

inline void print(vnl_real_npolynomial &poly)
{
  cout << poly.coefficients() << endl
       << poly.polyn() << endl;
}

static const double tolerance = 1e-4;
static bool aboutEqual(double x, double y)
{ return abs(y-x) < tolerance;
}

// this is used in making polynomials
class VariableAssignment : public unary_function<string,unsigned int>
{
public:
  virtual ~VariableAssignment() {}
  // use operator() to convert a variable name into column index
  virtual unsigned int operator()(const string &) = 0;
  virtual unsigned int size() = 0; // number of columns
};

class AlgebraExpression;
typedef map<string,AlgebraExpression> VariableSubstitution;

// various kinds of AlgebraElements -
//  variables, scalars, operators, functions
class AlgebraElement
{
public:
  virtual ~AlgebraElement() {}
  virtual AlgebraElement *normalize()
  { return this; }
  virtual AlgebraElement *deepcopy() const = 0;

  virtual AlgebraElement *substitute(VariableSubstitution&)
  { return deepcopy(); }
  
  // for debugging
  virtual string prettyPrint() const { return asString(); }
  // for <<
  // this writes the expression in algebra notation
  // when the 'code'  flag is set, it writes c, i.e. with * for multiply
  virtual string asString(bool code=false) const = 0;
  // for compatibility with some vnl math that shouldn't be
  // seriously used
  virtual double asDouble() const { return 0; }
  // making a vnl object for polynomial evaluation
  virtual vnl_real_polynomial as_polynomial() = 0;
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
    = 0;
  bool operator==(const AlgebraElement &that) const
  { return asString() == that.asString();
  }
  friend ostream&operator<<(ostream&,AlgebraElement&);
  // for sorting arguments to an operator
  virtual inline bool operator<(const AlgebraElement&b) const
  { return sortKey() < b.sortKey(); }
  virtual string sortKey() const { return ""; }
  virtual string sortKey_recursive() const { return ""; }
};

class ScalarElement : public AlgebraElement
{
public:
  double _value;
  ScalarElement(double n) : _value(n) {}
  double value() const
  { return _value; }
  string asString(bool code=false) const
  { std::ostringstream out;
    out << _value;
    return out.str();
  }
  string prettyPrint() const
  { return asString(); }
  double asDouble() const
  { return value(); }
  AlgebraElement *normalize();
  virtual AlgebraElement *deepcopy() const
  { return new ScalarElement(_value); }
  virtual ~ScalarElement() {}
  ScalarElement operator+(const ScalarElement &that)
  { return *new ScalarElement(_value + that._value);
  }
  ScalarElement operator*(const ScalarElement &that)
  { return *new ScalarElement(_value * that._value);
  }
  // numbers come before everything else
  //  and don't get compared to each other
  string sortKey() const { return "0 "; }

  virtual vnl_real_polynomial as_polynomial()
  { return vnl_real_polynomial(_value);
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { vnl_vector<double> coeff(1,_value);
    vnl_matrix<unsigned int> exps(1,va.size(),0);
    return vnl_real_npolynomial(coeff,exps);
  }
};

class VariableElement : public AlgebraElement
{ 
public:
  string _name;
  VariableElement(const string &n) : _name(n) {}
  string asString(bool code=false) const
  { return _name; }
  // variables come after numbers, before everything else
  string sortKey_recursive() const
  { return _name; }
  string sortKey() const
  { return "1 "+sortKey_recursive(); }
  virtual AlgebraElement *deepcopy() const
  { return new VariableElement(_name); }

  virtual AlgebraElement *substitute(VariableSubstitution&sub);
  
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { unsigned int index = va(_name);
    vnl_vector<double> coeffs(1,1);
    vnl_matrix<unsigned int> exps(1,va.size());
    exps[0][index] = 1;
    return vnl_real_npolynomial(coeffs,exps);
  }
  // this kind of polynomial assumes 1 unknown, so we don't distinguish
  virtual vnl_real_polynomial as_polynomial()
  { const double _varcoeffs[] = { 1, 0 }; 
    const vnl_real_polynomial varpoly(_varcoeffs,2);
    return varpoly;
  }
  virtual ~VariableElement() {} 
};

class OperatorElement : public AlgebraElement
{
public:
  virtual unsigned precedence() const = 0;
  // operator expressions get sorted after variables
  // ... except for products
  string sortKey() const
  { return "2 "+sortKey_recursive(); }
  virtual ~OperatorElement() {} 
};

class UnaryOperatorElement : public OperatorElement
{
public:
  AlgebraElement *argument;
  UnaryOperatorElement() : argument(0) {}
  void addArgument(const AlgebraElement*arg)
  { if (argument) delete argument;
    argument = arg? arg->deepcopy() : 0;
  } 
  virtual AlgebraElement *substitute(VariableSubstitution&sub)
  { UnaryOperatorElement *n
      = static_cast<UnaryOperatorElement*>(deepcopy());
    n->addArgument(n->argument->substitute(sub));
    return n;
  }
  string sortKey_recursive() const
  { return argument->sortKey_recursive(); } 
  string sortKey() const
  { return argument->sortKey(); } 
  virtual ~UnaryOperatorElement()
  { if (argument) delete argument; } 
};

class NegateElement : public UnaryOperatorElement
{
public:
  NegateElement(const AlgebraElement*arg=0)
  { if (arg) addArgument(arg);
  }
  NegateElement(const AlgebraElement &arg)
  { addArgument(&arg);
  }
  unsigned precedence() const
  { return 1; } 
  string asString(bool code=false) const
  { return "- " + argument->asString(code); }
  string prettyPrint() const
  { return "Negate( "+argument->prettyPrint()+" )"; }
  AlgebraElement *normalize();
  virtual AlgebraElement *deepcopy() const
  { return new NegateElement(argument->deepcopy()); }
  virtual vnl_real_polynomial as_polynomial()
  { return argument->as_polynomial() * -1.0;
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { return argument->as_npolynomial(va) * -1.0; }
  virtual ~NegateElement() {}
};

// for sorting arguments
struct compareArguments
  : public binary_function< const AlgebraElement*,
			    const AlgebraElement*, bool >
{
  bool operator()( const AlgebraElement* a,
		   const AlgebraElement* b ) const;
};

class NaryOperatorElement : public OperatorElement
{
public:
  list<AlgebraElement*> arguments;
  void addArgument(const AlgebraElement*arg)
  { arguments.push_back(arg->deepcopy()); }
  void clearArguments(void)
  { while (arguments.size() > 0)
    { AlgebraElement*f = arguments.front();
      arguments.pop_front();
      delete f;
    }
  }
  virtual string operatorString(bool code=false) const = 0;
  virtual string operatorName() const
  { return operatorString(false); }
  string asString(bool code=false) const
  {
    string as;
    for (list<AlgebraElement*>::const_iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
    {
      if (as.size() > 0)
	as += operatorString(code);
      OperatorElement *oarg = dynamic_cast<OperatorElement*>(*ai);
      // if it's an operator, maybe put it in parens
      if (oarg && oarg->precedence() < precedence())
	as += '(' + oarg->asString(code) + ')';
      else
	as += (*ai)->asString(code);
    }
    return as;
  }
  string prettyPrint() const
  { string pp = operatorName() + "( ";
    for (list<AlgebraElement*>::const_iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
    {
      if (ai != arguments.begin())
	pp += ", ";
      pp += (*ai)->prettyPrint();
    }
    pp += " )";
    return pp;
  }
  // this assumes an associative operator
  AlgebraElement *normalize()
  {
  //cout << "AlgebraElement::normalize(): from " << prettyPrint() << endl;
    for (list<AlgebraElement*>::iterator it = arguments.begin();
	 it != arguments.end(); )
    {
      *it = (*it)->normalize();
      NaryOperatorElement*nop = dynamic_cast<NaryOperatorElement*>(*it);
      // if it's the same kind of operator, splice them together
      if (nop && nop->operatorName() == operatorName())
      {
	list<AlgebraElement*>::iterator nit = it++;
	// take all its arguments
	arguments.splice(it,nop->arguments);
	delete nop;
	arguments.erase(nit);
      }
      else
	++it;
    }
    arguments.sort(compareArguments());
    //cout << " to " << prettyPrint() << endl;
    return this;
  }
  virtual AlgebraElement *substitute(VariableSubstitution&sub)
  { NaryOperatorElement *n
      = static_cast<NaryOperatorElement*>(deepcopy());
    n->clearArguments();
    for (list<AlgebraElement*>::iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
      n->addArgument((*ai)->substitute(sub));
    return n;
  }
  string sortKey_recursive() const
  { for (list<AlgebraElement*>::const_iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
    { string ski = (*ai)->sortKey_recursive();
      if (ski.length() > 0)
	return ski;
    }
    return "";
  }
  virtual ~NaryOperatorElement()
  { clearArguments(); }
};

class AddElement : public NaryOperatorElement
{
public:
  AddElement() {}
  AddElement(const AlgebraElement &t0, const AlgebraElement &t1)
  {
    addArgument(&t0);
    addArgument(&t1);
  }
  unsigned precedence() const
  { return 1; } 
  string operatorString(bool code=false) const
  { return " + "; }
  virtual string operatorName() const
  { return "Add"; }
  // special case because of '+' and '-'
  string asString(bool code=false) const
  {
    string as;
    for (list<AlgebraElement*>::const_iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
    {
      if (as.size() > 0)
      {
	NegateElement *narg = dynamic_cast<NegateElement*>(*ai);
	// if it's a negative, skip the + 
	if (narg)
	  as += ' ';
	else
	  as += operatorString(code);
      }
      // if it's an operator, maybe put it in parens
      OperatorElement *oarg = dynamic_cast<OperatorElement*>(*ai);
      if (oarg && oarg->precedence() < precedence())
	as += '(' + oarg->asString(code) + ')';
      else
	as += (*ai)->asString(code);
    }
    return as;
  }
  AlgebraElement *normalize()
  {
    // this normalizes + sorts arguments.
    // if it ever doesn't return this, we'll have trouble
    NaryOperatorElement::normalize();
    //cout << "AddElement::normalize(): from " << prettyPrint() << endl;
    if (arguments.size() > 0) // should be nonzero!
    {
      // special dealing with numbers
      double number = 0;
      // consolidate numbers into one
      for (list<AlgebraElement*>::iterator it = arguments.begin();
	   it != arguments.end(); )
      {
	ScalarElement *z;
	NegateElement *n;
	if ((z = dynamic_cast<ScalarElement*>(*it)))
	{
	  number += z->_value;
	  arguments.erase(it++);
	}
	else if ((n = dynamic_cast<NegateElement*>(*it)) &&
		 (z = dynamic_cast<ScalarElement*>(n->argument)))
	{
	  number -= z->_value;
	  arguments.erase(it++);
	}
	else
	  ++it;
      }
      // at this point all numbers have been removed from the list
      // if there is nothing else, return the number
      if (arguments.size() == 0)
      {
	delete this;
	return new ScalarElement(number);
      }
      // else if the number is nonzero, add it in
      if (!aboutEqual(0,number))
	arguments.push_front(new ScalarElement(number));

      // if there's only one term, return it.
      if (arguments.size() == 1)
      {
	AlgebraElement *el = arguments.front();
	arguments.clear();
	//cout << " to " << el->prettyPrint() << endl;
	delete this;
	return el;
      }
    }
    //cout << " to " << prettyPrint() << endl;
    return this;
  }
  AlgebraElement *deepcopy() const
  { AddElement *ne = new AddElement;
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      ne->arguments.push_back((*it)->deepcopy());
    return ne;
  }
  virtual vnl_real_polynomial as_polynomial()
  {
    vnl_real_polynomial p(0.0);
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      p = p + (*it)->as_polynomial();
    return p;
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { vnl_real_npolynomial p = ScalarElement(0).as_npolynomial(va);
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
    { vnl_real_npolynomial pp = p + (*it)->as_npolynomial(va);
      p.set(pp.coefficients(),pp.polyn()); // no += operator
    }
    //cout << asString() << endl; print(p);
    return p;
  }
  virtual ~AddElement() {} 
};

class MultiplyElement : public NaryOperatorElement
{
public:
  MultiplyElement() {} 
  MultiplyElement(const AlgebraElement&f0, const AlgebraElement&f1)
  {
    addArgument(&f0);
    addArgument(&f1);
  }
  unsigned precedence() const
  { return 3; } 
  string operatorString(bool code=false) const
  { return code ? " * " : " "; }
  virtual string operatorName() const
  { return "Multiply"; }
  AlgebraElement *normalize();
  string sortKey() const
  { return "1 "+sortKey_recursive(); }
  virtual AlgebraElement *deepcopy() const
  { MultiplyElement *nme = new MultiplyElement;
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      nme->arguments.push_back((*it)->deepcopy());
    return nme;
  }
  virtual vnl_real_polynomial as_polynomial()
  {
    vnl_real_polynomial p(1.0);
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      p = p * (*it)->as_polynomial();
    return p;
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { vnl_real_npolynomial p = ScalarElement(1.0).as_npolynomial(va);
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
    { vnl_real_npolynomial pp = p * (*it)->as_npolynomial(va);
      p.set(pp.coefficients(),pp.polyn()); // no *= operator
    }
    //cout << asString() << endl; print(p);
    return p;
  }
  virtual ~MultiplyElement() {} 
};

class BinaryOperatorElement : public NaryOperatorElement
{
public:
  void setArguments(const AlgebraElement &a0, const AlgebraElement &a1)
  { clearArguments();
    NaryOperatorElement::addArgument(&a0);
    NaryOperatorElement::addArgument(&a1);
  }
  void addArgument(const AlgebraElement*arg)
  { cerr << "Error: don't call BinaryOperatorElement::addArgument() !!"
	 << endl;
  }
  AlgebraElement *&first() 
  {
    return *arguments.begin();
  }
  AlgebraElement *&second() 
  {
    return *++arguments.begin();
  }
  virtual AlgebraElement *substitute(VariableSubstitution&sub)
  { BinaryOperatorElement *n
      = static_cast<BinaryOperatorElement*>(deepcopy());
    n->setArguments(*(first()->substitute(sub)),
		    *(second()->substitute(sub)));
    return n;
  }
  virtual ~BinaryOperatorElement()
  {} 
};

// first / second
class DivideElement : public BinaryOperatorElement
{
public:
  DivideElement() {}
  DivideElement(const AlgebraElement &num, const AlgebraElement &denom)
  { setArguments(num,denom);
  } 
  unsigned precedence() const
  { return 2; } 
  string operatorString(bool code=false) const
  { return " / "; } 
  virtual string operatorName() const
  { return "Divide"; }
  // do not expect too much smartness about simplifying here
  AlgebraElement *normalize()
  { // do not sort the two arguments
    // resolve fractions of fractions
    first() = first()->normalize();
    DivideElement *nde = dynamic_cast<DivideElement*>(first());
    if (nde)
    {
      first() = nde->first();
      second() = new MultiplyElement(*nde->second(),*second());
    }
    second() = second()->normalize();
    DivideElement *dde = dynamic_cast<DivideElement*>(second());
    if (dde)
    {
      second() = dde->first();
      first() = new MultiplyElement(*dde->second(),*first());
      first() = first()->normalize();
    }
    // now cancel any common terms
    list<AlgebraElement*>tops, bottoms;
    MultiplyElement *nme = dynamic_cast<MultiplyElement*>(first());
    if (nme)
      tops = nme->arguments;
    else tops.push_back(first());
    MultiplyElement *dme = dynamic_cast<MultiplyElement*>(second());
    if (dme)
      bottoms = dme->arguments;
    else bottoms.push_back(second());
    bool cancelled = false;
    for (list<AlgebraElement*>::iterator ti = tops.begin();
	 ti != tops.end();)
    {
      for (list<AlgebraElement*>::iterator bi = bottoms.begin();
	   bi != bottoms.end();)
      {
	if ((**ti) == (**bi))
	{
	  tops.erase(ti++);
	  bottoms.erase(bi++);
	  cancelled = true;
	  continue;
	}
	++bi;
      }
      ++ti;
    }
    if (cancelled)
    {
      if (tops.empty())
	first() = new ScalarElement(1);
      else
      { MultiplyElement *nnme = new MultiplyElement();
        for (list<AlgebraElement*>::iterator ti = tops.begin();
	     ti != tops.end(); ++ti)
	  nnme->addArgument(*ti);
	delete first();
        first() = nnme->normalize();
      }
      if (bottoms.empty())
	second() = new ScalarElement(1);
      else
      { MultiplyElement *ndme = new MultiplyElement();
        for (list<AlgebraElement*>::iterator bi = bottoms.begin();
	     bi != bottoms.end(); ++bi)
	  ndme->addArgument(*bi);
	delete second();
        second() = ndme->normalize();
      }
    }
    // now if the bottom has reduced to a number...
    ScalarElement *dse = dynamic_cast<ScalarElement*>(second());
    if (dse)
    {
      MultiplyElement *pme
	= new MultiplyElement(*first(),
			      ScalarElement(1.0/dse->value()));
      delete this;
      return pme->normalize();
    }
    return this;
  }
  virtual AlgebraElement *deepcopy() const
  { DivideElement *ne = new DivideElement;
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      ne->arguments.push_back((*it)->deepcopy());
    return ne;
  }
  virtual vnl_real_polynomial as_polynomial()
  { cerr << "Polynomial division not implemented!!" << endl;
    return first()->as_polynomial();
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { cerr << "Polynomial division not implemented!!" << endl;
    return first()->as_npolynomial(va);
  }
  virtual ~DivideElement()
  {} 
};

// first ^ second
class PowerElement : public BinaryOperatorElement
{
public:
  PowerElement() {}
  PowerElement(const AlgebraElement &base, const AlgebraElement &exponent)
  { setArguments(base,exponent);
  } 
  unsigned precedence() const
  { return 4; } 
  string operatorString(bool code=false) const
  { return " ^ "; } 
  virtual string operatorName() const
  { return "Power"; }
  AlgebraElement *normalize()
  { first() = first()->normalize();
    second() = second()->normalize();
    // distribute over products
    MultiplyElement *mb = dynamic_cast<MultiplyElement*>(first());
    if (mb)
    {
      for (list<AlgebraElement*>::iterator fit = mb->arguments.begin();
	   fit != mb->arguments.end(); ++fit)
      {
	AlgebraElement*factor = new PowerElement(**fit,*second());
	*fit = factor->normalize();
      }
      return mb;
    }
    // simplify case of 2 numbers
    ScalarElement *sb = dynamic_cast<ScalarElement*>(first());
    ScalarElement *se = dynamic_cast<ScalarElement*>(second());
    if (sb && se)
    { ScalarElement *ns = new ScalarElement(pow(sb->value(),se->value()));
      delete this;
      return ns;
    }
    return this;
  }
  virtual AlgebraElement *deepcopy() const
  { PowerElement *ne = new PowerElement;
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      ne->arguments.push_back((*it)->deepcopy());
    return ne;
  }
  virtual vnl_real_polynomial as_polynomial()
  { // only can do positive int powers of single variables
    VariableElement*vb = dynamic_cast<VariableElement*>(first());
    ScalarElement  *se = dynamic_cast<ScalarElement*>(second());
    if ( !(vb && se) )
    { cerr << "Complicated powers in polynomials not implemented!!"
	   << endl;
      return ScalarElement(0).as_polynomial();
    }
    vnl_real_polynomial p = vb->as_polynomial();
    vnl_real_polynomial accum = p;
    for (unsigned i = unsigned(se->value()); i > 1; --i)
      accum = accum * p;
    return accum;
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { // only can do int powers of single variables
    VariableElement*vb = dynamic_cast<VariableElement*>(first());
    ScalarElement  *se = dynamic_cast<ScalarElement*>(second());
    if ( !(vb && se) )
    { cerr << "Complicated powers in polynomials not implemented!!"
	   << endl;
      return ScalarElement(0).as_npolynomial(va);
    }
    // there's only one exponent in the matrix; multiply it
    vnl_real_npolynomial p = vb->as_npolynomial(va);
    p.polyn() *= unsigned(se->value());
    return p;
  }
  virtual ~PowerElement()
  {} 
};

class FunctionElement : public NaryOperatorElement
{
public:
  string _name;
  FunctionElement() {}
  FunctionElement(string fn) : _name(fn) {}
  FunctionElement(string fn, const AlgebraElement&a0) : _name(fn)
  { addArgument(&a0); }
  FunctionElement(string fn,
		  const AlgebraElement&a0,
		  const AlgebraElement&a1) : _name(fn)
  { addArgument(&a0);
    addArgument(&a1);
  }
  unsigned precedence() const
  { return 10; } 
  string operatorString(bool code=false) const
  { return _name; 
  }
  string asString(bool code=false) const
  {
    string as;
    for (list<AlgebraElement*>::const_iterator ai = arguments.begin();
	 ai != arguments.end(); ++ai)
    {
      if (as.size() > 0)
	as += ", ";
      as += (*ai)->asString(code);
    }
    return operatorString(code) + '(' + as + ')';
  }
  AlgebraElement *normalize()
  { return this;
  }
  virtual AlgebraElement *deepcopy() const
  { FunctionElement *ne = new FunctionElement;
    for (list<AlgebraElement*>::const_iterator it = arguments.begin();
	 it != arguments.end(); ++it)
      ne->arguments.push_back((*it)->deepcopy());
    return ne;
  }
  string sortKey() const
  { return "3 "+sortKey_recursive(); }
  virtual vnl_real_polynomial as_polynomial()
  { cerr << "No functions in polynomials! For Christ's sake!" << endl;
    return vnl_real_polynomial(0.0);
  }
  virtual vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { cerr << "No functions in polynomials! For Christ's sake!" << endl;
    return ScalarElement(0).as_npolynomial(va);
  }
  virtual ~FunctionElement() {}
};

// inline AlgebraElement &operator+(const AlgebraElement&a,const AlgebraElement&b)
// {
//   return *new AddElement(a,b);
// //   AddElement sum(a,b);
// //   return sum.normalize();
// }

// inline AlgebraElement &operator-(const AlgebraElement&a)
// {
//   NegateElement neg(a);
//   return neg.normalize();
// }

// inline AlgebraElement &operator-(const AlgebraElement&a,const AlgebraElement&b)
// {
//   return a + (-b);
// }

// inline AlgebraElement &operator*(const AlgebraElement&a,const AlgebraElement&b)
// {
//   MultiplyElement product(a,b);
//   return product.normalize();
// }

// inline AlgebraElement &operator/(const AlgebraElement&a,const AlgebraElement&b)
// {
//   DivideElement quotient(a,b);
//   return quotient.normalize();
// }

// a number of specializations of this compare

inline ostream&operator<<(ostream&o, AlgebraElement&el)
{
  return o << el.asString();
}

// now the expression
class AlgebraExpression
{
public:
  AlgebraElement *element;
  void setElement(AlgebraElement *el)
  { if (element)
      delete element;
    element = el ? el->normalize() : 0;
  }
  AlgebraExpression(AlgebraElement *el=0)
    : element(el?el->normalize():0)  {}
  AlgebraExpression(AlgebraElement &el)
    : element(el.normalize()) {}
  AlgebraExpression(int i)
    : element(new ScalarElement(double(i)))    {}
  AlgebraExpression(unsigned int i)
    : element(new ScalarElement(double(i)))    {} 
  AlgebraExpression(double d)
    : element(new ScalarElement(d))            {} 
  AlgebraExpression(const string &s)
    : element(new VariableElement(s))          {}
  AlgebraExpression(const char *cp)
    : element(new VariableElement(string(cp))) {}
  AlgebraExpression(AlgebraExpression const &that)
    : element(that.element ? that.element->deepcopy() : 0)
  {}

  ~AlgebraExpression()
  { if (element) delete element; }

  AlgebraExpression &operator=(AlgebraExpression const&that)
  { if (that.element) setElement(that.element->deepcopy());
    return *this;
  }
  AlgebraExpression &operator=(double d)
  { setElement(new ScalarElement(d));
    return *this;
  }
  AlgebraExpression &operator=(int d)
  { setElement(new ScalarElement(double(d)));
    return *this;
  }
  AlgebraExpression &operator=(const string &s)
  { setElement(new VariableElement(s));
    return *this;
  }
  AlgebraExpression &operator=(const char *cp)
  { return operator=(string(cp)); }
  string asString() const
  {
    return element->asString();
  }
  string asCode() const
  {
    return element->asString(true);
  }
  string prettyPrint() const
  {
    return element->prettyPrint();
  }
  operator double() const
  { return element->asDouble(); }
  AlgebraExpression operator+(const AlgebraExpression &that) const
  {
    return AlgebraExpression(new AddElement(*element,*that.element));
  }
  AlgebraExpression operator-() const
  {
    return AlgebraExpression(new NegateElement(*element));
  }
  AlgebraExpression operator-(const AlgebraExpression &that) const
  {
    return *this + (-that);
  }
  AlgebraExpression operator*(const AlgebraExpression &that) const
  {
    return AlgebraExpression(new MultiplyElement(*element,*that.element));
  }
  AlgebraExpression operator/(const AlgebraExpression &that) const
  {
    return AlgebraExpression(new DivideElement(*element,*that.element));
  }
  AlgebraExpression &operator+=(const AlgebraExpression&r)
  { setElement(new AddElement(*element,*r.element));
    return *this;
  }
  AlgebraExpression &operator-=(const AlgebraExpression&r)
  { return (*this = *this-r);
  }
  AlgebraExpression &operator*=(const AlgebraExpression&r)
  { return (*this = *this * r);
  }
  AlgebraExpression &operator/=(const AlgebraExpression&r)
  { return (*this = *this/r);
  }
  bool operator==(const AlgebraExpression &that) const
  { return *element == *that.element;
  }
  bool operator!=(const AlgebraExpression &that) const
  { return !(*this==that);
  }
  bool operator!=(int i) const
  { return !(*this==AlgebraExpression(i));
  }
  bool operator<(const AlgebraExpression &that) const
  { return *element < *that.element;
  }
  /*
  bool operator<=(const AlgebraExpression &that) const
  { return !(that < *this);
  }
  bool operator>(const AlgebraExpression &that) const
  { return that < *this;
  }
  bool operator>=(const AlgebraExpression &that) const
  { return !(*this < that);
  }
  */
  AlgebraExpression substitute(VariableSubstitution&sub)
  { return AlgebraExpression(element->substitute(sub));
  }
  vnl_real_polynomial as_polynomial()
  { return element->as_polynomial(); }
  vnl_real_npolynomial as_npolynomial(VariableAssignment &va)
  { return element->as_npolynomial(va); }
};

inline ostream&operator<<(ostream &o, const AlgebraExpression &x)
{
  return o << x.asString();
}

inline istream&operator>>(istream &i, AlgebraExpression &x)
{ // forget it
  return i;
}

// you can use this with apply() (see vnl_extensions.h)
//  to convert a vnl_matrix<double> to vnl_matrix<AlgebraExpression>
//  and use it in algebra.
inline AlgebraExpression algebraize(const double &d)
{
  return AlgebraExpression(d);
}

// for vnl
VCL_DEFINE_SPECIALIZATION 
class vnl_numeric_traits<AlgebraExpression>
{
public:
//: Additive identity
  static const AlgebraExpression zero;
//: Multiplicative identity
  static const AlgebraExpression one;
//: Maximum value which this type can assume
  static const AlgebraExpression maxval;
//: Return value of abs()
  typedef AlgebraExpression abs_t;
//: Name of a type twice as long as this one for accumulators
//  and products.
  typedef AlgebraExpression double_t;
//as in vnl_rational.h:
//  Name of type which results from multiplying this type with a
//  double.  Note that this requires an explicit cast from double to
//  vnl_rational.  This must be a built-in type: do not set this to
//  vnl_rational, since that would require std::sqrt(vnl_rational)
//  etc., which is not allowed.
   typedef double real_t;
   //  typedef AlgebraExpression real_t;
};

// let's suppose these things are real, not complex
VCL_DEFINE_SPECIALIZATION struct vnl_complex_traits<AlgebraExpression>
{
  enum { isreal = true };
  static AlgebraExpression conjugate(const AlgebraExpression &x)
  { return x; }
  // here you go but by god, don't use it
  static vcl_complex<AlgebraExpression>
    complexify(const AlgebraExpression &x)
  { return vcl_complex<AlgebraExpression>(x, AlgebraExpression(0.0)); }
};

inline bool vnl_math_isfinite(const AlgebraExpression &x)
{ return true; } 

inline bool vnl_math_isnan(const AlgebraExpression &x)
{ return false; } 
/*
namespace std {
  //inline AlgebraExpression sqrt(const AlgebraExpression &x)
  inline AlgebraExpression sqrt(AlgebraExpression x)
  {
    return AlgebraExpression(
	     new FunctionElement(string("sqrt"),*x.element));
  }

//   inline AlgebraExpression pow(const AlgebraExpression &base,
// 			       const AlgebraExpression &exponent)
//   {
//     return AlgebraExpression(
// 	     new PowerElement(*base.element,*exponent.element));
//   }
  
  // can't actually do angles because vnl requires an actual double
  inline double acos(const AlgebraExpression &x)
  {
    return 0;
  }
}
*/

inline AlgebraExpression vnl_math_squared_magnitude(AlgebraExpression x)
{ return x*x; }

inline AlgebraExpression vnl_math_abs(AlgebraExpression x)
{
  return AlgebraExpression(
	   new FunctionElement(string("abs"),*x.element));
}

#endif//ALGEBRA_H
