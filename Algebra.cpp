#include "Algebra.h"

const AlgebraExpression
 vnl_numeric_traits<AlgebraExpression>::zero = 0.0;
const AlgebraExpression
 vnl_numeric_traits<AlgebraExpression>::one = 1.0;
const AlgebraExpression
 vnl_numeric_traits<AlgebraExpression>::maxval = HUGE;

AlgebraElement *ScalarElement::normalize()
{
  if (_value < 0)
  {
    _value = - _value;
    AlgebraElement *nel = new NegateElement(this);
    return nel->normalize();
  }
  return this;
}

AlgebraElement *VariableElement::substitute(VariableSubstitution&sub)
{ VariableSubstitution::iterator found = sub.find(_name);
  if (found == sub.end())
    return deepcopy();
  else
    return found->second.element->deepcopy();
}

AlgebraElement *NegateElement::normalize()
{
  argument = argument->normalize();
  if (NegateElement*narg = dynamic_cast<NegateElement*>(argument))
  {
    AlgebraElement *r = narg->argument;
    narg->argument = 0;
    delete narg;
    argument = 0;
    delete this;
    return r;
  }
  // if we want to distribute - over sums
  if (AddElement*sum = dynamic_cast<AddElement*>(argument))
  {
    for (list<AlgebraElement*>::iterator ti = sum->arguments.begin();
	 ti != sum->arguments.end(); ++ti)
      *ti = new NegateElement(**ti);
    argument = 0;
    delete this;
    return sum->normalize();
  }
  return this;
}

AlgebraElement *MultiplyElement::normalize()
{
  { // this sorts arguments
    AlgebraElement *nm = NaryOperatorElement::normalize();
    if (nm != this)
      return nm;
  }
  // collect all negative signs
  int sign = 1;
  for (list<AlgebraElement*>::iterator it = arguments.begin();
       it != arguments.end(); ++it)
  {
    NegateElement *n =
      dynamic_cast<NegateElement*>(*it);
    if (n)
    {
      sign *= -1;
      *it = n->argument;
      n->argument = 0;
      delete n;
    }
  }
  if (arguments.size() > 0) // should be nonzero!
  {
    // special dealing with numbers
    ScalarElement *no1 =
      dynamic_cast<ScalarElement*>(*arguments.begin());
    // this assumes if there are any numbers one of them
    // will be the first argument
    if (no1)
    {
      // consolidate numbers into one
      for (list<AlgebraElement*>::iterator it = ++arguments.begin();
	   it != arguments.end(); )
      {
	ScalarElement *z =
	  dynamic_cast<ScalarElement*>(*it);
	if (z)
	{
	  *no1 = (*no1 * *z);
	  arguments.erase(it++);
	}
	else
	  ++it;
      }
      // if the answer is one
      if (aboutEqual(no1->_value,1))
      { // suppress it if it multiplies something
	if (arguments.size() > 1)
	  arguments.erase(arguments.begin());
      }
      // if it's zero, zero is what you get.
      else if (aboutEqual(0,no1->_value))
      {
	delete this;
	return new ScalarElement(0);
      }
    }
    if (arguments.size() == 1)
    {
      AlgebraElement *el = arguments.front();
      arguments.clear();
      delete this;
      if (sign == -1)
	return new NegateElement(el);
      return el;
    }
    // if we want to distribute into fractions...
    if (arguments.size() > 0 && 1)
    {
      for (list<AlgebraElement*>::iterator fi = arguments.begin();
	   fi != arguments.end(); ++fi)
	if (DivideElement*frac = dynamic_cast<DivideElement*>(*fi))
	{
	  // take fraction out of list of factors
	  arguments.erase(fi);
	  // put all the remaining factors into the numerator
	  arguments.push_back(frac->first()->deepcopy());
	  frac->setArguments(*this,*frac->second()->deepcopy());
	  if (sign==-1)
	    return (new NegateElement(*frac))->normalize();
	  else
	    return frac->normalize();
	}
    }
    // if we want to distribute over sums...
    if (arguments.size() > 0 && 0)
    {
      for (list<AlgebraElement*>::iterator fi = arguments.begin();
	   fi != arguments.end(); ++fi)
	if (AddElement*sum = dynamic_cast<AddElement*>(*fi))
	{
	  // take sum out of list of factors
	  arguments.erase(fi);
	  // multiply the remaining factors into each term of sum
	  for (list<AlgebraElement*>::iterator ti
		 = sum->arguments.begin();
	       ti != sum->arguments.end(); ++ti)
	    *ti = new MultiplyElement(*deepcopy(),**ti);
	  if (sign==-1)
	    return new NegateElement(*sum->normalize());
	  else
	    return sum->normalize();
	}
    }
    if (sign == -1)
      return new NegateElement(this);
  }
  return this;
}

bool compareArguments::operator()( const AlgebraElement* a,
				   const AlgebraElement* b ) const
{ return a->sortKey() < b->sortKey();
}

#include <vnl/vnl_c_vector.txx>
VNL_C_VECTOR_INSTANTIATE_ordered(AlgebraExpression);
#include <vnl/vnl_vector.txx>
VNL_VECTOR_INSTANTIATE(AlgebraExpression);
#include <vnl/vnl_matrix.txx>
VNL_MATRIX_INSTANTIATE(AlgebraExpression);
