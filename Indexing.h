/* -*- C++ -*- */
#ifndef INDEXING_H
#define INDEXING_H
//
// Indexing and VectorAccess classes defined here
//
#include <map>
//using std::map;
#include <algorithm>
//using std::pair;
#include <vector>
//using std::vector;
#include <iostream>
using namespace std;
using std::ostream;

template<typename I>
class _Index;
template<typename I>
class _Indexing;

// typically we use int, but this could change.
typedef int index_type;
typedef _Index<index_type> Index;
typedef _Indexing<index_type> Indexing;

// unique key factory is common to all indexings regardless
//  of index type
class _Indexing_base
{
public:
  typedef int unique_key_type;
  static const unique_key_type flag_key = -1;
  static void stagger(int ini, int off)
  { init = ini; offset = off; }
protected:
  // stagger for possible mpi use
  static int init, offset;
  static unique_key_type new_key(void)
  { static unique_key_type next_key = -1;
    if (next_key == -1)
      next_key = init;
    unique_key_type ret = next_key;
    next_key += offset;
    return ret;
  }
};

// this is a particular indexing scheme
//  associates numbers to unique keys and v.v.
template <typename I>
class _Indexing : public _Indexing_base
{
public:
  typedef I index_type;
  static const I flag_index;
protected:
  map<I,unique_key_type> map_i_k;
  map<unique_key_type,I> map_k_i;
public:
  _Indexing() {}
  // lookup functions
  _Index<I> lookup(I i)
  { typename map<I,unique_key_type>::iterator ii = map_i_k.find(i);
    if (ii != map_i_k.end())
      return _Index<I>(ii->first,this);
    else
      return _Index<I>(flag_index,this);
  }
  unique_key_type key(const _Index<I> &i)
  { return i.key();
  }
  unique_key_type key(I i)
  { typename map<I,unique_key_type>::iterator ii = map_i_k.find(i);
    return (ii!=map_i_k.end()) ? ii->second : flag_key;
  } 
  I index(unique_key_type k)
  { typename map<unique_key_type,I>::iterator ki = map_k_i.find(k);
    return (ki!=map_k_i.end()) ? ki->second : flag_index;
  } 
  // this is used when the indexing is given an index of a different type
  template<class II>
  I index(const _Index<II> &i)
  {
    return index(i.key());
  }
  // this specializes the above when the indexing and index have same type
//   I index(const _Index<I> &i)
//   { if (i.indexing == this || i.i == flag_index)
//       return i.i;
//     else
//       return index(i.key());
//   }
  _Index<I> convert(const _Index<I> &i)
  { return i;
  }

  template<class II>
  bool knowsIndex(const _Index<II> &i)
  { return (index(i.key()) >= 0);
  }
//   { if (i.indexing == this)
//       return i;
//     else
//       return _Index<I>(index(i),*this);
//   }
  // setup functions
  // insert an index, with or without associations to other indexings
  void registerIndex(I i, unique_key_type k = new_key())
  { typename map<I,unique_key_type>::iterator ii = map_i_k.find(i);
    if (ii != map_i_k.end())
    { if (ii->second != k)
        cerr << "registerIndex(): " << i
	     << " already associated with " << ii->second
	     << " -- not " << k << endl;
      return;
    }
    map_i_k[i] = k;
    map_k_i[k] = i;
    //print(map_i_k);
  }
  template<typename II>
  void registerIndex(I i, const _Index<II> &ix)
  { registerIndex(i, ix.key());
  }
  // associate indexes in 2 indexings together
  template<typename II>
  static void associate(I i0, _Indexing<I> &ix0, II i1, _Indexing<II> &ix1)
  {
    unique_key_type k0 = ix0.key(i0), k1 = ix1.key(i1);
    if (k0 != flag_key && k1 != flag_key)
    { if (k0 != k1)
      {
	cerr << "can't assocate " << i0 << " (indexing " << &ix0
	     << ") with " << i1 << " (indexing " << &ix1
	     << ") -- already mapped in both" << endl;
	return;
      } 
    }
    else if (k0 != flag_key)
      ix1.registerIndex(i1,k0);
    else if (k1 != flag_key)
      ix0.registerIndex(i0,k1);
    else
    {
      ix0.registerIndex(i0);
      ix1.registerIndex(i1,ix0.key(i0));
    } 
  }
  template<typename II>
  void associate(I i0, II i1, _Indexing<II> &ix1)
  { associate(i0, *this, i1, ix1);
  }
  void reindex(I was, I is)
  { typename map<I,unique_key_type>::iterator ii = map_i_k.find(was);
    if (ii == map_i_k.end())
    { cerr << "can't reindex " << was << " -- not indexed" << endl;
      return;
    } // else
    if (map_i_k.find(is) != map_i_k.end())
    { cerr << "can't reindex " << was << " to " << is << " -- already in use"
	   << endl;
      return;
    }
    if (was == is) return;
    unique_key_type k = ii->second;
    map_i_k.erase(ii);
    map_i_k.insert(make_pair(is,k));
    typename map<unique_key_type,I>::iterator ki = map_k_i.find(k);
    if (ki==map_k_i.end())
      cerr << was << " not found in " << this << "map_k_i" << endl;
    ki->second = is;
  }
  void removeIndex(I was)
  { typename map<I,unique_key_type>::iterator ii = map_i_k.find(was);
    if (ii == map_i_k.end())
    { cerr << "can't remove " << was << " -- not indexed" << endl;
      return;
    } // else
    unique_key_type k = ii->second;
    map_i_k.erase(ii);
    map_k_i.erase(k);
  }
  unsigned int size()
  {
    return map_i_k.size();
  }
  void clear()
  {
    map_i_k.clear();
    map_k_i.clear();
  }
};

// altered implementation of Index, that stores key instead of index,
// introduced 9/25/2006 LW, because in the other implementation an
// index can get out of date, for instance when the Integrator reindexes
template <typename I>
class _Index
{
public:
  typedef _Indexing_base::unique_key_type unique_key_type;
  unique_key_type _key;
  //I i;
  //_Indexing<I> *indexing;

  // default constructor gives 'null index'
  _Index(): _key(-1) {} 

  // construction from a given key shouldn't be implicit
  explicit _Index(unique_key_type k): _key(k) {} 

  _Index(I ii, _Indexing<I> &xx) : _key(xx.key(ii)) {}
  _Index(I ii, _Indexing<I> *xx) : _key(xx->key(ii)) {}

  _Index(const _Index &other) : _key(other._key) {}
  _Index &operator=(const _Index &other)
  { _key = other._key; return *this; }

  static _Index nullIndex()
  { return _Index(-1); } 

  bool isValid()
  { return (_key >= 0);
  }

  typename _Indexing<I>::unique_key_type key(void) const
  { return _key;
  }

  bool operator==(const _Index<I>&other) const
  { return key()==other.key(); } 
  bool operator!=(const _Index<I>&other) const
  { return !(*this==other); }
};

// virtual base class for VectorAccess'es of various kinds
template<typename T>
class VectorAccess
{
public:
  virtual ~VectorAccess() {}
  virtual T &operator[] (int i)=0;
  virtual const T &operator[] (int i) const=0;
  virtual unsigned int size(void) const = 0;
  virtual void resize(unsigned int ns) = 0;
  void print(ostream&os=cout) const
  { os << *this << endl; }
};

template <typename T>
void copyTo(const VectorAccess<T> *ivp, VectorAccess<T> *vp)
{
  if (ivp == vp) return;
  unsigned int n = ivp->size();
  if (n != vp->size()) // some classes don't want resize called
    vp->resize(n);
  for (unsigned int i = 0; i < n; i++)
    (*vp)[i] = (*ivp)[i];
}

template <typename T>
class c_vectorAccess : public VectorAccess<T>
{
private:
  int n;
  T *x;
public:
  c_vectorAccess(int nn, T *xx) : n(nn), x(xx) {}
  virtual ~c_vectorAccess() {}
  virtual T &operator[] (int i)
  { return x[i]; }
  virtual const T &operator[] (int i) const
  { return x[i]; }
  virtual unsigned int size(void) const
  { return n; }
  virtual void resize(unsigned int ns)
  { cerr << "error: c_vectorAccess doesn't resize" << endl;
    return;
  }
};

typedef c_vectorAccess<double> c_dvectorAccess;

template <typename T>
class stl_vectorAccess : public VectorAccess<T>
{
private:
  vector<T> &_v;
public:
  stl_vectorAccess(vector<T>&sv) : _v(sv) {}
  virtual ~stl_vectorAccess() {}
  // if one of these gives you a 'returning reference to temporary' error
  // it's because you don't have -D_BVECTOR_H in your CFLAGS
  T &operator[] (int i)
  { return _v[i]; }
  const T &operator[] (int i) const
  { return _v[i]; }
   unsigned int size(void) const
  { return _v.size(); }
  void resize(unsigned int ns)
  { return _v.resize(ns); }
  // Used by IndexedVector.  You should not use.
  vector<T> &v()
  { return _v; }  
};

typedef stl_vectorAccess<double> stl_dvectorAccess;

// IndexedVector uses this weird thing
template<typename T>
class VectorAccessAccess : public VectorAccess<T>
{
protected:
public:
  VectorAccess<T> &va;
public:
  VectorAccessAccess(VectorAccess<T>&inner) : va(inner) {}
  virtual ~VectorAccessAccess() {}
  virtual T &operator[] (int i)
  { return va[i]; }
  virtual const T &operator[] (int i) const
  { return va[i]; }
  virtual unsigned int size(void) const
  { return va.size(); }
  virtual void resize(unsigned int ns)
  { return va.resize(ns); }
};

template <typename T, typename I=index_type>
class IndexedVector : public VectorAccessAccess<T>
{
  using VectorAccessAccess<T>::va;
public:
  _Indexing<I> &indexing;
private:
  static VectorAccess<T> &createV(void)
  { vector<T>*vv = new vector<T>;
    stl_vectorAccess<T>*vx = new stl_vectorAccess<T>(*vv);
    return *vx;
  }
  bool ownv;
  //using VectorAccessAccess<T>::va;
public:
  IndexedVector(VectorAccess<T> &x, _Indexing<I> &i) :
    VectorAccessAccess<T>(x), indexing(i), ownv(false) {}
  IndexedVector(_Indexing<I> &i) :
    VectorAccessAccess<T>(createV()), indexing(i), ownv(true) {}

  template<class AccessClass, class VectorClass>
  IndexedVector(VectorClass&v, _Indexing<I> &i) :
    VectorAccessAccess<T>(AccessClass(v)), indexing(i), ownv(false) {} 

  virtual ~IndexedVector()
  { if (ownv)
    { stl_vectorAccess<T>*vx = static_cast<stl_vectorAccess<T>*>(&va);
      delete &(vx->v());
      delete vx;
    }
  }
  virtual T &operator[](_Index<I> i)
  { return va[indexing.index(i)]; }
  virtual const T &operator[](_Index<I> i) const
  { return va[indexing.index(i)]; }
  // shouldn't need to repeat these but ... ?
  virtual T &operator[] (I i)
  { return va[i]; }
  virtual const T &operator[] (I i) const
  { return va[i]; }
};

// for MatrixAccess and IndexedMatrix see Matrix.h

// diagnostic
template<typename T>
ostream&operator<<(ostream&o,const VectorAccess<T>&vx)
{
  for (int i = 0; i < (int)vx.size(); i++)
    o << (i>0?" ":"") << vx[i];
  return o;
}

template <class _Key, class _Tp>
void print(map<_Key,_Tp> &h)
{
  for(typename map<_Key,_Tp>::iterator
	it = h.begin(); it != h.end(); it++)
  {
    pair<_Key,_Tp> p = *it;
    cout << p.first << ' ' << p.second << endl;
  }
}

template <class _Key, class _Tp>
void print(map<_Key,_Tp> *h)
{ print(*h); }

#endif//INDEXING_H
