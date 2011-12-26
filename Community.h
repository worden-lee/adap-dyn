/* -*- C++ -*-
 *
 * Base class for community object
 *
 * 7/1/2000  Lee Worden
 *
 * Derived Community classes should provide:
 *  initialize()   sets up initial community - call parent function first
 *  calcNextState() or calcDerivatives()
 *  insertNewSpecies() if you use a derived species class
 *  reindexSpecies()
 *
 * rewritten for adap-dyn2,
 * 1/2005    Lee Worden
 * 
*/
#ifndef BASE_COMMUNITY_H
#define BASE_COMMUNITY_H

#include "Indexing.h"
#include <stdio.h>
#include <iostream>
using std::cout;

//#include "Site.h"
class Site;

class Community
{
protected:
  int nX;		// total # of dynamic variables, indexed 0 .. nX-1
  int _nSpecies;	// # of species, may or may not be 1-1 with variables
  int _nDead;		// <= _nSpecies

  Indexing _speciesIndexing;
  vector<unsigned char> _alive; // vector of bool, but the bool thing 
  stl_vectorAccess<unsigned char> _alive_a; // is implemented all tricky
  IndexedVector<unsigned char> alive;
  unsigned int _mutationTries;
public:
  Site *site;

  Community(); // sets up the initial community

  Community(Community&other); // copy constructor
  
  virtual ~Community() {}

  virtual void siteIs(Site *);
  virtual void initialize(void); // call after all objects in site are ready

  Indexing& speciesIndexing()
  { return _speciesIndexing; } 

  // for counting from 0 to nX-1
  //  sharedCommunity type things should provide
  virtual Indexing& variableIndexing();
  
  // adds a variant species spawned from an existing one
  virtual void doSpeciation(void);
  Index chooseParentForSpeciation(void);
  // override this one to add specifics
  virtual Index speciate(const Index &parent);

  // doImmigration like doSpeciation introduces a new species, but
  //  in this case it's not derived from a parent
  void doImmigration(VectorAccess<double> *);
  // override this one
  virtual Index createImmigrant(void);

  // this is how many species are currently represented, including
  // dead ones
  int nSpecies(void) 
  { return _nSpecies; }

  // this might or might not match nSpecies, depending on whether there
  // are non-species variables
  virtual int nVars(void) const
  { return nX; }

  // nSpecies() counts both alive and dead; this counts only alive species
  virtual int speciesCount(void) // { return nX - nDead; }
  { long nc = 0;
    for (long i = 0; i < nX; i++)
      if (isVariableAPopulation(i) && isVariableInUse(i))
        nc++;
    return nc;
  }
  virtual bool allDead(void);
  virtual bool isVariableInUse(const Index &n);
  virtual bool isVariableAPopulation(const Index &n);
  virtual bool isVariableASpeciationCandidate(const Index &n)
  { return isVariableAPopulation(n); }
  // int forms assume variable indexing, for convenience
  bool isVariableInUse(int n);
  bool isVariableAPopulation(int n);
  bool isVariableASpeciationCandidate(int n);

  void printCommunity(void) { cout << *this; cout.flush(); }

  virtual double totalPopulation(const VectorAccess<double>*x)
  { double total = 0;
    for (int i = 0; i < nX; i++)
      if (isVariableInUse(i) && isVariableAPopulation(i))
	total += (*x)[i];
    return total;
  }
  
  // for use with Integrator (continuous time dynamics) 
  virtual void calcDerivatives(double t, const VectorAccess<double>*x,
			       VectorAccess<double>*dxdt);

  // for use with Iterator (discrete time dynamics) 
  virtual void calcNextState(double t,  const VectorAccess<double>*x,
			     VectorAccess<double>*nextx);

  virtual Index insertNewSpecies(bool _alive=true);
  virtual Index insertNewSpecies(Index existingIndex, bool _alive=true);

  // should a given variable be set to zero when it gets small?
  virtual bool shouldGoExtinct(const Index &i, const VectorAccess<double> &x);

  //  for spatial models
  void calcDiffusion(const VectorAccess<double>*x, VectorAccess<double>*d);
  double amountOfDiffusion(const Index &index,
			   const VectorAccess<double>*x);

  // notification from integrator
  virtual void extinction(const Index &index);

  // called at end of doSpeciation
  virtual void postSpeciation(const Index &parent, const Index &daughter);

  virtual void possiblyReindex(void);

  //virtual void printForMathematica(ostream &o);
  virtual void recordCommunity(ostream &o);

  virtual void equilibrium(double t) {}
  virtual void tiredOfWaitingForEquilibrium(double t) {}

  virtual void addToCurrentState(ostream &o) {}

  static double perturb(double orig, double step, bool bounded);
  
  friend ostream& operator<< (ostream &o, Community &comm);

protected:
  // this should begin with ",\n" and end just before "\n};"
  virtual void addToPrintForMathematica(ostream &o) {}

  // override this to make sure your vectors are big enough
  virtual void checkAllocation();

  // override one or the other of these
  virtual void reindex(void);
  virtual void reindexSpecies(int olds, int news);
  virtual void retireSpecies(int olds);
};

#endif //BASE_COMMUNITY_H
