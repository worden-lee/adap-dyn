/*
 * Lotka-Volterra evolutionary ecology program
 *
 * 5/9/2000  Lee Worden
*/
#include "Community.h"
#include "Communicator.h"
#include "Integrator.h"
#include "OutputController.h"
#include "Site.h"
#include "Node.h"
#include "util.h"
#include "rand.h"
#include "Parameters.h"
#include <stdlib.h>
#include <string.h> // memset

/* default constructor: set up initial community
*/
Community::Community(void)
  : nX(0), _nSpecies(0), _nDead(0),
    _alive_a(_alive), alive(_alive_a,_speciesIndexing),
    _mutationTries(0), site(0)
{}

Community::Community(Community &other)
  : nX(other.nX), _nSpecies(other._nSpecies), _nDead(other._nDead),
    _speciesIndexing(other._speciesIndexing), _alive(other._alive),
    _alive_a(_alive), alive(_alive_a,_speciesIndexing),
    _mutationTries(0), site(0)
{}

// ensure all vectors are allocated to fit _nSpecies, nX
void Community::checkAllocation(void)
{
  unsigned int as = alive.size();
  if (_nSpecies>0 && as == 0)
    as = 1;
  while ( as < (unsigned)_nSpecies )
    as *= 2;
  if (as != alive.size())
    alive.resize(as);
}

// assumes this is called after the whole site is set up
//  this is generally going to be overridden
void Community::initialize(void)
{
  nX = _nSpecies = 0;
  checkAllocation();
}

void Community::siteIs(Site *s)
{
  site = s;
}

Indexing& Community::variableIndexing()
{ return site->integrator->variableIndexing(); }

bool Community::allDead(void)
{ if (speciesCount() == 0)
  { site->node->outputcontroller->log("No more species\n");
    return true;
  }
  return false;
}

// given x values at time t, provide the derivatives
void Community::calcDerivatives(double t, const VectorAccess<double>*x,
				VectorAccess<double>*dxdt)
{
}

void Community::calcNextState(double t, const VectorAccess<double>*x,
			      VectorAccess<double>*nextx)
{
}

bool Community::shouldGoExtinct(const Index &i, const VectorAccess<double> &x)
{
  int ii = site->integrator->variableIndexing().index(i);
  if (isVariableAPopulation(i) && isVariableInUse(i) &&
      x[ii] < parameters.extinctionThreshold())
    return true;
  else
    return false;
}

Index Community::insertNewSpecies(bool __alive)
{
  _speciesIndexing.registerIndex(_nSpecies);
  return insertNewSpecies(Index(_nSpecies,_speciesIndexing), __alive);
}

Index Community::insertNewSpecies(Index existingIndex, bool __alive)
{
  int nS = _nSpecies++;
  //cerr << "insertNewSpecies " << nS << endl;
  nX++;
  checkAllocation();
  _speciesIndexing.registerIndex(nS,existingIndex.key());
  alive[existingIndex] = __alive;
  if (__alive && site->integrator)
    site->integrator->newVariable(existingIndex);
  return existingIndex;
}

// when a local extinction occurs the species might reappear
// due to immigration from another site
void Community::extinction(const Index &n)
{
  int ni = speciesIndexing().index(n);
  if ( ni < 0 && site->outputcontroller )
    site->outputcontroller->
      log("Community::extinction passed invalid index\n");
  else if ( alive[ni] )
  {
    alive[ni] = false;
    _nDead++;
  }
}

void Community::possiblyReindex(void)
{
  // arbitrary rule but seems pretty good
  if ( _nSpecies >= 8 && _nDead > _nSpecies / 4 )
    reindex();
}

// reindexes only species
void Community::reindex(void)
{
  //site->outputcontroller->log("reindex community\n");
  int nAlive = 0;
  for ( int n = 0; n < _nSpecies; n++ )
  {
    Index ni(n,_speciesIndexing);
    if (isVariableAPopulation(ni))
    {
      if(alive[n])
      {
	if (nAlive != n)
	  reindexSpecies(n,nAlive);
	nAlive++;
      }
      else
	retireSpecies(n);
    }
  }
  nX += nAlive - _nSpecies;
  _nSpecies = nAlive;
  _nDead = 0;
  checkAllocation(); // shrink vectors, if you want to
}

void Community::reindexSpecies(int olds, int news)
{
  _speciesIndexing.reindex(olds,news);
  alive[news] = alive[olds];
}

void Community::retireSpecies(int olds)
{
  _speciesIndexing.removeIndex(olds);
}

bool Community::isVariableInUse(const Index &n)
{ // assuming it's a population
  int i = speciesIndexing().index(n);
  return (i!=-1) && alive[i];
}

bool Community::isVariableAPopulation(const Index &n)
{
  return true; // they are all, until further notice
}

bool Community::isVariableInUse(int n)
{ return isVariableInUse(Index(n,variableIndexing())); } 
bool Community::isVariableAPopulation(int n)
{ return isVariableAPopulation(Index(n,variableIndexing())); }
bool Community::isVariableASpeciationCandidate(int n)
{ return isVariableASpeciationCandidate(Index(n,variableIndexing())); }

void Community::doSpeciation(void)
{
  Index parent, daughter;
  bool goodenough = true;
  bool ateq = site->integrator->atEquilibrium();
  unsigned maxtries = parameters.mutationTriesBeforeESS();
  
  do
  {
    parent = chooseParentForSpeciation();
    if (!parent.isValid())
      return;
    // convert parent index because variable indexes might change
    //  during speciation?
    parent = speciesIndexing().convert(parent);
      
    site->integrator->increaseEvolutionaryTime();
    // this modifies the state
    daughter = speciate(parent);
    ++_mutationTries;

    if (parameters.useInvasionProbability())
    { // this doesn't work away from equilibrium
      if (!site->integrator->atEquilibrium())
	continue;
      // now that it's constructed, figure out whether
      //  it's good enough to invade
      double I = site->integrator->invasionFunction(daughter);
      double threshold =
	parameters.invasionProbabilityFactor() * I;
      VectorAccess<double> &x = site->integrator->state(); 
      int pi = site->integrator->variableIndexing().index(parent);
      if (parameters.usePopulationInInvasionProbability())
	threshold *= x[pi];
      //if (parameters.outputODEMessages)
      //site->outputcontroller->log("invasion function %g\n", I);
      //site->outputcontroller->log("invasion probability %g\n", threshold);
      goodenough = (unif_distn() < threshold);
      if (!goodenough)
      {
	//site->outputcontroller->log("reject mutant\n");
	int di = site->integrator->variableIndexing().index(daughter);
	x[pi] += x[di];
	x[di] = 0;
	//site->integrator->setState(t,&x); 
	// community may notify the outputcontroller
	extinction(daughter);
	site->integrator->extinction(daughter);
      }
    }  // if called when we're not at equilibrium, we try once,
       // and keep count in _mutationTries
  }  while(!goodenough && ateq && _mutationTries < maxtries);

  if (goodenough)
  {
    _mutationTries = 0;
    //    if (!site->integrator->currentlyAtESS())
    postSpeciation(parent,daughter);
  }
  else if (_mutationTries >= maxtries)
  {
    site->outputcontroller->log("couldn't find a viable mutant "
				"after trying %d candidate%s\n",
				_mutationTries,
				_mutationTries>1?"s":"");

    _mutationTries = 0;
    site->outputcontroller->log("Reached ESS\n");
    //site->outputcontroller->logCurrentState();
    //site->outputcontroller->recordCommunity();
    site->integrator->reachedESS();
    site->outputcontroller->ess(site->integrator->time());
  }
}

Index Community::chooseParentForSpeciation(void)
{
  double biomass;
  int i;
  VectorAccess<double> &x = site->integrator->state(); 
  // choose parent species i by biomass
  for ( i = 0, biomass = 0;
	i != (signed)x.size(); ++i )
    if ( isVariableInUse(i)
	 && isVariableASpeciationCandidate(i)
	 && ( x[i] > parameters.introductionSize() ) )
      biomass += x[i];
  if ( biomass <= parameters.introductionSize() )
  {
    cerr << "can't speciate -- not enough biomass!\n";
    return Index(-1,_speciesIndexing);
  }
  double z = unif_distn() * biomass;
  for ( i = 0; i != (int)x.size(); ++i )
    if ( isVariableInUse(i)
	 && isVariableASpeciationCandidate(i)
	 && ( x[i] > parameters.introductionSize() ) )
    {
      z -= x[i];
      if ( z <= 0 )
	break;
    }
  if ( i == (int)x.size() )
  {
    cerr << "Error in speciation!\n";
    site->outputcontroller->log(
      "Error in speciation (chose %i'th variable as parent)!\n", i);
  }
  return Index(i,site->integrator->variableIndexing());
}

void Community::postSpeciation(const Index &parent, const Index &daughter)
{
  site->outputcontroller->speciation(site->integrator->time(),
				     parent, daughter);
  site->outputcontroller->recordCommunity();
}

void Community::doImmigration(VectorAccess<double> *x)
{
  Index migrant = createImmigrant();
  if ( site->integrator )
  {
    site->integrator->immigration(migrant);
    site->outputcontroller->immigration(site->integrator->time(),
				   migrant);
  }
}

// this should be overridden, it's just to make the compiler happy
Index Community::createImmigrant()
{
  return insertNewSpecies();
}

Index Community::speciate(const Index &parent)
{
  Index daughter = insertNewSpecies();
  return daughter;
}

double Community::perturb(double orig, double step, bool bounded)
{
  double repl;
  double mut;
  
  if (bounded)
    mut = uniform_double(-step, step);
  else
    mut = norml_distn(0,step);

  repl = orig + mut;
  
  return repl;
}

/* diffusion of a given value doesn't happen if the amount to diffuse
  is below the 'grain size' -- to prevent underflow (?) 
*/
double Community::amountOfDiffusion(const Index &index,
				    const VectorAccess<double>*x)
{
  double total = (*x)[site->integrator->variableIndexing().index(index)];
  double drip = total * parameters.diffusionConstant();
  if ( drip >= parameters.extinctionThreshold() )
    return drip;
  else
    return 0;
}    

void Community::calcDiffusion(const VectorAccess<double> *x,
			      VectorAccess<double> *d)
{
  d->resize(x->size());
  for ( int i = 0; i != (int)x->size(); ++i )
    (*d)[i] = amountOfDiffusion(Index(i,site->integrator->variableIndexing()),
				x);
}

void print(Community &comm)
{
  cout << comm;
}

void Community::recordCommunity(ostream &o)
{
  o << "cx =\n{ ";
  if ( site->integrator )
  {
    o << "t -> " << site->integrator->time() << ",\n  ";
    o << "t_ev -> " << site->integrator->evolutionaryTime() << ",\n  ";
  }
  o << "variables -> {";
  { bool fi = true;
    for ( int i = 0; i < nX; i++ )
      if ( isVariableInUse(i) )
	o << (fi? ((fi=false),""):", ")
	  << site->outputcontroller->basename(Index(i,variableIndexing()));
  }
  o << "}";
  if ( site->integrator )
  {
    VectorAccess<double> &x = site->integrator->state();
    o << ",\n  state -> {";
    {
      bool fi = true;
      for ( int i = 0; i != (int)x.size(); ++i )
	if ( isVariableInUse(i) )
	{ o << (fi ? ((fi=false),""): ", ") << x[i]; }
    }
    o << "}";
  }
  addToPrintForMathematica(o);
  o << "\n}\n";
}

ostream& operator<< (ostream &sr, Community &comm)
{
  sr << "variables: ";
  bool comma=false;
  for ( long i = 0; i < comm.nX; i++ )
    if ( comm.isVariableInUse(i) ) 
      sr << (comma? (comma=true),"":", ")
	 << comm.site->outputcontroller->
	      basename(Index(i,comm.variableIndexing()))
	 << " (" << i << ")";
  sr << ")";
  return sr;
}

