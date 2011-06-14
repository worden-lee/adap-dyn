#include "Integrator.h"
#include "Parameters.h"
#include "Site.h"
#include "rand.h"

Integrator::Integrator()
 : site(0), indexing(), t(0), evoTime(0),
   noChangeSince(0), lastEventWas(0), timeToSpeciate(false),
   exploded(false), equil(false), tiredOfWaiting(false), atESS(false)
{}

Integrator::Integrator(const Integrator &in)
  : site(0), indexing(in.indexing), t(0), evoTime(0),
    noChangeSince(0), lastEventWas(0), timeToSpeciate(false),
    exploded(false), equil(false), tiredOfWaiting(false), atESS(false)
{}

void Integrator::siteIs(Site *c)
{
  site = c;
}

// bool Integrator::atEquilibrium(void)
// {
//   return ( t - noChangeSince >= parameters.equilibriumTime() );
// }

void Integrator::increaseEvolutionaryTime() 
{
  evoTime++;
  if (site->outputcontroller)
    site->outputcontroller->step(t);
}

void Integrator::newVariable(Index new_index)
{
  if (!indexing.knowsIndex(new_index))
  { int ni = nextIndex();
    indexing.registerIndex(ni,new_index);
    //  set(new_index, 0);
  }
}

void Integrator::speciation(Index pi, Index di)
{
  //newVariable(di);
  set( di, parameters.introductionSize() );
  set( pi, state(pi) - parameters.introductionSize() );
}

void Integrator::immigration(Index new_index)
{
  //  newVariable(new_index);
  set( new_index, parameters.introductionSize() );
}

// sum of all variables that are known to be populations
//  and that are not extinct
//double Integrator::totalBiomass(void)
//{ 
//  VectorAccess<double> &s = state();
//  return site->community->totalPopulation(&s);
//}

double Integrator::invasionFunction(Index &ix)
{
  int i = indexing.index(ix);
  VectorAccess<double> &x = state(); 
  vector<double> dxdt(x.size());
  stl_dvectorAccess da(dxdt);
  site->community->calcDerivatives(time(),&x,&da);
  return dxdt[i] / x[i];
}

#if VNL
#include "vnl_extensions.h"
vnl_matrix<double> Integrator::jacobian(const VectorAccess<double> *ixp)
{
  vector<double> x;
  stl_dvectorAccess xp(x);
  copyTo(ixp,&xp);
  int n = x.size();
  vnl_matrix<double> J(n,n);
  vnl_vector<double> f0(n), f1(n);
  vector<double> f0v(n), f1v(n);
  vnl_dvectorAccess f0x(f0), f1x(f1);
  stl_dvectorAccess f0vx(f0v), f1vx(f1v);
    // @@ bug here - if it discovers a resource, it
    // wants to resize the fx0 - should add an exception
    // maybe?
    // hopefully addressed by using stl vectors here
    // and double checking
  site->community->calcDerivatives(t,&xp,&f0vx);
  if (f0vx.size() != f0x.size())
    return vnl_matrix<double>(0,0);
  copyTo(&f0vx,&f0x);
  for(int j=0;j<n;j++)
  {
    const double eps=0.01;
    x[j]+=eps;
    // @@ here too
    site->community->calcDerivatives(t,&xp,&f1vx);
    if (f1vx.size() != f1x.size())
      return vnl_matrix<double>(0,0);
    copyTo(&f1vx,&f1x);
    x[j]-=eps;
    J.set_column(j,(f1-f0)/eps);
  }
  return J;
}
#endif

void Integrator::resetSpeciationThreshold(void)
{
  if ( parameters.doSpeciation() )
  {
    timeToSpeciate = false;
    if ( parameters.speciateAtEquilibrium() )
      noChangeSince = t;
    else
      speciationThreshold = exponl_distn( 1/parameters.speciationRate() );
    // speciationRate is in events per unit biomass per unit time
  }
}

void Integrator::resetNextImmigration(void)
{
  if ( parameters.doImmigration() )
  {
    nextImmigration = t + exponl_distn( 1/parameters.immigrationRate() );
    // immigrationRate is in events per unit time
  }
}
