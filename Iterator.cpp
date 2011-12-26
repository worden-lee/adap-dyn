#include "Iterator.h"
#include "util.h"
#include "Node.h"
#include "Indexing.h"
#include <errno.h>

Iterator::Iterator() : stateaccess(statevector)
{
}
  
VectorAccess<double> &Iterator::state()
{ return stateaccess; }

void Iterator::setState(double t0, VectorAccess<double> *s)
{ t = t0;
  copyTo(s, &stateaccess);
}

void Iterator::integrateNonstop(double t1)
{
  //  double dtsav = parameters.outputTimeStep;
  //  long nok, nbad;

  if ( site->community->nSpecies() > 0 )
  {
    // do speciation, whatever necessary to get to far end
    while ( t < t1 )
    {
      double eps = parameters.speciateAtEquilibrium() ?
	DMIN(1e-6,parameters.equilibriumThreshold()/100.) :
	1e-6;

      //      odeint_mut(t1, eps, 100, 0,
      //		 EXTINCTION_THRESHOLD, &nok, &nbad, dxsav);

      /*------------replace integration with iteration-----------------*/
      
      iterate_to(t1, eps, parameters.extinctionThreshold()/*, dtsav*/);
      
#ifndef MPI
      if ( site->community->allDead() )
      {
        break;
      }
#endif
      if ( parameters.doSpeciation() )
      {
	if ( timeToSpeciate ) // ( speciationThreshold <= 0 )
	{
	  // includes checkSize()
	  site->community->doSpeciation();
	  lastEventWas = t;
	  resetSpeciationThreshold();
	}
      }    // end if doSpeciation
    }   // end while t < t1
    site->outputcontroller->flush();
  }  // end if totalBiomass > 0
}

void Iterator::iterate_to(double t2, double eps, 
			  double extinctionthreshold/*,
			  double dtsav*/)
{
  double t1 = t;
  errno = 0;
  stl_dvectorAccess sx(statevector);
  while(t<t2)
  {
    vector<double> prev = statevector;
    const stl_dvectorAccess px(prev);
    site->community->calcNextState(t, &px, &sx);
    t++;

    // check for extinctions and explosions, and figure out when to
    // speciate: either by counting down speciationThreshold for poisson
    // process, or by looking for equilibrium
    {
      bool expl = false;
      // for speciateAtEquilibrium
      bool stillChanging = false;
      double max_rel = 0;
      //      site->community.calc_derivs(t,statevector,dxdt);
      
      for(int i = 0; i != (int)statevector.size(); ++i)
      {
	if ( site->community->isVariableInUse(i) 
	     && site->community->isVariableAPopulation(i) )
	{
	  if ( parameters.doExtinction() &&
	       ( statevector[i] < extinctionthreshold ) &&
	       site->community->isVariableAPopulation(i) )
	  {
	    Index ix(i,indexing);
	    site->community->extinction(ix);
	    site->outputcontroller->extinction(t,ix); 
	    statevector[i] = 0;
	    if ( parameters.speciateAtEquilibrium() )
	      lastEventWas = t;
	  }
	  else if ( statevector[i] > parameters.explosionCeiling() )
	  {
	    Index ix(i,indexing);
	    site->node->outputcontroller->logWithCommunity(
	      "Explosion of species %s -- removed\n",
	      site->outputcontroller->basename(ix).c_str());
	    site->community->extinction(ix);
	    statevector[i] = 0;
	    expl = true;
	  if ( parameters.speciateAtEquilibrium() )
	    lastEventWas = t;
	  }
	  if ( parameters.doSpeciation() &&
	       !parameters.speciateAtEquilibrium() )
	    speciationThreshold -= statevector[i] /* * hdid */;
	} // end of if ( i.isAPopulation() )

	if ( parameters.doSpeciation() && parameters.speciateAtEquilibrium() )
	{ // check for motion of all variables
	  double rel = fabs((statevector[i] - prev[i]) / statevector[i]);
	  if (rel > max_rel) max_rel = rel;
	  if ( rel > parameters.equilibriumThreshold() )
	    stillChanging = true;
	}
      } // end of for
      
      if ( parameters.doSpeciation() )
      {
	if ( parameters.speciateAtEquilibrium() )
	{
	  // site->node->outputcontroller->log << "t = " << t << ", max_rel = "
	  //	 << max_rel << ", stillChanging = " << stillChanging << "\n";
	  if ( stillChanging )
	    noChangeSince = t;
	  else if ( t - noChangeSince >= parameters.equilibriumTime() )
	  {
	    site->outputcontroller->logEquilibrium();
	    timeToSpeciate = true;
	  }
	  if ( t - lastEventWas >= parameters.maxWaitForEquilibrium() )
	  {
	    site->node->outputcontroller->logWithCommunity(
		 "Tired of waiting for equilibrium...\n");
	    site->outputcontroller->logEquilibrium();
	    timeToSpeciate = true;
	  }
	}   // end of if speciateAtEquilibrium
        else
	{
	  if ( speciationThreshold <= 0 )
	    timeToSpeciate = true;
	}
      }   // end of if doSpeciation
    
      if (expl)
	site->outputcontroller->recordCommunity();
    }   // end of bool expl;

    //    if (dtsav > 0 && fabs(t-tsav) >= fabs(dtsav)) 
    {
      site->outputcontroller->recordValues(&state(),t);
      //      tsav=t;
    }

    if ((t-t2)*(t2-t1) >= 0.0 || timeToSpeciate)
      break;
  } // end of while t < t1
  return;
}

// don't use
void Iterator::integratePartially(double t1)
{}

double Iterator::invasionFunction(Index &ix)
{
  int i = indexing.index(ix);
  vector<double> ns(statevector.size());
  stl_dvectorAccess nsx(ns), sx(statevector);
  site->community->calcNextState(site->integrator->time(),&sx,&nsx);
  return (ns[i] - statevector[i]) / statevector[i]; // dt = 1
}
