#include "NRIntegrator.h"
#include "Site.h"
#include "Node.h"
#include "Community.h"
#include "OutputController.h"
#include "util.h"
#include "rand.h"
#include <errno.h>
#include <stdlib.h>
#include <math.h>

void NRIntegrator::setState(double tx, VectorAccess<double> *vp)
{
  VectorAccess<double> &v = *vp;
  t = tx;
  copyTo(&v,&svx);
  resetSpeciationThreshold();
  resetNextImmigration();
  exploded = false;
}

VectorAccess<double> &NRIntegrator::state(void)
{
  return svx;
}

void NRIntegrator::integratePartially(double t1)
{
  //  double dtsav = parameters.outputTimeStep;
  long nok, nbad;
  errno = 0;
  
  double eps = parameters.speciateAtEquilibrium() ?
    DMIN(1e-6,parameters.equilibriumThreshold()/100.) :
    1e-6;

  // odeint returns early if it's time for a speciation or immigration
  odeint(t1, eps, 100, 0, &nok, &nbad/*, dtsav*/);
}

// does not handle diffusion, calling function must do that after t1
void NRIntegrator::integrateNonstop(double t1)
{
  //  if ( totalBiomass() > 0 )
  {
    // do speciation, whatever necessary to get to far end
    while ( t < t1 )
    {
      integratePartially(t1);

      if (site->outputcontroller)
        site->outputcontroller->flush();
#ifndef MPI
      if ( site->community->allDead())
      {
        break;
      }
      if ( hadExplosion() )
      {
        break;
      }
#endif
      if ( parameters.doSpeciation() && timeToSpeciate )
	 // ( speciationThreshold <= 0 )
      {
	site->community->doSpeciation();
	lastEventWas = t;
	resetSpeciationThreshold();
      }    // end if doSpeciation
      if ( parameters.doImmigration() && t >= nextImmigration )
      {
	site->community->doImmigration(&svx);
	lastEventWas = t;
	resetNextImmigration();
      }
    }   // end while t < t1
  }  // end if biomass > 0
}

// void print(const StateVector &x)
// {
//   cout << "(";
//   for (VariableIndex i = x.firstIndex(); i != x.finalIndex(); ++i )
//     cout << (i != x.firstIndex() ? ", " : "") 
// 	 << (site->outputcontroller ? site->outputcontroller->basename(i) : "X(?)")
// 	 << "=" << x[i];
//   cout << ")\n";
// }

/*
 * odeint()
 *
 * This is from Numerical Recipes:
 *  integrate ODE from x1 to x2
 *
 * But this one includes extinction, speciation and immigrations
 *
 * modified to use double instead of float
 * and use Community object to calculate derivatives instead of function *
 * and StateVector objects instead of arrays of doubles
 * and member variables to record state instead of lots of reference arguments
 *
 * also xp, yp et al. aren't provided, they're allocated.
 */
#define NRANSI
#define MAXSTP 10000
#define TINY 1.0e-30

// if we are speciating by poisson process, we decrease speciationThreshold
// by the number of accumulated mass-time units, and flag timeToSpeciate
// when we hit zero.  If we are speciating at equilibrium, we flag 
// timeToSpeciate when we have been within EQUILIBRIUM_THRESHOLD for
// EQUILIBRIUM_TIME.  We return with either ( timeToSpeciate == true )
// or ( t == t2 ).  t always records the current time of the integrator's
// state.
void NRIntegrator::odeint(double t2, double eps, 
			double h1, double hmin, 
			long *nok, long *nbad/*,
			double dtsav*/)
{
  int i;
  long nstp;
  double hnext,hdid,h;
  long nvar = statevector.size();
  static vector<double> dxdt, xscal;
  stl_dvectorAccess dxdtx(dxdt), xsx(xscal);
  double t1 = t;

  dxdt.resize(nvar);
  xscal.resize(nvar);

  // there is code in this function that assumes t2 > t1 (h > 0)
  h=SIGN(h1,t2-t1);
  *nok = (*nbad) = 0;
  
  //if (dtsav > 0) tsav=t-dtsav*2.0;
  for (nstp=0;nstp<MAXSTP;nstp++) {
    site->community->calcDerivatives(t,&svx,&dxdtx);
    for (i = 0; i != (int)statevector.size();++i)
      xscal[i]=fabs(statevector[i])+fabs(dxdt[i]*h)+TINY;
    //    if (dtsav > 0 && fabs(t-tsav) > fabs(dtsav)) 
    {
      if (site->outputcontroller)
	site->outputcontroller->recordValues(&svx,t);
      //      tsav=t;
    }
    if ((t+h-t2)*(t+h-t1) > 0.0) h=t2-t;
    if ( fabs(h) > parameters.outputTimeStep() )
      h = copysign(parameters.outputTimeStep(),h);
    // rkqs changes t
    rkqs(statevector,dxdt,h,eps,xscal,&hdid,&hnext);
    // check for extinctions and explosions, and figure out when to
    // speciate: either by counting down speciationThreshold for poisson
    // process, or by looking for equilibrium
    {
      bool stillChanging = false;
      //if ( parameters.doSpeciation && parameters.speciateAtEquilibrium )
      {
	site->community->calcDerivatives(t,&svx,&dxdtx);
      }

      // personhours-=hdid; //Not accounting for biomass
      for(i = 0; i != (int)statevector.size(); ++i)
      {
	if ( site->community->isVariableAPopulation(i)
	     && site->community->isVariableInUse(i) )
	{
	  if ( parameters.doSpeciation() &&
	       !parameters.speciateAtEquilibrium() )
	  {
	    speciationThreshold -= statevector[i] * hdid;
	  }
      
          if ( parameters.doExtinction() &&
	       ( statevector[i] < parameters.extinctionThreshold() ) )
	  {
	    Index ix(i,indexing);
	    site->community->extinction(ix);
	    if (site->outputcontroller)
	      site->outputcontroller->extinction(t,ix); 
	    extinction(ix);
	    statevector[i] = 0;
	    if ( parameters.speciateAtEquilibrium() )
	    {
	      lastEventWas = t;
	    }
	  }
	  // This should be outside of the check for Population.
	  else if ( statevector[i] > parameters.explosionCeiling() )
	  {
	    Index victim(i,indexing);
	    if (site->outputcontroller)
	      site->outputcontroller->logWithCommunity(
	        "%g Explosion of species %s\n", t,
	        site->outputcontroller->basename(victim).c_str());
	    exploded = true;
	    if (site->outputcontroller)
	      site->outputcontroller->explosion(t, victim);
	    return;
	  }
	  // This check also should be outside of the Population check.
	  //if ( parameters.doSpeciation && parameters.speciateAtEquilibrium )
	  { // check for motion of all variables
	    double rel = fabs(dxdt[i] / statevector[i]);
	    if ( rel > parameters.equilibriumThreshold() )
	      stillChanging = true;
	  }
	} // end of if ( isAPopulation() and isVariableInUse() )
      } // end of for

      if ( parameters.doSpeciation() )
      { 
	if ( parameters.speciateAtEquilibrium() )
	{
	  if ( stillChanging )
	    noChangeSince = t;
	  else if ( t - noChangeSince >= parameters.equilibriumTime() )
	  {
	    if (site->outputcontroller)
	      site->outputcontroller->equilibrium(t);
	    timeToSpeciate = true;
	    noChangeSince = t;
	  }
	  if ( t - lastEventWas >= parameters.maxWaitForEquilibrium() )
	  {
	    if (site->outputcontroller)
	    { site->outputcontroller->
		logWithCommunity("%g Tired of waiting for equilibrium...\n", t);
	      // the right thing to do?
	      site->outputcontroller->logEquilibrium();
	    }
	    timeToSpeciate = true;
	    noChangeSince = t;
	  }
	}
	else
	{
	  if ( speciationThreshold <= 0 )
	    timeToSpeciate = true;
	}
      }
    }   // end of bool expl;
    if (hdid == h) ++(*nok); else ++(*nbad);
      // stop if at max time or threshold of personhours.
    if ((t-t2)*(t2-t1) >= 0.0 ||
	(parameters.doSpeciation() && timeToSpeciate) ||
	(parameters.doImmigration() && t >= nextImmigration))
      break;
    if (fabs(hnext) <= hmin) ierror("Step size too small in odeint");
    h=hnext;
  }
  if (nstp >= MAXSTP)
    ierror(""); // just print the time
  //ierror("Too many steps in routine odeint");

  //  if (dtsav > 0)
  //  {
  //    for (i=0;i<nR;i++) outputcontroller.recordValue(true,i,t,x[i]);
  //    for (i=0;i<nN;i++) outputcontroller.recordValue(false,i,t,x[i+nR]);
  //  }

  return;
}
#undef MAXSTP
#undef TINY
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 3^03. */

/*
 * rkqs()
 *
 * This is fifth-order Runge-Kutta code (with adaptive stepsize)
 * pp. 719ff. of Numerical Recipes in C (2nd Ed.)
 *
 * changed all floats to double, arrays to vectors,
 *  all names from (y,x) to (x,t)
 *
 * Lee Worden
 */
#include <math.h>
#define NRANSI
#include "Integrator.h"
#include "util.h"

#define SAFETY 0.9
#define PGROW -0.2
#define PSHRNK -0.25
#define ERRCON 1.89e-4

void NRIntegrator::ierror(const char *er)
{
  if (site->outputcontroller)
    site->outputcontroller->log("%g %s\n",t,(char*)er);
}
 
void NRIntegrator::rkqs(vector<double> &x, vector<double> &dxdt, 
		      double htry, double eps,
		      vector<double> &xscal, double *hdid, double *hnext)
{
  int i;
  double errmax,h,htemp,tnew;
  //  StateVector xerr(x.nVars()), xtemp(x.nVars());
  static vector<double> xerr, xtemp;
  xerr.resize(x.size());
  xtemp.resize(x.size());
  
  h=htry;
  for (;;) {
    rkck(x,dxdt,h,xtemp,xerr);
    errmax=0.0;
    for ( i=0; i != (int)x.size(); ++i )
      errmax=DMAX(errmax,fabs(xerr[i]/xscal[i]));
    errmax /= eps;
    if (errmax > 1.0) {
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h=(h >= 0.0 ? DMAX(htemp,0.1*h) : DMIN(htemp,0.1*h));
      tnew=t+h;
      if (tnew == t) ierror("stepsize underflow in rkqs");
      continue;
    } else {
      if (errmax > ERRCON)
	*hnext=SAFETY*h*pow(errmax,PGROW);
      else *hnext=5.0*h;
      t += (*hdid=h);
      for ( i=0; i != (int)x.size(); ++i )
	x[i]=xtemp[i];
      break;
    }
  }
}
#undef SAFETY
#undef PGROW
#undef PSHRNK
#undef ERRCON
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 3^03. */

/*
 * rkck()
 *
 * This is fifth-order Runge-Kutta code (with adaptive stepsize)
 * pp. 719ff. of Numerical Recipes in C (2nd Ed.)
 *
 * I changed arrays to StateVectors, changed all floats to double,
 * I think I changed all names from (y,x) to (x,t)
 *
 * Lee Worden
 */
#define NRANSI
#include "Integrator.h"
#include "Community.h"
#include "Site.h"

void NRIntegrator::rkck(vector<double> &x, vector<double> &dxdt, 
		      double h, vector<double> &xout, vector<double> &xerr)
{
  int i;
  static double a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
    b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2,
    b51 = -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
    b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
    b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
    c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0,
    dc5 = -277.00/14336.0;
  double dc1=c1-2825.0/27648.0,dc3=c3-18575.0/48384.0,
    dc4=c4-13525.0/55296.0,dc6=c6-0.25;
  //  StateVector ak2(x.nVars()), ak3(x.nVars()), ak4(x.nVars()),
  //    ak5(x.nVars()), ak6(x.nVars()), xtemp(x.nVars());
  static vector<double> ak2, ak3, ak4, ak5, ak6, xtemp; // lots faster
  static stl_dvectorAccess
    ak2x(ak2), ak3x(ak3), ak4x(ak4), ak5x(ak5), ak6x(ak6), xtx(xtemp);
  long nv = x.size();
  ak2.resize(nv);
  ak3.resize(nv);
  ak4.resize(nv);
  ak5.resize(nv);
  ak6.resize(nv);
  xtemp.resize(nv);

  for ( i=0; i != (int)x.size(); ++i )
    xtemp[i]=x[i]+b21*h*dxdt[i];
  site->community->calcDerivatives(t+a2*h,&xtx,&ak2x);
  for ( i=0; i != (int)x.size(); ++i )
    xtemp[i]=x[i]+h*(b31*dxdt[i]+b32*ak2[i]);
  site->community->calcDerivatives(t+a3*h,&xtx,&ak3x);
  for ( i=0; i != (int)x.size(); ++i )
    xtemp[i]=x[i]+h*(b41*dxdt[i]+b42*ak2[i]+b43*ak3[i]);
  site->community->calcDerivatives(t+a4*h,&xtx,&ak4x);
  for ( i=0; i != (int)x.size(); ++i )
    xtemp[i]=x[i]+h*(b51*dxdt[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
  site->community->calcDerivatives(t+a5*h,&xtx,&ak5x);
  for ( i=0; i != (int)x.size(); ++i )
    xtemp[i]=x[i]+h*(b61*dxdt[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
  site->community->calcDerivatives(t+a6*h,&xtx,&ak6x);
  for ( i=0; i != (int)x.size(); ++i )
    xout[i]=x[i]+h*(c1*dxdt[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
  for ( i=0; i != (int)x.size(); ++i )
    xerr[i]=h*(dc1*dxdt[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
}
#undef NRANSI
/* (C) Copr. 1986-92 Numerical Recipes Software 3^03. */

