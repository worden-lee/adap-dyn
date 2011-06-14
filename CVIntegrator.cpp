#include "CVIntegrator.h"
#include "util.h"
#include "Node.h"
#include <exception>
#include <errno.h>
#include "CVodeWithStoppingCondition.h"

class DenseMatAccess : public MatrixAccess<realtype>
{
public:
  DenseMatAccess(DenseMat &m) : d(m), data(m->data) {}
  virtual ~DenseMatAccess() {}
  unsigned int rows(void) const
  { return d->M; } 
  unsigned int cols(void) const
  { return d->N; }
  void resize(unsigned nr, unsigned nc)
  {} // this is not the resize you are looking for
  realtype &entry(const long &r, const long &c)
  { return DENSE_ELEM(d,r,c);
  }
private:
  DenseMat &d;
  realtype **data;
};

// exception that gets raised when the integration routine needs to
// be interrupted and restarted, i.e. when the number of dynamic variables
// changes.
class CVIntegrationException : public exception {};

// to be precise: when it's time for an extinction, immigration or
// speciation, check() tells l_cvode() to return early (don't throw
// the exception because it destroys the integrating we've done so
// far). When there's a 'discovery' of a variable that needs to start
// being integrated explicitly, calcDerivatives() is allowed to go
// ahead and add that variable to the integrator, then we throw and
// redo some integrating if necessary.
static int f(realtype t, N_Vector y, N_Vector ydot, void *f_data)
{
  N_VectorAccess yn(&y), ydn(&ydot);
  ((CVAccess*)f_data)->calcDerivatives(&yn, &ydn);
//   if (((CVAccess*)f_data)->check(&yn,t))
//     return 1;
// does throwIfNeeded() throw an exception, or return nonzero to
//  stop integration?
  if (((CVAccess*)f_data)->throwIfNeeded())
    return -1;
  return 0;
}  

/*  currently not enabled
// this is a CVDenseJacFn
static void Jac(long int N, DenseMat J, realtype t,
		N_Vector y, N_Vector fy, void *jac_data,
		N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static void Jac(long int N, DenseMat J, RhsFn f, void *f_data, realtype t,
                N_Vector y, N_Vector fy, N_Vector ewt,
		realtype h, realtype uround,
                void *jac_data, long int *nfePtr, N_Vector vtemp1,
                N_Vector vtemp2, N_Vector vtemp3)
{
  N_VectorAccess ya(&y);
  DenseMatAccess ja(J);
  ((CVAccess*)f_data)->jacobian(&ya, f, &ja);
}
*/

static int Check(N_Vector x, realtype t, void *f_data)
{ N_VectorAccess xx(&x);
  return ((CVAccess*)f_data)->check(&xx,t);
}

void CVAccess::calcDerivatives(const VectorAccess<double> *x,
			       VectorAccess<double> *dxdt)
{
  in.site->community->calcDerivatives(in.time(),x,dxdt);
}

// trying a change for new cvode:
//  instead of throw, return negative signifying an unrecoverable
//  error
int CVAccess::throwIfNeeded(void)
{
  if (in.discovery)
  { //cout << "throw\n" << flush;
    throw CVIntegrationException();
    return -1;
  }
  return 0;
}

int CVAccess::check(const VectorAccess<double> *xp, realtype t)
{
  return in.check(xp,t);
}

// this shouldn't modify outside objects, call to outputcontroller, etc
// except recordValues, I guess.
int CVIntegrator::check(const VectorAccess<double> *xp, realtype t0)
{
  t = t0;
//  cout << t << endl << *xp << endl;

  bool stillChanging = false;
  bool forgo = false;
  const VectorAccess<double> &x    = *xp;
  //const VectorAccess<double> &dxdt = *dxdtp;
  vector<double> dxdt;
  stl_dvectorAccess dxdtx(dxdt);
  unsigned int xs = x.size();
  site->community->calcDerivatives(t,xp,&dxdtx);
  if (x.size() != xs)
    return 1;

  // discovery can be set by calcDerivatives()
  if (discovery)
    forgo = true;

  // go on and check for extinctions before returning

  for (unsigned int i = 0; i < x.size(); ++i)
  {
    if ( site->community->isVariableInUse(i) )
    {
      bool ispop = site->community->isVariableAPopulation(i);

      if (ispop)
      {
	if ( parameters.doSpeciation() && !parameters.speciateAtEquilibrium() )
	{
	  // this stuff doesn't work without further work
	  //speciationThreshold -= statevector[i] * hdid;
	}
      } // end of if ( isAPopulation() )

      // allow it to consider extinction for all, not just populations
      if ( parameters.doExtinction() &&
	   ( x[i] < parameters.extinctionThreshold() ) &&
	   /*(x[i] > 0) &&*/ (dxdt[i] <= 0) &&
	   site->community->shouldGoExtinct(Index(i,indexing), x))
      {
	//cout << "check: " << t << " extinc" << endl;
	noChangeSince = t;
	// doing an extinction will disrupt the state, so we have to
	// return nonzero to pop us out of CVode() and then do it
	extnctn = true;
	forgo = true;
      }

      if ( nvx[i] > parameters.explosionCeiling() )
      {
	//cout << "check: " << t << " expl" << endl;
	exploded = true;
	forgo = true;
      }

      // check for motion of all variables
      // could be more sophisticated (2 thresholds, for instance)
      double slope = (ispop ? dxdt[i] / x[i] : dxdt[i]);
      if ( fabs(slope) > parameters.equilibriumThreshold() )
	stillChanging = true;
    } // end of if ( isVariableInUse() )
  } // end of for
  
  if ( parameters.doSpeciation() )
  { 
    if ( parameters.speciateAtEquilibrium() )
    {
      if ( stillChanging )
      {
	noChangeSince = t;
	equil = false;
      }
      else if ( t - noChangeSince >= parameters.equilibriumTime() )
      {
	//cout << "check: " << t << " equil" << endl;
	equil = true;
	timeToSpeciate = true;
	noChangeSince = t;
	forgo = true;
      }

      if ( !equil &&
	   t - lastEventWas >= parameters.maxWaitForEquilibrium() )
      {
	//cout << "check: " << t << " tired" << endl;
// 	site->outputcontroller->
// 	  logWithCommunity("(last event was %g)\n", lastEventWas); 
	equil = false;
	tiredOfWaiting = true;
	timeToSpeciate = true;
	noChangeSince = t;
	forgo = true;
      }
    }
    else
    {
      if ( speciationThreshold <= 0 )
      { timeToSpeciate = true;
	noChangeSince = t;
        forgo = true;
      }
    }
  }
  if (site->outputcontroller)
    site->outputcontroller->recordValues(xp,t);

  return forgo? 1 : 0; // nonzero if should stop integrating
}

void CVIntegrator::setState(double nt, VectorAccess<double> *vx)
{
  t = nt;
  copyTo(vx,&nvx);
  noChangeSince = nt;
  lastEventWas = nt;
  timeToSpeciate = exploded = equil = atESS = false;
  restartCVODE = true;
}

VectorAccess<double> &CVIntegrator::state(void)
{
  return nvx;
}

void CVIntegrator::integratePartially(double t1)
{
  extnctn = false;
  exploded = false;
  discovery = false;
  atESS = false;
  errno = 0;

  ++integrating;

  // have we resized? if so invalidate the cvode_mem
  if (nvx.size()!=cvode_mem_N || restartCVODE)
  { cvode_mem_N = nvx.size();
    if (cvode_mem)
    {
      CVodeFree(&cvode_mem);
      cvode_mem = 0;
      restartCVODE = false;
    }
  } 
  // init or reinit cvode_mem if needed
  if (!cvode_mem) 
  {
    if (reltol==0)
    { // if eps is too small everything goes to hell (?)
      const double eps=1e-6;//1e-12;
      reltol = parameters.speciateAtEquilibrium() ?
	DMIN(eps, parameters.equilibriumThreshold()/100.) :
	eps;
      abstol = DMIN(parameters.extinctionThreshold()/100.,eps);
//    for (int k = 0; k < OPT_SIZE; k++)
// 	switch(k)
// 	{
// 	case MXSTEP:
// 	  iopt[k] = 100000;//1000;
// 	  break;
// 	default:
// 	  iopt[k] = 0;
// 	  break;
// 	}
    }
    if ( !(cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON)) )
    { // if there's an error, it printed a message
      return;
    }
//     cvode_mem = CVodeMalloc(cvode_mem_N, f, t, nv, BDF, NEWTON,
// 			    SS, &reltol, &abstol, &ax,
// 			    NULL, TRUE, iopt, NULL, NULL);
    switch ( CVodeMalloc(cvode_mem, f, t, nv, CV_SS, reltol, &abstol) )
    {
    case CV_MEM_NULL:
      printf("null pointer passed to CVodeMalloc.\n");
      return;
    case CV_MEM_FAIL:
      printf("CVodeMalloc failed.\n");
      return;
    case CV_ILL_INPUT:
      printf("CVodeMalloc failed: ill input.\n");
      return;
    case CV_SUCCESS:
    default:
      break;
    }    
    switch(CVodeSetMaxNumSteps(cvode_mem, 100000))
    {
    case CV_MEM_NULL:
      printf("null pointer passed to CVodeSetMaxNumSteps.\n");
      return;
    case CV_ILL_INPUT:
      printf("CVodeSetMaxNumSteps failed: ill input.\n");
      return;
    case CV_SUCCESS:
    default:
      break;
    }    
    // Tell CVODE to use the dense matrix linear solver
    // and optionally a user-supplied Jacobian routine
//     if (ax.doJacobian())
//       CVDense(cvode_mem, Jac, &ax);
//     else 
    switch(CVDense(cvode_mem, cvode_mem_N))
    {
    case CV_MEM_NULL:
      printf("null pointer passed to CVDense.\n");
      return;
    case CV_MEM_FAIL:
      printf("CVDense failed.\n");
      return;
    case CV_ILL_INPUT:
      printf("CVDense failed: ill input.\n");
      return;
    case CV_SUCCESS:
    default:
      break;
    }    
    // pointer to the CVAccess, to be passed into the f() function
    switch(CVodeSetFdata(cvode_mem,&ax))
    {
    case CV_MEM_NULL:
      printf("null pointer passed to CVodeSetFdata.\n");
      return;
    case CV_ILL_INPUT:
      printf("CVodeSetFdata failed: ill input.\n");
      return;
    case CV_SUCCESS:
    default:
      break;
    }    
    // don't write error messages to cerr, let me decide whether
    // to report them
    switch(CVodeSetErrFile(cvode_mem,0))
    {
    case CV_MEM_NULL:
      printf("null pointer passed to CVodeSetErrFile.\n");
      return;
    case CV_ILL_INPUT:
      printf("CVodeSetErrFile failed: ill input.\n");
      return;
    case CV_SUCCESS:
    default:
      break;
    }    
  }

  // // check init state, maybe go straight down to handling
  //if (check(&nvx,t)==0)
  try 
  {
    // int flag;
    double tn = t;
    
    // cout << "entering CVodeWithStoppingCondition():\nt = "
    //      << t << "\nnv = " << nvx << endl;
    
    // Check() keeps t updated, checks for extinctions, updates displays, etc.
    integrateFlag = CVodeWithStoppingCondition(cvode_mem, t1, nv, &tn,
					       CV_NORMAL, Check);
    //integrateFlag = CVode(cvode_mem, t1, nv, &tn, CV_NORMAL);

    // cout << "CVodeWithStoppingCondition() returned"
    //      << integrateFlag << ":\nt = "
    //      << t << "\nnv = " << nvx << endl;

    switch (integrateFlag)
    {
    case CV_SUCCESS: break;
#define MIN_HANDLE(code) \
 case code: printf("CVodeWithStoppingCondition returned "#code"\n");break;
      MIN_HANDLE(CV_ROOT_RETURN);
      MIN_HANDLE(CV_TSTOP_RETURN);
      MIN_HANDLE(CV_MEM_NULL);
      MIN_HANDLE(CV_NO_MALLOC);
      MIN_HANDLE(CV_ILL_INPUT);
    case CV_TOO_MUCH_WORK:
      printf("Too many steps in CVode\nCVIntegrator forges on...\n"); break;
      MIN_HANDLE(CV_TOO_MUCH_ACC);
      MIN_HANDLE(CV_ERR_FAILURE);
      MIN_HANDLE(CV_CONV_FAILURE);
      MIN_HANDLE(CV_LINIT_FAIL);
      MIN_HANDLE(CV_LSETUP_FAIL);
      MIN_HANDLE(CV_LSOLVE_FAIL);
    case -8:// this is what we get when we return 'unrecoverable'
      cout << "-8\n";
      break;
    default:
      printf("CVodeWithStoppingCondition failed, flag=%d.\n", integrateFlag);
      break;//return;
    }
  }
  catch (const CVIntegrationException& cex)
  {
//     cout << "catch! " << t << (extnctn ? " extnctn":"")
// 	 << (discovery? " discovery":"") << endl;
    //cout << "t = " << t << "\nnv = " << nvx << endl;
  }
  if ( discovery )
    restartCVODE = true;

  --integrating;
}

void CVIntegrator::integrateNonstop(double t1)
{
  if (site->outputcontroller)
    site->outputcontroller->recordValues(&nvx,t);

  // do speciation, whatever necessary to get to far end
  while ( t < t1 )
  {
    integratePartially(t1);

    if ( equil )
    { if (site->outputcontroller)
        site->outputcontroller->equilibrium(t);
      site->community->equilibrium(t);
      lastEventWas = t;
    }
    else if ( tiredOfWaiting )
    { site->community->tiredOfWaitingForEquilibrium(t);
      if (site->outputcontroller)
	site->outputcontroller->tiredOfWaitingForEquilibrium(t);
      tiredOfWaiting = false;
    }
    
    if ( extnctn )
    {
      vector<double> dxdt;
      stl_dvectorAccess dxdtx(dxdt);
      site->community->calcDerivatives(t,&nvx,&dxdtx);
      for (int i = 0; i < (int)nvx.size(); ++i)
      {
	Index ix(i,indexing);
	if (/*0 < nvx[i] &&*/ nvx[i] < parameters.extinctionThreshold() &&
	    dxdt[i] <= 0 &&
	    site->community->shouldGoExtinct(ix,nvx))
	{
	  if (site->outputcontroller)
	    site->outputcontroller->extinction(t,ix);
	  site->community->extinction(ix);
	  extinction(ix);
	  lastEventWas = t;
	  restartCVODE = true;
	}
      }
//       if (site->outputcontroller)
// 	site->outputcontroller->logCurrentState();
    }

    if ( timeToSpeciate && parameters.doSpeciation() )
    {
      // includes checkSize()

      // moved into Community::doSpeciation
      //increaseEvolutionaryTime();
							
      site->community->doSpeciation();
      lastEventWas = t;
      timeToSpeciate = false;
      resetSpeciationThreshold();	
      if (atESS && parameters.quitAfterReachingESS())
	break;
    }    // end if doSpeciation, timeToSpeciate

    if ( parameters.doImmigration() && t >= nextImmigration )
    {
      site->community->doImmigration(&nvx);
      lastEventWas = t;
      resetNextImmigration();
    }

#ifndef MPI
    if ( hadExplosion() )
    {
      for (int i = 0; i < (int)nvx.size(); ++i)
      {
	if ( nvx[i] > parameters.explosionCeiling() )
	{
	  Index victim(i,indexing);
	  if (site->outputcontroller)
	    site->outputcontroller->logWithCommunity(
              "%g Explosion of %s\n", t,
	      site->outputcontroller->basename(victim).c_str());
	  if (site->outputcontroller)
	    site->outputcontroller->explosion(t, victim);
	}
      }
      break;
    }
    if ( site->community->speciesCount() == 0 )
    { if (site->outputcontroller)
      {
	site->outputcontroller->log("No more species\n");
	site->outputcontroller->flush();      
      }
      break;
    }
#endif
    if (site->outputcontroller)
      site->outputcontroller->flush();      
  }   // end while t < t1
}

void CVIntegrator::extinction(Index departed)
{ //cout << nvx << endl;
  // move last variable into the open slot and shrink
  Integrator::extinction(departed); // x[i] = 0;
  int you = nvx.size() - 1;
  int dep = indexing.index(departed);
  //cout << "replace " << site->outputcontroller->basename(Index(dep,indexing))
  //     << " (" << dep << ") with ";
  //cout << site->outputcontroller->basename(Index(you,indexing)) << " ("
  //     << you << ")\n";
  if (dep != -1)
  { indexing.removeIndex(dep);
    if (you != dep)
    { indexing.reindex(you,dep);
      nvx[dep] = nvx[you];
    }
  }
  nvx.resize(you);
  //cout << nvx << endl;
}
