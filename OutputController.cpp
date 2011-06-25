#include "OutputController.h"
#include "Site.h"
#include "Node.h"
#include "Community.h"
#include "Communicator.h"
#include "Integrator.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdarg.h>
#include <unistd.h> // for ftruncate
//#include <glob.h>

//std::map <int, string> Color::colormap;
//int counter::globalCount = 0;
//string counter::countString = "";
int Color::_gpcolors[] =
  { 1, 2, 3, 4, 5, 8 }; 
string Color::_dotcolors[] =
  { "red", "green", "blue", "magenta",
    "cyan", "coral" };
int Color::_nextIndex=0;
int Color::_firstColor=0, Color::_lastColor=5;

SiteOutputController::SiteOutputController()
{
  site = NULL;
}

SiteOutputController::~SiteOutputController()
{
  while (controllers.size() > 0)
  {
    delete controllers.back();
    controllers.pop_back();
  }
}

void SiteOutputController::reopenMemberFiles(void)
{
  {  // try to create the output directory
    (void)mkdir(outputDirectory().c_str(),S_IRWXU|S_IRWXG|S_IRWXO);

    struct stat st;
    if (stat(outputDirectory().c_str(),&st))
      cerr << "[" << parameters.rank() << "] problem with directory '"
	   << outputDirectory() << "'\n";
    else if ( !S_ISDIR(st.st_mode) )
      cerr << "[" << parameters.rank() << "] '" << outputDirectory()
	   << "' is not a directory\n";
  }

  if ( parameters.outputCommunities() )
  {
    string comm(outputDirectory() + DIR_SEP + "communities");
    communities.open(comm.c_str());
    if ( !communities.is_open() )
    {
      cerr << "Failed to open communities file (" << comm << ")\n";
      cerr.flush();
    }
  }
}

void SiteOutputController::openMemberFiles(void)
{
  return reopenMemberFiles();
}

void SiteOutputController::possiblyReopenFiles(double t)
{
  // since the NodeOutputController may not have an integrator attached we tell it
  //  what time it is
  if ( site->node && site->node->outputcontroller
       && site->node->outputcontroller != this )
    site->node->outputcontroller->possiblyReopenFiles(t);
}

string SiteOutputController::basename(Index ix)
{
  ostringstream os;
  os << 'X' << ix.key();
  return os.str();
}

string SiteOutputController::outputDirectory(void)
{
  string od = parameters.outputDirectory();
  if (dir.length())
    od += DIR_SEP + dir;
  return od;
}

string SiteOutputController::outfilename(const Index &i)
{
  return outputDirectory() + DIR_SEP + basename(i);
}

void SiteOutputController::siteIs(Site *s)
{
  site = s;
  initDir();
}

void SiteOutputController::initDir(void)
{
  ostringstream os(dir);
  os << site->siteIndex;
  openMemberFiles();
    
  log("init site %i (row = %i, col = %i)\n",
      site->siteIndex, site->row, site->col);
  flush();
}

void SiteOutputController::writeTag(ostream &o)
{
  if ( parameters.totalGridSize() > 1 )
  {
    o << "[" << parameters.rank() << "][" << site->siteIndex << "] ";
  }
}

void SiteOutputController::logTag(void)
{
  writeTag(site->node->outputcontroller->logfile);
  if ( parameters.logToCout() )
    writeTag(cout);
}

void SiteOutputController::log(const char*fmt, ...)
{
  va_list va;
  va_start(va,fmt);
  logTag();
  logRawV(fmt,va);
  va_end(va);
}

void SiteOutputController::logRaw(const char*fmt, ...)
{
  va_list va;
  va_start(va,fmt);
  logRawV(fmt,va);
  va_end(va);
}

// to be overridden by NodeOutputController
void SiteOutputController::logRawV(const char*fmt, va_list va)
{
  site->node->outputcontroller->logRawV(fmt,va);
}

void SiteOutputController::logWithCommunity(const char*fmt, ...)
{
  va_list va;
  logTag();
  va_start(va,fmt);
  site->node->outputcontroller->logRawV(fmt,va);
  va_end(va);
  if ( parameters.outputCommunities() )
  { va_start(va,fmt);
    communities << vfstring(fmt,va);
    va_end(va);
  }
}

void SiteOutputController::recordValues(const VectorAccess<double>*x, double t)
{
  step(t);
}

void SiteOutputController::step(double t)
{
  possiblyReopenFiles(t);
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->step(t);
}
  
void SiteOutputController::installController(OutputController*d)
{
  controllers.push_back(d);
}

void SiteOutputController::recordCommunity(void)
{
  if ( site->community &&
       parameters.reallyWrite() && parameters.outputCommunities() )
  {
    //site->community->printForMathematica(communities);
    //site->community->recordCommunity(communities);
    logCurrentState(site->node->outputcontroller->logfile);
  }
}

void SiteOutputController::speciation(double t,
				      const Index& parent,
				      const Index& daughter)
{
  if ( parameters.outputSpeciations() )
  { log("%g speciation (%s", t, basename(parent).c_str());
    // repeated calls to basename have to be handled separately
    logRaw(" => %s)\n", basename(daughter).c_str());
  }
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->speciation(t, parent, daughter);
}

void SiteOutputController::extinction(double t, const Index& former)
{
  if ( parameters.outputExtinctions() )
    log("%g extinction (%s)\n", t, basename(former).c_str());
  logCurrentState();
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->extinction(t, former);
}

void SiteOutputController::immigration(double t, const Index& newcomer)
{
  log("%g immigration (%s)\n", t, basename(newcomer).c_str());
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->immigration(t, newcomer);
  recordCommunity();
}

void SiteOutputController::explosion(double t, const Index &victim)
{
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->explosion(t, victim);
  recordCommunity();
}

void SiteOutputController::equilibrium(double t)
{
  if ( parameters.outputCommunities() )
    communities << t << " equilibrium:\n";
  log("%g equilibrium\n",t);
  logEquilibrium();
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->equilibrium(t);
}

void SiteOutputController::tiredOfWaitingForEquilibrium(double t)
{
  logWithCommunity("%g tired of waiting for equilibrium...\n", t); 
  logEquilibrium();
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->tiredOfWaitingForEquilibrium(t);
}

void SiteOutputController::ess(double t)
{
  if (site->integrator)
    logWithCommunity("%g %d ess\n", t,
		     site->integrator->evolutionaryTime());
  else
    logWithCommunity("%g ess\n", t);
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->ess(t);
}

void SiteOutputController::logEquilibrium(void)
{
  recordCommunity();
  //  logCurrentState(site->node->outputcontroller->logfile);
  logCurrentState();
  flush();
}

void SiteOutputController::logCurrentState()
{
  logCurrentState(site->node->outputcontroller->logfile);
}

void SiteOutputController::logCurrentState(ostream&os)
{
  os << "{ ";
  if ( site->integrator )
  {
    VectorAccess<double> &x = site->integrator->state();
    /*
    os << "species -> {";
    { bool fi = true;
      for ( int i = 0; i < (int)x.size(); i++ )
	if ( site->community->isVariableInUse(i) &&
	     site->community->isVariableAPopulation(i) )
	  os << (fi? ((fi=false),""):", ") 
	    //<< site->community->variableIndexing().key(i);
	     << basename(Index(i,site->integrator->variableIndexing());
    }
    os << "},\n  ";
    */
    os << "t -> " << site->integrator->time();
    os << ",\n  t_ev -> " << site->integrator->evolutionaryTime();
    for ( int i = 0; i != (int)x.size(); ++i )
      if ( site->community->isVariableInUse(i) )
	os << ",\n  "
	   << basename(Index(i,site->integrator->variableIndexing()))
	   << " -> " << x[i];
  }
  site->community->addToCurrentState(os);
  os << "\n}\n";
  //os << "====================\n";
  os.flush();
}

void SiteOutputController::flush(void)
{
  if (site->node->outputcontroller != this)
    site->node->outputcontroller->flush();
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->flush();
  if ( parameters.reallyWrite() )
  {
    if (logfile.is_open())
      logfile.flush();
    if ( parameters.outputCommunities() )
      communities.flush();
  }
}

void SiteOutputController::finish(void)
{
  for ( int i = 0; i < (int)controllers.size(); ++i )
    controllers[i]->finish();
  //    logfile.close();
  if ( parameters.outputCommunities() )
    communities.close();
}
