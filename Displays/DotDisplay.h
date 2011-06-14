/* -*- C++ -*- */
#ifndef DOTDISPLAY_H
#define DOTDISPLAY_H

#include "Display.h"
#include "DisplayController.h"
#include <sys/types.h> // mkdir
#include <sys/stat.h>  // mkdir
#include <fstream> // ofstream
#include <signal.h> // kill
#include <sys/wait.h> // waitpid

class DotDisplay;

class DotGraphController : public DisplayController<DotDisplay>
{
public:
  DotGraphController(double period, double fperiod, Site *s);
  virtual ~DotGraphController();
  
  void createDisplay();
  
  virtual double time(void) = 0;
  virtual bool approve(void) 
  { return true; }
  virtual bool allowDisplaying(void)
  { return approve(); } 
  virtual bool allowRecording(void)
  { return true; }

  virtual vector<string> graphattrs(void)
  { return vector<string>();
  }
  virtual int setsofedges(void) = 0;
  virtual int nedges(int set) = 0;
  virtual bool doedge(int set, int i, int skip) = 0;
  virtual bool donodetitles(void)
  { return false; } 
  virtual string nodelabel(int i) 
  { return ""; }
  virtual string nodetitle(int i) 
  { return ""; }
  virtual string nodecolor(int i) 
  { return ""; }
  virtual string edgetitle(int set, int i) 
  { return ""; }
  virtual string edgecolor(int set, int i) 
  { return ""; }
  virtual bool doedgebackward(int set, int i)
  { return false; } 
  virtual bool boldedge(int set, int i) 
  { return false; }
  virtual bool dottededge(int set, int i) 
  { return false; }
  virtual bool weightedge(int set, int i) // for placement
  { return false; }
  virtual vector<string> otheredgeattrs(int set, int i)
  { return vector<string>();
  }
  virtual int from(int set, int i) = 0;
  virtual int to(int set, int i) = 0;
  virtual int nnodes(void) // no sets of nodes
  { return 0; } 
  virtual int rootnode(void)
  { return -1; } // this means no root
  virtual bool donode(int i, int skip=-1) 
  { return true; }
  virtual bool boldnode(int i) 
  { return false; }
  virtual int nodeindex(const Index &ix) = 0; // convert to node #
  virtual int edgeindex(const Index &ix, int set) = 0; // convert to edge #
  virtual bool equilibriumonly(void)
  { return false; } 
  virtual string title(void)
  { return "Community Graph"; }
  virtual string dotdir()
  { return outdir()+"/dot"; }
  virtual string graphfile(void)
  { return (dotdir()+"/graph.dot").c_str(); }
  virtual string epsfile(void)
  { return (dotdir()+"/graph.eps").c_str(); }
  
  void _update(void);
  void recordCommunity();
  void updateDisplay(void);
  void extinction(double t, const Index &decedent);
  void speciation(double t, const Index &arrival);
  void immigration(double t, const Index &arrival);
  void explosion(double t, const Index &victim);
  void equilibrium(double t);
  void tiredOfWaitingForEquilibrium(double t);
};

class DotDisplay: public Display
{
protected:
  DotGraphController *xs;
  ofstream os;
  char graphfile[1024];
  //  double every,counter;
  bool startedgv;
  int gvpid;
  bool wrote;
  
public:
  DotDisplay(DotGraphController *a)
    : xs(a), startedgv(false), gvpid(-1), wrote(false)
  {
    (void)mkdir(dotdir().c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
    //if(displayEvery>0)
      updateDisplay();
  }

  void recordCommunityWithout(const Index &skip)
  {
    //double time = xs->time();
    //os.seekp(0);
    //cerr<<"without "<<(long)wo << endl;
    sprintf(graphfile,"%s/%s",dotdir().c_str(),xs->graphfile().c_str());
    os.open(graphfile);
    bool bracket=false;
    os << "digraph \"" << xs->title() << "\" {\n";

    vector<string> gas = xs->graphattrs();
    for (vector<string>::iterator it = gas.begin();
	 it != gas.end(); ++it)
      os << (bracket?",":(bracket=true,"graph ["))
	 << *it;
    if (bracket)
      os << "];\n";
    
    bool dnt = xs->donodetitles();
    int rn = xs->rootnode();
    for ( long i = 0; i < xs->nnodes(); ++i )
      if (xs->donode(i,xs->nodeindex(skip)))
      {
	bracket = false;
	bool bn = xs->boldnode(i);
	string nc = xs->nodecolor(i);
	if (dnt || bn || nc.length() || i==rn)
	{
	  os << "  " << xs->nodelabel(i);
	  if (dnt)
	    os << (bracket?",":(bracket=true," ["))
	       << "label=\"" << xs->nodetitle(i) << "\"";		  
	  if (bn)
	    os << (bracket?",":(bracket=true," ["))
	       << "style=bold";
	  if (nc.length())
	    os << (bracket?",":(bracket=true," ["))
	       << "color=" << nc;
	  os << (bracket?"]":"") << ";\n";
	}
      }
    for ( long set = 0; set < xs->setsofedges(); set++ )
      for ( long i = 0; i < xs->nedges(set); i++ )
	if (xs->doedge(set,i,xs->edgeindex(skip,set)))
	{
	  bracket = false;
	  if (!xs->doedgebackward(set,i))
	  {
	    // have to be done in separate commands
	    os << "  " << xs->nodelabel(xs->from(set,i));
	    os << " -> " << xs->nodelabel(xs->to(set,i));
	  }
	  else
	  {
	    os << "  " << xs->nodelabel(xs->to(set,i));
	    os << " -> " << xs->nodelabel(xs->from(set,i));
	    os << " [dir=back";
	    bracket=true;
	  }
	  string el=xs->edgetitle(set,i);
	  if (el.length())
	    os << (bracket?",":(bracket=true," ["))
	       << "label=\"" << el << "\"";
	  string ec=xs->edgecolor(set,i);
	  if (ec.length())
	    if (ec[0]!='\0')
	      os << (bracket?",":(bracket=true," ["))
		 << "color=" << ec;
	  if (xs->dottededge(set,i))
	  {
	    if (xs->boldedge(set,i))
	      os << (bracket?",":(bracket=true," ["))
		 << "style=\"dotted,bold\"";
	    else
	      os << (bracket?",":(bracket=true," ["))
		 << "style=dotted";
	  }
	  else if (xs->boldedge(set,i))
	    os << (bracket?",":(bracket=true," ["))
	       << "style=bold";
	  if (xs->weightedge(set,i))
	    os << (bracket?",":(bracket=true," ["))
	       << "weight=2";
	  vector<string> oed = xs->otheredgeattrs(set,i);
	  for (vector<string>::iterator edi = oed.begin();
	       edi != oed.end(); ++edi)
	    os << (bracket?",":(bracket=true," ["))
	       << *edi;
	  os << (bracket?"]":"") << ";\n";
	}
    os << "}\n";
    os.close();
    wrote = true;
  }
  string dotdir(void)
  { return xs->dotdir(); }
  void ensureGV(void)
  { // start/restart child ghostview if needed
    if (parameters.disableDisplaying())
      return;

    // gvpid is -1 if we haven't forked yet
    // if we have, make sure it's still running
    if (gvpid > 0)
    {
      int checkpid = waitpid(gvpid,0,WNOHANG);
      switch (checkpid)
      {
      case  0: // it's still running
	break;
      case -1: // error
	cerr << "error in DotDisplay::ensureGV() : "
	     << strerror(errno) << endl;
	break;
      default: // if return value is gvpid, it died
	if (checkpid != gvpid)
	  cerr << "DotDisplay: why does waitpid return "
	       << checkpid << " instead of " << gvpid
	       << "?" << endl;
// 	cerr << "DotDisplay " << this << ": waitpid(" << gvpid << ")"
// 	     << " returned nonzero : restarting gv" << endl;
	gvpid = -1; // need to start a new one
	break;
      }
    }

    if (gvpid <= 0) // if we need to start a gv process
    {
      int pid = fork();
      if (pid == 0)
      { // child process in here
	const string &gvcomm = parameters.ghostviewCommand();
	string eps = dotdir() + "/" + xs->epsfile();
	//cout << "DotDisplay child process: execlp "
	//     << gvcomm << ' ' << buf << endl;
	execlp(gvcomm.c_str(), gvcomm.c_str(), eps.c_str(), NULL);
	// if we get here there's a problem
	cerr<<"DotDisplay had error starting "
	    << gvcomm << ' ' << eps << endl;
	exit(-1);
      }
      // parent process out here
      if (pid < 0)
	cerr << "DotDisplay" << this
	     << ": can't fork gv subprocess" << endl;
      else
      { gvpid = pid;
//         cerr << "DotDisplay " << this
// 	     << ": gv subprocess is " << gvpid << endl;
      }
    }
    else // if gv is already running
    { // tell it to refresh itself
//       cerr << "DotDisplay " << this
// 	   << ": send HUP to " << gvpid << endl;
      kill(gvpid,SIGHUP);
    }
  }

  void updateDisplay(void)
  {
    if (wrote)
    {
      //double time = xs->time();
      string dotcomm =
	string("/usr/bin/dot -Tps -o") + dotdir() + "/" + xs->epsfile()
	+ " " + graphfile;
      
      //cout << dotcomm << endl;
      if(system(dotcomm.c_str()))
	cerr << "DotDisplay had problem running " << dotcomm << endl;

      ensureGV();
    } 
  }

  ~DotDisplay()
  {
    if (gvpid > 0)
    {
//       cerr << "DotDisplay::~DotDisplay() : kill child process "
// 	   << gvpid << endl;
      kill(gvpid,SIGINT);
      //cout << "Waiting for gv to terminate" << endl;
      waitpid(gvpid,0,0);
//       cerr << "DotDisplay::~DotDisplay() : "
// 	   << gvpid << " terminated" << endl;
    }
  }
};

#endif // DOTDISPLAY_H
