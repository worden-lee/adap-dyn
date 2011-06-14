/* -*- C++ -*- */
#ifndef VALUEWRITER_H
#define VALUEWRITER_H

#include <stdio.h>
#include <iostream>
#include <fstream>
#include "Display.h"
#include "Parameters.h"
#include "Integrator.h"
#include "Communicator.h"
#include "Vector.h"
#include <stdarg.h>

#include <string>
#include <vector>
#include <list>
#include <map>

using namespace std;

class Site;
class Community;
class Communicator;

// base class for a thing that gets notifications of 
//  events and decides when and how to do output
// subclasses include the SiteOutputController that manages
//  output for a site, and the little DisplayControllers
//  that manage each gnuplot view.
class OutputController
{
public:
  OutputController() {}
  virtual ~OutputController() {}

  // some of these will need to be overridden so that
  //  something will happen
  virtual void step(double t) {} 
  virtual void speciation(double t, const Index &parent,
			  const Index &daughter) {}
  virtual void extinction(double t, const Index &decedent) {}
  virtual void immigration(double t, const Index &newcomer) {}
  virtual void equilibrium(double t) {}
  virtual void tiredOfWaitingForEquilibrium(double t) {}
  virtual void ess(double t) {}
  virtual void explosion(double t, const Index &victim) {}
  virtual void finish(void) {}
  virtual void flush(void) {}
};

class SiteOutputController
{
public:
  ofstream logfile;
  ofstream communities;
  // we manage a herd of little controllers for different displays
  vector<OutputController*> controllers;
  
  SiteOutputController();
  virtual ~SiteOutputController();

  virtual void siteIs(Site *);

  // check outputTimeStep, either do nothing, or call recordValuesForSure
  //  which can be overridden (or this can)
  virtual void recordValues(const VectorAccess<double>*x, double t);
  virtual void recordCommunity(void);

  // just the var. name. Default is "X@@@" where @@@ is the unique key
    // you need to use this string right away b/c the buffer gets reused
  virtual string basename(Index);

  virtual string outputDirectory(void);

  // these get passed on notifications
  // they have to be created with 'new' because we will delete them.
  void installController(OutputController*);

  // override this one, to set up your project's display controllers
  virtual void installControllers(void)
  { }

  virtual void step(double t);
  virtual void equilibrium(double t);
  virtual void tiredOfWaitingForEquilibrium(double t);
  virtual void ess(double t);
  virtual void speciation(double t, const Index &parent,
			  const Index &daughter);
  virtual void extinction(double t, const Index &decedent);
  virtual void immigration(double t, const Index &newcomer);
  virtual void explosion(double t, const Index &victim);
  virtual void finish(void);
  virtual void flush(void);

  // logging

  virtual void logTag(void);
  virtual void writeTag(ostream &);
    // log with tag
  void log(const char*fmt,...)  __attribute__ ((format (printf, 2,3)));
				// 2, 3 not 1, 2 because of 'this', I guess
    // log without tag
  void logRaw(const char*fmt,...)  __attribute__ ((format (printf, 2,3)));
  virtual void logRawV(const char*fmt,va_list va);

  // write to community and log
  virtual void logWithCommunity(const char*fmt,...)
    __attribute__ ((format (printf, 2, 3)));
  virtual void logEquilibrium(void);
  void logCurrentState(ostream&);
  void logCurrentState();
  
protected:
  string dir;
  Site *site;

  virtual void openMemberFiles(void);
  virtual void reopenMemberFiles(void);
  virtual void possiblyReopenFiles(double t);
  
  virtual void initDir(void);
  virtual string outfilename(const Index&); // full path (relative)
};

class NodeOutputController : public SiteOutputController
{
public:
  void logTag(void);
  virtual void logRawV(const char*fmt,va_list va);
  virtual void logWithAllCommunities(char*fmt,...)
    __attribute__ ((format (printf, 2, 3)));
protected:
  // just uses different filenames (and skips some of the files)
  virtual void reopenMemberFiles(void);
  virtual void initDir(void);

  // subclass might want this
  virtual void recordCommunity(void) {}

  void writeTag(ostream&);

  friend class PopAccess;
};

// Color class: lineages and such have color associated, and
//  this thing helps achieve the same color in dot and gnuplot
class Color
{
protected:
  // we only use the colors that are the same in gnuplot
  //  eps output as in gnuplot screen output.
  static int _gpcolors[];
  static string _dotcolors[];
  // _lastColor - _firstColor + 1 must be the number of entries in each
  // of those arrays there.
  static int _firstColor, _lastColor;
  // _nextIndex is a class-wide counter, for producing sequential colors
  static int _nextIndex;
  // myIndex holds an individual's index
  int myIndex;
public: 
  static int nextIndex(void)
  { int ni = _nextIndex;
    ++_nextIndex;
    if (_nextIndex > _lastColor)
      _nextIndex = _firstColor;
    return ni;
  }
  Color(int ix=-1)
  { if (ix != -1) // mod_correct is like % but guaranteed
      myIndex =   //  to return positive
	_firstColor + mod_correct(ix - _firstColor, _lastColor - _firstColor + 1);
    else // assign color ourself
      myIndex = nextIndex();
  }
  static int name_to_index(string cname)
  { for (int i = _firstColor; i <= _lastColor; ++i)
      if (_dotcolors[i] == cname)
	return i;
    cerr << "unknown color name " << cname << endl;
    return 0;
  }
  Color(string cname) : myIndex(name_to_index(cname))
  {} 
  int gpColor() {
    return _gpcolors[myIndex];
  }
  string dotColor() {
    return _dotcolors[myIndex];
  }
  int index() {
    return myIndex;
  }
  void setIndex(int newIndex) {
    myIndex = newIndex;
  }
  void setNextIndex() {
    // For incrementing indexes on a class-wide basis
    myIndex = nextIndex();
  }
};


class counter {
    protected:
	int myCount;                      // Change to local
	//static string countString;
    public:
	string countString;
	counter() {
		myCount = 0;
	}
	string getString() { return countString; }
	void increment() {
	    myCount++;
	}
	int value() { return myCount; }
	string precisionValue(int precision) {
		string temp = iterate(precision, myCount);
		return temp;
	}
	void calculate(int precision) {
		countString = iterate(precision, myCount);
	}
    private:
	string iterate(int precision, int localCount) {
		if (precision == 0)
			return "";
		else {
			string temp = iterate(precision-1, localCount/10) + lsb(localCount);
			return temp;
		}
	}
	string lsb(int localCount) {
		switch(localCount%10) {
			case 0: return "0";
			case 1: return "1";
			case 2: return "2";
			case 3: return "3";
			case 4: return "4";
			case 5: return "5";
			case 6: return "6";
			case 7: return "7";
			case 8: return "8";
			case 9: return "9";
			default: return "0";
		}
		return "";
	}
};


#endif //VALUEWRITER_H
