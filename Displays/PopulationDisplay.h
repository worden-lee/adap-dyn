/* -*- C++ -*- */
#ifndef POPULATIONDISPLAY_H
#define POPULATIONDISPLAY_H

#include "TimeSeriesDisplay.h"
#include "Parameters.h"
#include "Indexing.h"
#include "util.h"
#include <fstream>
#include <vector>
#include <list>
#include <map>

#include "OutputController.h"

class PopulationController;

typedef _Indexing_base::unique_key_type pop_key_t;
class PopulationDisplay :
  public TimeSeriesDisplay<pop_key_t>
{
protected:
  int nextcolor;
public:
  PopulationDisplay(PopulationController*con);
  virtual ~PopulationDisplay() { }
  void initialize(void);
  string title(timeseries_t &ts);
  virtual void writeGnuplotHeaders(void);
};

class PopulationController :
  public TimeSeriesController<pop_key_t>
{
protected:
  double wpart;
  double dfactor, rfactor;
public:
  PopulationController(Site *s, double period,
		       double recPeriod, double wind);
  virtual ~PopulationController()
  { if (display) delete display; }

  void createDisplay();

  virtual string title(void)
  { return "Micro-nature"; } 
  virtual VectorAccess<double> &x(void)
  { return site->integrator->state(); };
//   virtual double time(void)
//   { return site->integrator->time(); };
  virtual string basename(pop_key_t k)
  { return site->outputcontroller->basename(Index(k)); }
  virtual string description(pop_key_t k)
  { return ""; }
  virtual string filenamebase(void)
  { return "population";
  }
  virtual int color(pop_key_t k)
  { return 1;  }
//   virtual bool trackLineages(void)
//   { return false; } 
//   virtual int lineage(int i)
//   { return 0; }
//   virtual int key(int i)
//   { return indexing().key(i); }
  virtual Indexing &indexing(void)
  { return site->integrator->variableIndexing(); } 
  virtual double displayThreshold(void)
  { return parameters.extinctionThreshold(); } 
  // whether to use y2axis at all
  virtual bool axis2(void)
  { return false; } 
  // if so, whether to use it for a given variable
  virtual bool axis2(pop_key_t k)
  { return false; }

  void updateDisplay(void);
  void newSeries(pop_key_t k);
  void recordCommunityWithout(const Index &skip);
  void recordCommunity(void);
  void extinction(double t, const Index &decedent);
  void speciation(double t, const Index &from, const Index &to);
  void explosion(double t);
};

#endif // POPULATIONDISPLAY_H
