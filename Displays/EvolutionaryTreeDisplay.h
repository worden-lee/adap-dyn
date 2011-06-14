/* -*- C++ -*- */
#ifndef EVOLUTIONARYTREEDISPLAY_H
#define EVOLUTIONARYTREEDISPLAY_H
#include "TimeSeriesDisplay.h"
#include <fstream>
#include <iterator>

class EvolutionaryTreeController;

typedef int ev_key_t;

class EvolutionaryTreeDisplay: public TimeSeriesDisplay<ev_key_t>
{
protected:
  map< key_t,set<key_t> > descent;
  
public:
  EvolutionaryTreeDisplay(EvolutionaryTreeController *a);
  virtual ~EvolutionaryTreeDisplay() {}

  void initialize(void)
  { TimeSeriesDisplay<ev_key_t>::initialize();
  }

  virtual void descendsFrom(key_t daughter, key_t parent,
			    double birthTime)
  { descent[parent].insert(daughter);
    connectLineages(daughter,parent,birthTime); 
  }
  virtual void deactivate(key_t deceased, double t)
  {
    series_map_t::iterator deci = series.find(deceased);
    if (deci==series.end())
      return;
    // if parent branch expires at same time an offspring branch
    // is born, merge the two
    timeseries_t &dects = deci->second;
    set<key_t> &offs = descent[deceased];
    for (set<key_t>::iterator offi = offs.begin();
	 offi != offs.end(); ++offi)
    {
      series_map_t::iterator offmi = series.find(*offi);
      // this doesn't work if the thing goes extinct before
      // ever being plotted
      if (offmi == series.end())
	continue;
      timeseries_t &offts = offmi->second;
      if (offts.series[0].t == t)
      {
// 	cout << "series join " << deceased << " to " << *offi << endl;
// 	cout << "before:\nparent ";
// 	transform(dects.series.begin(),dects.series.end(),
// 		  ostream_iterator<double>(cout, " "),&get_t);
// 	cout << "\noffspring ";
// 	transform(offts.series.begin(),offts.series.end(),
// 		  ostream_iterator<double>(cout, " "),&get_t);
// 	cout << "\n";
	offts.series.insert(offts.series.begin(),
			    dects.series.begin(),dects.series.end()-1);
	offts.stoptimes.insert(offts.stoptimes.begin(),
			       dects.stoptimes.begin(),
			       dects.stoptimes.end());

	removeKeyFromMemory(deci->first);
	series.erase(deci);

// 	cout << "after: ";
// 	transform(offts.series.begin(),offts.series.end(),
// 		  ostream_iterator<double>(cout, " "),&get_t);
// 	cout << endl;
	return;
      }
    }
      //else    
    TimeSeriesDisplay<ev_key_t>::deactivate(deceased,t);
    descent.erase(descent.find(deceased));
  }
  void switchToDots(bool d)
  {
    if (!lines)
    {
      if (d)
	gnuplot << "set pointsize 0.1\n";
      else
	gnuplot << "set pointsize 2\n";
    }
  }

  virtual void writeGnuplotHeaders(void) {
    TimeSeriesDisplay<ev_key_t>::writeGnuplotHeaders();
    gnuplot << "set parametric\n";
    if (!lines)
      gnuplot << "set pointsize 2\n";
    gnuplot << "set nokey\n";
  }
  
};

// this object controls an EvolutionaryTreeDisplay 
//  and provides the actual phenotypic data.
// Whoever wishes to use the display object must supply
//  a subclass of this class that provides meaningful
//  instantiations of these functions.
class EvolutionaryTreeController
  : public TimeSeriesController<ev_key_t>
{
public:
  // for generality, there are up to count() counter_ts
  // and each gets converted to a key_t (if alive)
  // indices can be converted to counters
  typedef ev_key_t key_t;
  typedef int counter_t;

  EvolutionaryTreeController(double period, double filePeriod,
			     Site *s)
    : TimeSeriesController<ev_key_t>(period, filePeriod, s)
  { displayCounter = recordCounter = 0;  // don't record at t=0
  }
  virtual ~EvolutionaryTreeController()
  { if (display) delete display;
  }
  
  void createDisplay()
  { display = new EvolutionaryTreeDisplay(this);
    display->initialize();
  }

  timescale_t timescale()
  { return evolutionary_time;
  }
  
  // tired of vertical plots!!
//   virtual bool vertical()
//   { return true; }
  virtual counter_t count(void) = 0;
  virtual bool alive(counter_t i) = 0;
  virtual double s(key_t k) = 0;
  virtual string title(void)
  { return "Evolutionary Tree"; }
  virtual string filenamebase(void)
  { return "evolutionary-tree"; }
//   virtual const char *basename(int i) = 0;
  virtual bool trackLineages(void)
  { return false; } 
//   virtual int lineage(int i)
//   { return 0; }
  virtual int color(key_t i)
  { return 1; }
  virtual counter_t counter(const Index &i) = 0;
  virtual key_t key(counter_t i) = 0;

  void recordCommunity(void)
  { // just shove them at the end of the branches
    if (!display)
      createDisplay();
    double t = time();
    for ( counter_t i = 0; i < count(); i++ )
      if (alive(i))
      {
	key_t k = key(i);
	if (k >= 0)
	  display->record(k,t,s(k));
      }
  }
  void updateDisplay()
  {
    if (!display)
      createDisplay();
    if (time()/recordEvery >= 100)
      display->setLines(false);
    TimeSeriesController<ev_key_t>::updateDisplay();
  }
  void speciation(double t, const Index &parent,
		  const Index &daughter)
  { key_t dk = key(counter(daughter)), pk = key(counter(parent));
    if (dk >= 0 && pk >= 0 && display)
      display->descendsFrom(dk,pk,time());
  }
  void extinction(double t, const Index &victim)
  {
    key_t k = key(counter(victim));
    if (k >= 0 && display)
      display->deactivate(k,time());
  }
  void explosion(double t, const Index &victim)
  {
    if (displayEvery > 0)
      update(site->integrator->evolutionaryTime());
    TimeSeriesController<ev_key_t>::explosion(t,victim);
  }
};

#endif // EVOLUTIONARYTREEDISPLAY_H
