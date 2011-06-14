/* -*- C++ -*- */
#ifndef INTERACTIONLINEAGEDISPLAY_H
#define INTERACTIONLINEAGEDISPLAY_H

#include "TimeSeriesDisplay.h"
#include "Integrator.h"
#include <utility> // pair

template <class key_t>
class InteractionLineageController;

template <class K>
class InteractionLineageDisplay :
  public TimeSeriesDisplay< pair<K,K> >
{
public:
  typedef pair<K,K> key_t;

  InteractionLineageDisplay(InteractionLineageController<K>*con)
    : TimeSeriesDisplay< pair<K,K> >(con)
  {} 

  void deactivate_all(K i, double t)
  { typedef typename
      TimeSeriesDisplay< pair<K,K> >::series_map_t::iterator
      it_t;
    for (it_t it = TimeSeriesDisplay<key_t>::series.begin();
	 it != TimeSeriesDisplay<key_t>::series.end();
	 ++it)
      if (it->first.first == i || it->first.second == i)
	TimeSeriesDisplay<key_t>::deactivate(it->first, t);
  }
    
  void descendsFrom(key_t daughter, key_t parent, double birthTime)
  { TimeSeriesDisplay<key_t>::connectLineages(daughter, parent, birthTime);
  }

//   // can't cast a pair to an int, here's a unique diagonal mapping
//   //  (at least if you can cast K to int)
//   lineage_t defaultLineage(key_t k)
//   { K sum = k.first + k.second;
//     return (sum * (sum + 1))/2 + k.second;
//   }
};

template<class K>
class InteractionLineageController :
  public TimeSeriesController< pair<K,K> >
{
public:
  typedef pair<K,K> key_t;
  using DisplayController< TimeSeriesDisplay<key_t> >::site;
  using DisplayController< TimeSeriesDisplay<key_t> >::display;
  using DisplayController< TimeSeriesDisplay<key_t> >::time;
  using DisplayController< TimeSeriesDisplay<key_t> >::outdir;
  using DisplayController< TimeSeriesDisplay<key_t> >::displayEvery;
        
  InteractionLineageController(double period, double filePeriod,
			       Site *s)
    : TimeSeriesController<key_t>(period, filePeriod, s)
  { // don't record at t=0
    this->displayCounter = this->recordCounter = 0;
  }

  virtual ~InteractionLineageController()
  { if (display) delete display;
  }
  
  void createDisplay()
  { display = new InteractionLineageDisplay<K>(this);
    display->initialize();
  }

  timescale_t timescale()
  { return evolutionary_time;
  }

  // user has to provide these
  virtual int count(void) = 0;
  virtual double interaction(int i, int j) = 0;
  virtual int key(const Index &i) = 0;
  virtual int key(int i) = 0;
  virtual bool alive(int i)
  { return true; } 
  virtual bool alive(int i, int j)
  { return true; } 

  virtual bool vertical()
  { return false; }
  virtual string title(void)
  { return "Evolution of Interactions"; }
  virtual string filenamebase(void)
  { return "interactions"; }
  virtual bool trackLineages(void)
  { return false; } 
  virtual int lineage(int i)
  { return 0; }

  void recordCommunity(void)
  { double t = time();
    if (!display)
      createDisplay();
    for ( long i = 0; i < count(); i++ )
      if (alive(i))
      { int ki = key(i);
	if (ki >= 0)
	  for ( long j = 0; j < count(); j++ )
	  { int kj = key(j);
	    if (kj >= 0)
	    { if (alive(j) && alive(i,j))
		display->record(make_pair(ki,kj),t,
				interaction(i,j));
	      else
		display->deactivate(make_pair(ki,kj),t);
	    }
	  }
      }
  }

  void speciation(double t, const Index &parent,
		  const Index &daughter)
  { int dk = key(daughter), pk = key(parent);
    if (dk >= 0 && pk >= 0)
    {
      display->descendsFrom(make_pair(dk,dk),make_pair(pk,pk),t);
      for ( long i = 0; i < count(); i++ )
	if (alive(i))
	{ int ki = key(i);
 	  if (ki >= 0)
	  {
	    display->descendsFrom(make_pair(ki,dk),
				  make_pair(ki,pk),t);
	    display->descendsFrom(make_pair(dk,ki),
				  make_pair(pk,ki),t);
	  }
	}
    }
  }
//   void extinction(double t, const Index &victim)
//   { int k = key(victim);
//     if (k >= 0)
//       static_cast<InteractionLineageDisplay<K>*>(display)
// 	->deactivate_all(key(victim),time());
//   }
  void explosion(double t, const Index &victim)
  { if (displayEvery > 0)
      update(site->integrator->evolutionaryTime());
      //update(evtime);
    TimeSeriesController<key_t>::explosion(t,victim);
  }

  void newSeries(pair<K,K> k)
  {
    enum { ON_COLOR = -1, OFF_COLOR = 1 }; // black/red
    if (k.first == k.second)
      display->setColor(k, ON_COLOR);
    else
      display->setColor(k, OFF_COLOR);
  }
  
};

#endif // INTERACTIONLINEAGEDISPLAY_H
