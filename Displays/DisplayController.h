/* -*- C++ -*- */
#ifndef DISPLAYCONTROLLER_H
#define DISPLAYCONTROLLER_H
#include "OutputController.h"
#include "GnuplotDisplay.h"
#include <math.h>
#include "CVIntegrator.h"
//#include "CVIntegrator.cpp"
#include "Site.h"

// whether to use the integrator's time, or # of mutations
typedef enum { ecological_time, evolutionary_time } timescale_t;

template<class DISPLAY>
class DisplayController : public OutputController
{
protected:
  DISPLAY *display;
  Site *site;
  //timescale_t timescale;
  double recordEvery, recordCounter;
  double displayEvery, displayCounter;
  //int evtime;
  
public:
  // your subclass needs to assign the site pointer, either
  // in its constructor or thereafter.
  DisplayController(double displayPeriod, double filePeriod, Site *s=0/*,
		    timescale_t ts=ecological_time*/)
    : display(0), site(s),// timescale(ts),
      recordEvery(filePeriod), recordCounter(-HUGE),
      displayEvery(displayPeriod), displayCounter(-HUGE)
    // set ...Counter=0 to skip plotting at time=0
  {}
  virtual ~DisplayController() 
  { //cout << "~DisplayController()" << endl;
    //if (display) delete display;
  }

  // your subclass constructor should provide this, to allocate
  //  and initialize an appropriate Display object,
  //  and assign it to display.
  // this gets called from update() so be warned if you don't use that.
  virtual void createDisplay() = 0;
  
  void period(int period)
  { displayEvery = recordEvery = period; }
//   virtual void setTimescale(timescale_t ts)
//   { timescale = ts; }
  virtual timescale_t timescale()
  {
    return ecological_time;
  }
  
  virtual double time()
  { switch(timescale())
    {
    case evolutionary_time:
      return site->integrator->evolutionaryTime();
    default:
    case ecological_time:
      return site->integrator->time();
    }
  }

  // notifications from SiteOutputController
  // subclasses - override all these at will
  virtual void step(double t)
  { //if (timescale() == ecological_time)
      update(time());
  }
  virtual void speciation(double t, const Index& parent,
			  const Index& daughter)
  {}
  virtual void extinction(double t, const Index& former) {}
  virtual void immigration(double t, const Index& newcomer) {}
  virtual void explosion(double t, const Index& victim) {}
  virtual void equilibrium(double t)
  { if (timescale() == evolutionary_time)
      update(time());
  }
  virtual void tiredOfWaitingForEquilibrium(double t)
  { if (timescale() == evolutionary_time)
      update(time());
  }
  virtual void flush(void)
  { if (display) display->flush();
  }
  virtual void finish(void) {}

  virtual string outdir(void)
  { return site->outputcontroller->outputDirectory();
  }

  virtual void updateDisplay(void)
  { if (!display)
      createDisplay();
    display->updateDisplay();
  }
  virtual void recordCommunity(void)
  { if (!display)
      createDisplay();
    display->recordCommunity();
  }

  virtual void update(double t)
  {
    if (!display)
      createDisplay();
    if (recordEvery > 0 && recordCounter+recordEvery <= t)
    {
      recordCommunity();
      recordCounter = t;
    }
    if (displayEvery > 0 && displayCounter+displayEvery <= t)
    {
      updateDisplay();
      displayCounter = t;
    }
  }
};

#endif // DISPLAYCONTROLLER_H
