#include "PopulationDisplay.h"
#include <set>

PopulationController::PopulationController(Site *s, double period,
					   double recPeriod, double wind=HUGE)
  : TimeSeriesController<pop_key_t>(period, recPeriod, s, wind)
{
  recordCounter = displayCounter = -HUGE;
  if (timescale() == ecological_time)
  {
    // points: 0, ip, ip + ip*f, ip + ip*f + ip*f^2,
    //  window = ip + ip*f + ... + ip*f^n-1
    //  want ip*f^n = window / wpart, say ...
    // window / ip = 1 + f + ... + f^n-1 = 1 + (1 - f^n)/(1 - f)
    // window / ip = wpart * f^n
    // f^n = (w/ip)/wpart
    // (1-f) = (1 - f^n)/[w/ip - 1]
    wpart = 100000;//100;//20;
    dfactor = 1 - (1 - (window/period)/wpart)/((window/period) - 1);
    rfactor = 1 - (1 - (window/recPeriod)/wpart)/((window/recPeriod) - 1);
  }
}

void PopulationController::createDisplay()
{ display = new PopulationDisplay(this); 
  display->initialize();
}

void PopulationController::updateDisplay()
{
  TimeSeriesController<pop_key_t>::updateDisplay();
  if (timescale() == ecological_time)
    if (displayEvery < window/wpart)
      displayEvery *= dfactor;
}

// this overrides the TimeSeriesController one
void PopulationController::newSeries(pop_key_t k)
{
  if (!display)
    createDisplay();
  display->setBasename(k,basename(k));
  display->setDescription(k,description(k));
  display->setColor(k,color(k));
  display->setAxis2(k,axis2(k));
}

void PopulationController::recordCommunityWithout(const Index &skip)
{
//   // @@@ this is a temporary hack @@@ if it even works @@@
//   truncate(display->gnuplotLogFile().c_str(),0);

  //DisplayController<PopulationDisplay>::recordCommunity();
  double t = time();
  VectorAccess<double> &x = this->x();
  for ( unsigned i = 0; i < x.size(); i++ )
  {
    Index ix(i,site->integrator->variableIndexing());
    if (ix != skip && site->community->isVariableInUse(ix))
    { pop_key_t k = ix.key();
      display->record(k,t,x[i]);
      //display->setColor(k,color(k));// over+ over
    }
  }
  if (timescale() == ecological_time)
    if (recordEvery < window/wpart)
      recordEvery *= rfactor;
}

void PopulationController::recordCommunity()
{
  if (timescale() == evolutionary_time &&
      !site->integrator->atEquilibrium())
    return;
  recordCommunityWithout(Index::nullIndex());
}

void PopulationController::extinction(double t, const Index &decedent)
{
  if (display)
    display->deactivate(decedent.key(),time());
  if (timescale() == ecological_time)
    recordCommunityWithout(decedent);
}

void PopulationController::speciation(double _t,
				      const Index &from,
				      const Index &to)
{ double t = time();
  pop_key_t fk = from.key();
  pop_key_t tk = to.key();
  if (!display)
    createDisplay();
  if (timescale() == ecological_time)
  {
    VectorAccess<double> &x = this->x();
    int fi = site->integrator->variableIndexing().index(fk);
    int ti = site->integrator->variableIndexing().index(tk);
    display->record(fk,t,x[fi]);
    display->record(tk,t,x[ti]);
  }
  else
  {
  //newSeries(tk);
  //display->descendsFrom(tk, fk, t);
    display->connectLineages(tk, fk, t);
  }
}

void PopulationController::explosion(double t)
{
  recordCommunity();
  updateDisplay();
}

PopulationDisplay::PopulationDisplay(PopulationController *con)
  : TimeSeriesDisplay<pop_key_t>(con),
    nextcolor(1)
{}

void PopulationDisplay::initialize(void)
{
  TimeSeriesDisplay<pop_key_t>::initialize();
  if (static_cast<PopulationController&>(xs).axis2())
    gnuplot << "set y2tics\n";
}

void PopulationDisplay::writeGnuplotHeaders(void) {
  TimeSeriesDisplay<pop_key_t>::writeGnuplotHeaders();
  if (static_cast<PopulationController&>(xs).axis2())
    gnuplot << "set y2tics\n";
}

string PopulationDisplay::title(timeseries_t &ts)
{
  bool useTitle = true;
//   if (!ts.series.empty() && ts.series.back().t == t
//       && ts.series.back().x >= xs.displayThreshold())
//     useTitle = true;
//   if (useTitle && !ts.stoptimes.empty() && ts.stoptimes.back() == t)
//     useTitle = false;
  if (!ts.stoptimes.empty() &&
      ts.stoptimes.back() >= ts.series.back().t)
    useTitle = false;
  if (useTitle)
    return ts.basename
      + (ts.description.length() > 0 ? " " + ts.description : "");
  else
    return "";
}
