#include "DotDisplay.h"

DotGraphController::DotGraphController(double dperiod, double fperiod, Site*s)
  : DisplayController<DotDisplay>(dperiod, fperiod, s)
{}

void DotGraphController::createDisplay()
{ display = new DotDisplay(this); }

void DotGraphController::_update(void)
{
  if (allowRecording())
    recordCommunity();
  if (allowDisplaying())
    DisplayController<DotDisplay>::updateDisplay();
}
void DotGraphController::recordCommunity()
{
  if (allowRecording())
  { if (!display)
      createDisplay();
    display->recordCommunityWithout(Index::nullIndex());
  }
}

void DotGraphController::updateDisplay(void)
{
  if (!equilibriumonly())
    _update();
}

// this update, and the SiteOutputController::logCurrentState(), should
// both happen at all the same times

void DotGraphController::extinction(double t, const Index &decedent)
{
  //site->outputcontroller->log("%g Dot extinction\n", t);
  if ((displayEvery > 0) && (!equilibriumonly()))
  {
    if (!display)
      createDisplay();
    if (allowRecording())
      display->recordCommunityWithout(decedent);
    if (allowDisplaying())
      DisplayController<DotDisplay>::updateDisplay();
  }
}
void DotGraphController::speciation(double t, const Index &arrival)
{
  //site->outputcontroller->log("%g Dot speciation\n", t);
  if ((displayEvery > 0)&&(!equilibriumonly()))
    _update();
}
void DotGraphController::immigration(double t, const Index &arrival)
{
  //site->outputcontroller->log("%g Dot immigration\n", t);
  if ((displayEvery > 0)&&(!equilibriumonly()))
    _update();
}
void DotGraphController::explosion(double t, const Index &victim)
{
  //site->outputcontroller->log("%g Dot explosion\n", t);
  if ((displayEvery > 0)&&(!equilibriumonly()))
    _update();
}
void DotGraphController::equilibrium(double t)
{
  //site->outputcontroller->log("%g Dot equilibrium\n", t);
  if ((displayEvery > 0))//&&(equilibriumonly()))
    _update();
}
void DotGraphController::tiredOfWaitingForEquilibrium(double t)
{
  //site->outputcontroller->log("%g Dot tiredOfWaiting\n", t);
  if ((displayEvery > 0)&&(equilibriumonly()))
    _update();
}

DotGraphController::~DotGraphController()
{ delete display; }
