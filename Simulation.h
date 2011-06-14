/* -*- C++ -*- */
/*
 * responsibilities of application writer:
 * provide a Simulation subclass which provides createNode() and
 * populateSite(), and one object of that class to fill the extern
 * below.
*/
#include "Parameters.h"

#ifndef SIMULATION_H
#define SIMULATION_H

class Node;
class Site;

class Simulation
{
public:
  Node *node;
  double start_clock, end_clock;
  
  virtual void setup(void);
  virtual void doSimulation(double endtime = parameters.runLength());
  virtual void finish(void);

  virtual void reachedESS(void);
  
  virtual bool finished(void);
  
  // these must be provided by subclass
  virtual void createNode() = 0;
  virtual void populateSite(Site *client) = 0;

  virtual ~Simulation();
protected:
  bool _finished;
};

extern Simulation &simulation;

#endif //SIMULATION_H
