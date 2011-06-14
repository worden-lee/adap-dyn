/* -*- C++ -*- */
#ifndef DISPLAY_H
#define DISPLAY_H
#include "Integrator.h"
#include <math.h>

class OutputController;

class Display
{
public:
  virtual ~Display()
  {// cout << "~Display()" << endl;
  }
  virtual void initialize(void) {}
  virtual void flush() {}
  // override one or both of these
  virtual void updateDisplay(void) {}
  virtual void recordCommunity(void) {}
};

#endif // DISPLAY_H
