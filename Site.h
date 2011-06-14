/* -*- C++ -*- */
#ifndef SITE_H
#define SITE_H

class Node;
#include "Communicator.h"
#include "Integrator.h"
#include "OutputController.h"
#include "Community.h"
#include "Parameters.h"

class Site
{
 public:
  int siteIndex;
  int row, col;
  Node *node;
  Community *community;
  Integrator *integrator;
  Communicator *communicator;
  SiteOutputController *outputcontroller;

  Site() : siteIndex(0), row(0), col(0),
	   node(0), community(0), integrator(0),
	   communicator(0), outputcontroller(0)
  {}
  
  void siteIndexIs(int si)
    {				
      siteIndex = si;
      col = siteIndex % parameters.rowLength();
      row = (siteIndex - col) / parameters.rowLength();
    }
  void nodeIs(Node *n)
    {
      node = n;
    }
  virtual void finishInitialize(void)
    {
      if ( communicator )
        communicator->siteIs(this);
      if ( community )
        community->siteIs(this); 
      if ( integrator )
        integrator->siteIs(this);    
      if ( outputcontroller )
        outputcontroller->siteIs(this);
      if ( community )
        community->initialize();
      if ( outputcontroller )
        outputcontroller->installControllers();
    }
  virtual bool isANode(void)
  { return false; }
  virtual ~Site()
  {
    //cout << "~Site()" << endl;
    if (outputcontroller)
      delete outputcontroller;
    if (community)
      delete community;
    if (integrator)
      delete integrator;
    if (communicator)
      delete communicator;
  }
};

#endif //SITE_H
