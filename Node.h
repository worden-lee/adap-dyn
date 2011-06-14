/* -*- C++ -*-
 *
 * There should be only one node object per process!
 */

#ifndef NODE_H
#define NODE_H

/* should we make Node a kind of Site with communicator, community, etc? */
#include "Site.h"
#include "Simulation.h"

class Node : public Site
{
public:
  Simulation *simulation;
  // these have to be pointers because a subclass of Simulation
  // will populate them with a subclass of Site
  vector<Site*> sites;
  int nSites;
  Node() 
  { node = this; }
  void assignSites(vector<Site*>&ss)
  {
    sites = ss;
    nSites = sites.size();
    for ( int i = 0; i < nSites; i++ )
    {
      simulation->populateSite(sites[i]);
      sites[i]->nodeIs(this);
    }
  }
  virtual void finishInitialize(void)
  {
      // set up initial community, initial conditions,
      //  get unique IDs of global resources and species right
    Site::finishInitialize();
    for ( int i = 0; i < nSites; i++ )
    {
      sites[i]->siteIndexIs(parameters.rank() * nSites + i);
      sites[i]->finishInitialize();
    }
  }
  virtual bool isANode(void)
  { return true; }
  virtual ~Node() 
  { for ( int i = 0; i < nSites; i++ )
      delete sites[i];
  } 
};

#endif //NODE_H
