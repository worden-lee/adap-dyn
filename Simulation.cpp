#include "Simulation.h"
#include "Parameters.h"
#include "Node.h"
#include "rand.h"
#include <stdlib.h>
#include <iostream>
#include <signal.h>
#include <sys/resource.h>

 // handler for Ctrl-C - end of this file
void handler_terminate(int sig_num);
 // handler for Ctrl-\ - end of this file
void handler_dumpncont(int sig_num);
 // handler for seg violation - end of ...
void handler_abort(int sig_num);    

void Simulation::setup()
{
  /*--------- first do mpi setup stuff -------------------*/
  int mpi_rank, mpi_size;                 /* of MPI_COMM_WORLD */
  
#ifdef MPI
  MPI_Init( &argc, &argv );
  MPI_Comm_rank( MPI_COMM_WORLD, &mpi_rank );
  MPI_Comm_size( MPI_COMM_WORLD, &mpi_size );
  start_clock = MPI_Wtime();
#else
  mpi_rank = 0;
  mpi_size = parameters.rowLength();
  start_clock = time(0);
#endif

  parameters.rankAndSizeAre(mpi_rank, mpi_size);

  setvbuf(stdout, NULL, _IONBF, 0);

  // seed the random numbers (before init'ing community, etc)
  long rand_seed = parameters.randSeed() + 2*mpi_rank;
  sgenrand2(rand_seed);

  // now create the node, sites, and all associated objects
  createNode();
  node->finishInitialize();

  parameters.finishInitialize();

  if ( parameters.outputODEMessages() )
  {
    node->outputcontroller->log("rand_seed %li\n", rand_seed);
  }
  
  // want core file!
  { struct rlimit rl = { RLIM_INFINITY, RLIM_INFINITY };
    if (setrlimit(RLIMIT_CORE,&rl))
      perror("error in setrlimit");
  }
  
#ifdef __sun__
#define SIGNAL sigset
#else
#define SIGNAL signal
#endif //__sun__
  SIGNAL(SIGQUIT, handler_dumpncont); // ^\ to dump out and continue
  SIGNAL(SIGINT, handler_terminate);  // ^C to dump out and quit
  SIGNAL(SIGPIPE, handler_terminate); // if pipe broken, dump + quit
  SIGNAL(SIGSEGV, handler_abort);

  _finished = false;
}

bool Simulation::finished(void)
{
  return _finished;
}

void Simulation::doSimulation(double endtime)
{
  double t = 0;

  _finished = true;
  for ( int i = 0; i < node->nSites; i++ )
  { if ( node->sites[i]->integrator )
      node->sites[i]->integrator->
	setState(t,&node->sites[i]->integrator->state());
    if ( ! node->sites[i]->community->allDead() )
      _finished = false;
  }
  
  // now integrate all sites for a long time
  while ( (t < endtime || endtime == -1) && !_finished )
  {
    double t1 = t + parameters.diffusionTimeStep();

    bool alldead = true;
    bool exploded = false;
    bool finalEquilibrium = true;
    if (parameters.doSpeciation())
      finalEquilibrium = false;

     // integrate each site one beat,
     //  then do synchronous diffusion between sites
    for ( int i = 0; i < node->nSites; i++ )
    {
#ifdef MPI
      ((NodeCommunicator*)(node->communicator))->handleRequests();
#endif //MPI
      // this should be fixed
      node->sites[i]->integrator->integrate(t1);
      if ( ! node->sites[i]->community->allDead() )
	alldead = false;
      if ( node->sites[i]->integrator->hadExplosion() )
	exploded = true;
      else if ( finalEquilibrium &&
		!node->sites[i]->integrator->atEquilibrium() )
	finalEquilibrium = false;
    }
    _finished = (alldead || exploded || finalEquilibrium);
    if (_finished)
      break;
    
    if (node->nSites > 1)
      t = t1;
    else 
      t = node->sites[0]->integrator->time();

    if (node->nSites == 1 &&
	node->sites[0]->integrator->currentlyAtESS())
    {
      //node->sites[0]->outputcontroller->logEquilibrium();
      reachedESS();
      if ( parameters.quitAfterReachingESS() )
	break;
    }
    
    if ( parameters.outputMPIMessages() )
      node->outputcontroller->log("t = %g\n", t);

    for ( int i = 0; i < node->nSites; i++ )
    {
      node->sites[i]->community->possiblyReindex();
      node->sites[i]->communicator->doDiffusion();
    }
  } // end while (t < endtime...)
}

void Simulation::reachedESS(void)
{
}

void Simulation::finish(void)
{
#ifdef MPI
  end_clock = MPI_Wtime();
#else
  end_clock = time(0);
#endif

  for ( int i = 0; i < node->nSites; i++ )
  { node->sites[i]->outputcontroller->recordCommunity();
    node->sites[i]->outputcontroller->logCurrentState();
  }
  
  if ( parameters.outputTiming() )
    node->outputcontroller->log("running time = %g sec\n",
				end_clock - start_clock);
    
  for ( int i = 0; i < node->nSites; i++ )
    node->sites[i]->outputcontroller->finish();

#ifdef MPI
  MPI_Finalize();
#endif //MPI
}

Simulation::~Simulation()
{ //cout << "~Simulation()" << endl;
  delete node;
}

void handler_dumpncont(int sig_num)
{
  simulation.node->outputcontroller->log("interrupted!\n");
  // go down anyway, but record the community structure first
  for ( int i = 0; i < simulation.node->nSites; i++ )
  {
    simulation.node->sites[i]->outputcontroller->log("interrupted!\n");
    simulation.node->sites[i]->outputcontroller->recordCommunity();
    if ( simulation.node->sites[i]->integrator &&
	 parameters.outputMPIMessages() )
      cout << "[" << i << "] t = "
	   << simulation.node->sites[i]->integrator->time() << "\n";
  }
  //cout << "interrupted!\n";
  cout.flush();
}

void handler_terminate(int sig_num)
{
  handler_dumpncont(sig_num);
  if ( !simulation.node )
  {
    cerr << "Error - global node pointer is null" << endl;
  }
  else
    for ( int i = 0; i < simulation.node->nSites; i++ )
    {
      simulation.node->sites[i]->outputcontroller->finish();
    }
#ifdef MPI
  MPI_Finalize();
#endif //MPI
  exit(-1);
}

void handler_abort(int sig_num)
{
  cerr << "Received SIGSEGV or other abort(), attempting to dump core..."
       << endl;
  abort();
}
