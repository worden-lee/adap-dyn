#include "Communicator.h"
#include "Community.h"
#include "Integrator.h"
#include "Site.h"
#include "Node.h"
#include "Indexing.h"
#include <stdlib.h>
#include <iostream>
#include <string>

#ifdef MPI
#include <mpi.h>
#endif

              /* at each diffusion time we have to exchange: 
		 ids of all things diffusing;
		 amount of each thing diffusing;
		 and maybe specifications of things
	      */
#define MIN_BLIP_SIZE  65536

inline Communicator *Communicator::communicatorForSite(int nsite)
{
  int idx = nsite - ( parameters.sitesPerProcessor() * parameters.rank() );
  return site->node->sites[idx]->communicator; 
}

Communicator::Communicator()
{
  site = NULL;
}

void Communicator::siteIs(Site *i)
{
  site = i;
}

void Communicator::doDiffusion()
{
  int row = site->row, col = site->col;
  const VectorAccess<double> &xVector = site->integrator->state();
  vector<double> dVector(xVector.size()); // make it the same size

  if (parameters.totalGridSize() == 1)
    return;

  if (parameters.outputMPIMessages())
  {
    char buf[10000];
    sprintf(buf,"%g beginning diffusion, variables = (",
	    site->integrator->time());
    bool comma = false;
    for (int vi = 0; vi != site->community->nVars(); vi++ )
      if (xVector[vi] > 0)
      {
	Indexing::unique_key_type vk
	  = site->community->variableIndexing().key(vi);
	sprintf(strchr(buf,'\0'),
		(comma? ", %i" : ((comma=true)," %i")), vk);
      }
    strcat(buf," )\n");
    site->node->outputcontroller->log("%s",buf);
  }

  stl_dvectorAccess dva(dVector);
  site->community->calcDiffusion(&xVector, &dva);
  // have to use blip for communicating even locally because everything
  // has to be indexed by uniqueIDs
  const char *out_blip = writeDiffusionBlip(&dva);

  // forget about parity, use buffered send then polling for receive
  sendDiffusionToNeighbor(out_blip, row, col-1);
  sendDiffusionToNeighbor(out_blip, row+1, col);
  sendDiffusionToNeighbor(out_blip, row, col+1);
  sendDiffusionToNeighbor(out_blip, row-1, col);

  receiveDiffusionFromNeighbor(row, col-1);
  receiveDiffusionFromNeighbor(row+1, col);
  receiveDiffusionFromNeighbor(row, col+1);
  receiveDiffusionFromNeighbor(row-1, col);

  if ( parameters.outputMPIMessages() )
  {
    char buf[10000];

    // not needed since it's a ref now?
    //xVector = site->integrator->state();
    sprintf(buf, "%g finishing diffusion, variables = (",
	    site->integrator->time());
    bool comma = false;
    for ( long vi = 0; vi != site->community->nVars(); vi++ )
      if ( xVector[vi] > 0 )
      {
	Indexing::unique_key_type vk
	  = site->community->variableIndexing().key(vi);
	sprintf(strchr(buf,'\0'),
		(comma? ", %i" : ((comma=true)," %i")), vk);
      }
    strcat(buf," )\n");
    site->node->outputcontroller->log("%s",buf);
  }
}

void Communicator::sendDiffusionToNeighbor(const char *out_blip,
					   int nrow, int ncol)
{
  int nsi = siteIndexOf(nrow, ncol);
  if ( nsi != -1 )
  {
    int nr = rankOfSite(nsi);
	  // take out the amount diffused away
    // local indexes may have changed
    subtractDiffusion( out_blip, parameters.rank() );
    if ( nr == parameters.rank() )
    {   // if the other site is on the same processor
      if ( parameters.outputMPIMessages() )
      {
	site->node->outputcontroller->log(
		"give diffusion directly to [%i][%i]\n", nr, nsi);
      }
      Communicator *nc = communicatorForSite(nsi);
      nc->addDiffusion(out_blip, nr);
    }
#ifdef MPI
    else
    {
      if (parameters.outputMPIMessages)
      {
	site->node->outputcontroller->log(
          "send diffusion to [%i][%i]\n", nr, nsi);
      }
      sendDiffusionBlip(out_blip, nr);
    }
#endif //MPI
  } // end if not off the edge of the grid
}

void Communicator::receiveDiffusionFromNeighbor(int nrow, int ncol)
{
  int nsi = siteIndexOf(nrow, ncol);
  if ( nsi != -1 )
  {
    int nr = rankOfSite(nsi);
    if ( nr == parameters.rank() )
    {   // if the other site is on the same processor
	// it has already given it to us directly
    }
#ifdef MPI
    else
    {
      if (parameters.outputMPIMessages)
      {
	site->node->outputcontroller->log(
                "receive diffusion from [%i][%i]\n", nr, nsi);
      }
      // will these come in the right order from the sites on node nr?
      const char *in_blip = receiveDiffusionBlip(nr);
      addDiffusion(in_blip, nr);
    }
#endif //MPI
  } // end if not off the edge of the grid
}

/*000 the sparc needs structure alignment so I have to waste a little
 *    space in the structure
 *  Very possible this ifdef should check something else
 *    or I should fix it with a more intelligent structure definition
 */
#ifdef __sun__
#define calcDiffusionBlipSize(nstructs) \
        ((1+nstructs) * sizeof(struct diffusion_blip_struct))
#else
#define calcDiffusionBlipSize(nstructs) \
        (sizeof(long) + (nstructs) * sizeof(struct diffusion_blip_struct))
#endif

const char *Communicator::writeDiffusionBlip(VectorAccess<double> *dva)
{
  static char *buf = 0;
  if ( parameters.outputMPIMessages() )
  {
    if ( buf == 0 )
      buf = new char[1000];
    strcpy(buf,"Species to send in diffusion: (");
  }
  char *blip = (char *)malloc(calcDiffusionBlipSize(dva->size()));
  long *countp = (long *)blip;
  *countp = dva->size();
#ifdef __sun__  // alignment problems
  struct diffusion_blip_struct *dsp = 
    ((struct diffusion_blip_struct *)blip) + 1;
#else
  struct diffusion_blip_struct *dsp = 
    (struct diffusion_blip_struct *)(countp + 1);
#endif
  
  int i, v;
  for (v = i = 0;
       v != (int)dva->size(); ++v)
  {
    Index vi(v,site->integrator->variableIndexing());
    if ( site->community->isVariableInUse(vi) && (*dva)[v] > 0 )
    {
      dsp[i].uid = vi.key();
      dsp[i].value = (*dva)[v];
      
      if ( parameters.outputMPIMessages() )
      {
	sprintf(strchr(buf,'\0'),
		(i == 0 ? " %i": ", %i"), dsp[i].uid);
      }
      i++;
    }
    else
      (*countp)--;
  }
  if ( parameters.outputMPIMessages() )
  {
    strcat(buf," )\n");
    site->node->outputcontroller->log("%s",buf);
  }
  return blip;
}

void Communicator::addOrSubtractDiffusion(const char *blip, int sign,
					  int sent_from)
{
  VectorAccess<double> &state = site->integrator->state();
  long *countp = (long *)blip;
#ifdef __sun__
  struct diffusion_blip_struct *dsp = 
    ((struct diffusion_blip_struct *)blip) + 1;
#else
  struct diffusion_blip_struct *dsp = 
    (struct diffusion_blip_struct *)(countp + 1);
#endif
  for ( long i = 0; i < *countp; i++ )
  {
    int vi = site->community->variableIndexing().index(dsp[i].uid);
    if ( vi != -1 )
      state[vi] += sign * dsp[i].value;
    else
      cerr << "received unknown unique_key " << dsp[i].uid
	   << " -- don't know what to do" << endl;
    
  }
}

void Communicator::subtractDiffusion(const char *blip, int sent_from)
{
  addOrSubtractDiffusion(blip, -1, sent_from);
}

void Communicator::addDiffusion(const char *blip, int sent_from)
{
  addOrSubtractDiffusion(blip, +1, sent_from);
}

#ifdef MPI
void Communicator::sendDiffusionBlip(const char *out_blip, int other_rank)
{
  // any outstanding requests here can deadlock like the grim reaper
  site->node->communicator->handleRequests();

  long *countp = (long *)out_blip;
  unsigned int charcount = calcDiffusionBlipSize(*countp);

  if (parameters.outputMPIMessages)
  {
    site->node->outputcontroller->log( 
             "Sending diffusion to %d, message size = %u\n", 
	     other_rank, charcount);
    site->node->outputcontroller->logfile.flush();
  }

  // have to use buffered send so I can keep checking for data requests
  Communicator::ensure_MPI_send_buffer_size(charcount);
  MPI_Bsend((void*)out_blip, charcount, MPI_UNSIGNED_CHAR,
	    other_rank, DIFFUSION_TAG, MPI_COMM_WORLD);

  //  if (parameters.outputMPIMessages)
  //  {
  //    site->node->outputcontroller->log.form( "[%i] Sent\n", parameters.rank);
  //    site->node->outputcontroller->log.flush();
  //  }
}

const char *Communicator::receiveDiffusionBlip(int other_rank)
{
  static char *in_blip = NULL;
  static size_t in_blip_size = 0;
  MPI_Status  status;

  if (parameters.outputMPIMessages)
  {
    site->node->outputcontroller->log(
             "Probing for diffusion message from %i\n", other_rank);
    site->node->outputcontroller->logfile.flush();
  }

  int flag;
  do
  {
    site->node->communicator->handleRequests();
    //    if (parameters.outputMPIMessages)
    //    {
    //      fputc('.', stdout);
    //      fflush(stdout);
    //    }
    MPI_Iprobe(other_rank, DIFFUSION_TAG, MPI_COMM_WORLD, &flag, &status);
  }  while(!flag);

  int recv_size;
  MPI_Get_count(&status, MPI_CHAR, &recv_size);
  
  if (parameters.outputMPIMessages)
  {
    site->node->outputcontroller->log(
             "Probe says diffusion message size = %d\n", recv_size);
    site->node->outputcontroller->logfile.flush();
  }

  if ( in_blip_size < recv_size*sizeof(char) )
  {
    while ( in_blip_size < recv_size*sizeof(char) )
      in_blip_size += MIN_BLIP_SIZE;
    in_blip = (char*)realloc(in_blip,in_blip_size);
  }
  
  if (parameters.outputMPIMessages)
  {   
    site->node->outputcontroller->log( "Receiving diffusion from %d\n",
	     other_rank);
    site->node->outputcontroller->logfile.flush();
  }

  MPI_Recv(in_blip, recv_size, MPI_CHAR,
	   other_rank, DIFFUSION_TAG, MPI_COMM_WORLD, &status);

  //  if (parameters.outputMPIMessages)
  //  {
  //    site->node->outputcontroller->log.form( "[%i] Received diffusion\n", parameters.rank);
  //    site->node->outputcontroller->log.flush();
  //  }
  return in_blip;
}

void Communicator::ensure_MPI_send_buffer_size(int atLeast)
{
  static char *mpi_send_buf = NULL;
  static int mpi_send_buf_size = 0;
  atLeast += MPI_BSEND_OVERHEAD;
  if ( mpi_send_buf_size < atLeast )
  {
    if (mpi_send_buf_size > 0)
      MPI_Buffer_detach(&mpi_send_buf, &mpi_send_buf_size);
    mpi_send_buf = 
      (char *)realloc(mpi_send_buf, (mpi_send_buf_size = atLeast));
    if ( mpi_send_buf == 0 )
      (cerr << "Error reallocating mpi_send_buf!\n").flush();
    //    if (parameters.outputMPIMessages)
    //    {
    //      cout.form(
    //               "[%i] Attaching MPI_Bsend buffer, size = %i\n",
    //	       parameters.rank, mpi_send_buf_size);
    //      cout.flush();
    //    }
    MPI_Buffer_attach(mpi_send_buf, mpi_send_buf_size);
  }
}

#endif //MPI
