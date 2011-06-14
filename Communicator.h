/* -*- C++ -*- */
#ifndef COMMUNICATOR_H
#define COMMUNICATOR_H
/*
 * every site (every community object) gets a communicator
 * which can talk to other communicators whether on same processor or not
 */
#include "Integrator.h"
#include "Parameters.h"
class Site;

class Communicator
{
 public:
  Communicator();
  virtual ~Communicator() {}

  //  static int rank, size;
  //  static void rankAndSizeAre(int rk, int sz);
  virtual void siteIs(Site *);
  
  void doDiffusion(void);

#ifdef MPI
  enum mpi_message_tag { DIFFUSION_TAG, FIRST_UNUSED_TAG };
  virtual void handleRequests(void) {};
#endif //MPI

protected:
  Site *site;

  inline static int siteIndexOf(int row, int col)
    {
      if ( row < 0 || row >= parameters.nRows()
	   || col < 0 || col >= parameters.rowLength() )
	return -1;
      int si = row * parameters.rowLength() + col;
      if ( si < parameters.totalGridSize() )  // in case of incomplete row
	return si;
      else return -1;
    }
  inline static int rankOfSite(int site)
    {
      return 
	(int)(site / parameters.sitesPerProcessor());
    }

  // communicatorForSite is null if site is on some other node
  //   i.e. rankOfSite(site) != our own rank
  virtual Communicator *communicatorForSite(int site);

  /*
   * blip format is:
   *  1. (long) # of structs
   *  2. that many structs
   *
   */
  struct diffusion_blip_struct
  { Indexing::unique_key_type uid; double value; };

  virtual const char *writeDiffusionBlip(VectorAccess<double> *dva);
  char *realloc_blip_maybe(const char *blip, char *far_end);

  //  void doDiffusionWithNeighbor(const char *out_blip, int nrow, int ncol);
  void sendDiffusionToNeighbor(const char *out_blip, int nrow, int ncol);
  void receiveDiffusionFromNeighbor(int nrow, int ncol);

    // these two directly modify integrator state
  void subtractDiffusion(const char *blip, int sent_from);
  void addDiffusion(const char *blip, int sent_from);
  virtual void addOrSubtractDiffusion(const char *blip, int sign,
				      int sent_from);

#ifdef MPI
  void sendDiffusionBlip(const char *blip, int other_rank);
  const char *receiveDiffusionBlip(int other_rank);

  static void ensure_MPI_send_buffer_size(int atLeast);
#endif //MPI
};

#endif //COMMUNICATOR_H
