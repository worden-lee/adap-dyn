/* -*- C++ -*- */
#ifndef VECTOR2DDISPLAY_H
#define VECTOR2DDISPLAY_H
/* this is obsolete, but I keep it around should I need to update it */
#include "Display.h"
#include "Integrator.h"
#include "vnl/vnl_matrix.h"
#include <fstream>

#define VECTORFILE     "out/vectors.gp"
class VectorAccess2D
{
public:
  virtual int count(void) = 0;
  virtual bool alive(int i) = 0;
  virtual vnl_vector<double> s(int i) = 0;
  virtual vnl_matrix<double> &M(void) = 0;
  virtual double X(int i) = 0;
  virtual double dimension(void) = 0;
  virtual double time(void) = 0;
  virtual ~VectorAccess2D() {}
};

class Vector2DDisplay: public Display
{
protected:
  VectorAccess2D *cx;
  ofstream os;
  
public:
  Vector2DDisplay(VectorAccess2D *a, double period, double filePeriod)
    : Display(period,filePeriod), cx(a)
  {
    closeGnuplot();
    startGnuplotWithOutputFile("out/vector-log.gp");
    os.open(VECTORFILE);
    writeToGnuplot("set title 's plane'\n");
    writeToGnuplot("set nokey\n");
    writeToGnuplot("set parametric\n");
  }
  virtual ~Vector2DDisplay() {}
  void updateDisplayUnconditionally(void)
  {
    //os.seekp(0,ios::beg);
    // draw circle,
    writeToGnuplot("plot [0:6.28] [-5:5] [-5:5] cos(t), sin(t) w l 3,");
    // plot past phenotypes from file, current phenotypes from input, selection gradient from input
    //writeToGnuplot(" '"VECTORFILE"' w p 1, '-' w p 3\n");
    writeToGnuplot(" '"VECTORFILE"' w p 2, '-' w p 1, '-' w p 3\n");
    // supply the phenotypes
    for ( long i = 0; i < cx->count(); i++ )
      if ( cx->alive(i) )
      {  // plot each vector
	vnl_vector<double> si = cx->s(i);
	os << si[0] << ' ' << si[1] << endl;
	printfToGnuplot("%g %g %%%g\n", si[0], si[1], cx->X(i));
      }
    os << endl;
    os.flush();
    writeToGnuplot("e\n");
    // specify G
    vnl_matrix<double> &M = cx->M();
    vnl_vector<double> G(2,0);
    for ( long i = 0; i < cx->count(); i++ )
      if ( cx->alive(i) )
	G += cx->X(i) * M * cx->s(i);
    printfToGnuplot("%g %g\ne\n", G[0], G[1]);
    flushGnuplot();
  }
  void equilibrium(void)
  {
    updateDisplay();
  }
  void explosion(double t)
  {
    updateDisplayUnconditionally();
  }
};

#endif // VECTOR2DDISPLAY_H
