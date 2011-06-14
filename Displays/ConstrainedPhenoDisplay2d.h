/* -*- C++ -*- */
#ifndef CONSTRDISPLAY_H
#define CONSTRDISPLAY_H

#include "DisplayController.h"
#include "GnuplotDisplay.h"
#include "Integrator.h"
#include "vnl/vnl_matrix.h"
#include <fstream>

class ConstrainedPhenoDisplay2d;

//#define VECTORFILE     "out/vectors.gp"
class Constraint2dController
  : public DisplayController<ConstrainedPhenoDisplay2d>
{
public:
  Constraint2dController(double period, double filePeriod, Site *s)
    : DisplayController<ConstrainedPhenoDisplay2d>(period,filePeriod,s)
  {}
  
  virtual ~Constraint2dController();

  void createDisplay();

  timescale_t timescale()
  { return evolutionary_time;
  }
  
  virtual string title(void)
  { return ""; }
  virtual string tmin(void)
  { return stringf("%g",mint()); }
  virtual string tmax(void)
  { return stringf("%g",maxt()); }
  virtual string xmin(void)
  { return stringf("%g",minx()); }
  virtual string xmax(void)
  { return stringf("%g",maxx()); }
  virtual string ymin(void)
  { return stringf("%g",miny()); }
  virtual string ymax(void)
  { return stringf("%g",maxy()); }
  virtual double mint(void)
  { return 0; }
  virtual double maxt(void)
  { return 1; }
  virtual double minx(void)
  { return 0; }
  virtual double maxx(void)
  { return 1; }
  virtual double miny(void)
  { return 0; }
  virtual double maxy(void)
  { return 1; }

  virtual bool approve(void) // whether to draw now or skip it
  { return true; }
  virtual int setsOfCurves(void) { return 0; }
  virtual int nCurves(int set) { return 0; }
  virtual int curveColor(int set, int i) { return 1; }
  virtual bool toDoCurve(int set, int no, int skip)// whether
  { return true; } 
  // give x,y as a function of t
  // example: gd.gnuplot << "t, 2*t*t";
  virtual void drawCurve(int set, int no, GnuplotDisplay&gd)
  {}
  // if the curve is given by a list of literal data points, put "'-'"
  // above and list them here, separated by '\n' and terminated by
  // 'e\n'.  don't get confused: there are 'points' below, which
  // appear disconnected, and there are 'curves'.  This function
  // specifies a curve by providing a list of connected points.
  virtual void plotPointsOfCurve(int set, int no, GnuplotDisplay&gd)
  {} 
  virtual int setsOfArrows(void){ return 0; }
  virtual int nArrows(int set) { return 0; }
  virtual int arrowColor(int set, int no) { return 1; }
  virtual bool toDoArrow(int set, int no, int skip)// whether
  { return true; }
  virtual vnl_vector<double> arrowBegin(int set, int no)
  { static vnl_vector<double> origin2d(2,0);
    return origin2d;
  }
  virtual vnl_vector<double> arrowLength(int set, int no)
  { static vnl_vector<double> origin2d;
    return origin2d;
  }
  virtual void customizeArrow(int set, int no, GnuplotDisplay&gw){}
  virtual int setsOfPoints(void) { return 0; }
  virtual int nPoints(int set) { return 0; }
  virtual int pointColor(int set) { return 3; }
  virtual bool toDoPoint(int set, int i, int skip) { return true; }
  virtual vnl_vector<double> point(int set,int i)
  { static vnl_vector<double> origin2d;
    return origin2d;
  }
  virtual bool keepHistory(void) { return false; }
  virtual vnl_vector<double> phenotype(int i)
  { static vnl_vector<double> origin2d;
    return origin2d;
  }
  virtual int index(const Index &ix)
  { return 0; } 
  virtual string logFile(void)
  { return outdir()+"/constraint-anim.gp"; } 
  virtual string historyFile(void) 
  { return outdir()+"/pheno-history.dat"; } 
//   virtual double time(void)
//   { return 0; }
};

class ConstrainedPhenoDisplay2d : public GnuplotDisplay
{
protected:
  Constraint2dController *cx;
  ofstream os;
  
public:
  ConstrainedPhenoDisplay2d(Constraint2dController *con)
    : cx(con)
  {}

  virtual ~ConstrainedPhenoDisplay2d() {}

  void initialize(void)
  { GnuplotDisplay::initialize();
    //updateDisplay();
  }

  string gnuplotLogFile(void)
  { return cx->logFile(); } 
  //  void updateDisplayUnconditionally(void)
  void updateDisplay(void)
  { updateDisplayWithout(-1); }
  void updateDisplayWithout(int skip)
  {
    if (!cx->approve())
    {
      clear();
      return;
    }

    rewindStream();
    writeGnuplotHeaders();
    gnuplot << "# t = " << cx->time() << "\n";
    gnuplot << "pause delay\n";
    gnuplot << "set noarrow\n";
    for (int set = 0; set < cx->setsOfArrows(); set++)
      for ( long i = 0; i < cx->nArrows(set); i++ )
	if ( cx->toDoArrow(set,i,skip) )
	{
	  vnl_vector<double> pt = cx->arrowBegin(set, i),
	    ar = cx->arrowLength(set, i);
	  vnl_vector<double> end = pt + ar;
// CLIMATE PROJECT MAY NEED THIS CODE
// 	  if (end[0] < 0)
// 	  { ar *= - pt[0] / ar[0];
// 	    end = pt + ar;
// 	  }
// 	  if (end[0] > 100)
// 	  { ar *= (100 - pt[0]) / ar[0];
// 	    end = pt + ar;
// 	  }
// 	  if (end[1] < 0)
// 	  { ar *= - pt[1] / ar[1];
// 	    end = pt + ar;
// 	  }
// 	  if (end[1] > 1)
// 	  { ar *= (1 - pt[1]) / ar[1];
// 	    end = pt + ar;
// 	  }	
	  gnuplot << "set arrow from " << pt[0] << ", " <<  pt[1]
		  << " to " << end[0] << ", " << end[1]
		  << " lt " << cx->arrowColor(set, i);
	  cx->customizeArrow(set, i, *this);
	  gnuplot << "\n";
	  gnuplot << "# length " << ar[0] << ", " << ar[1] << "\n";
	}
    bool comma=false;
    gnuplot << "plot [" << cx->tmin() << ":" << cx->tmax()
	    << "] [" << cx->xmin() << ':' << cx->xmax()
	    << "] [" << cx->ymin() << ':' << cx->ymax() << "] ";

    if (cx->keepHistory() && os.is_open())
    {
      if(comma) gnuplot << ", ";
      gnuplot << " '" << cx->historyFile() << "' w p 1";
      comma = true;
    }

    // put in the curves
    for (int set = 0; set < cx->setsOfCurves(); set++)
      for ( long i = 0; i < cx->nCurves(set); i++ )
      {
	if ( cx->toDoCurve(set,i,skip) )
	{
	  if(comma) gnuplot << ", ";
	  comma=true;
	  cx->drawCurve(set,i,*this);
	  gnuplot << " w l " << cx->curveColor(set,i);
	}
      }
    // first '-' is the history points, if any
    ostringstream collect;
    for (int set = 0; set < cx->setsOfPoints(); set++)
      for ( long i = 0; i < cx->nPoints(set); i++ )
      { bool todo = false;
	if ( cx->toDoPoint(set,i,skip) )
	{
	  vnl_vector<double> pt = cx->point(set,i);
	  collect << pt[0] << ' ' << pt[1] << '\n';
	  todo = true;
	}
	if (todo)
	{
	  collect << "e\n";
	  if(comma) gnuplot << ", ";
	  comma=true;
	  gnuplot << "'-' w p " << cx->pointColor(set);
	}
      }
    gnuplot << "\n" << collect.str();

    // drawCurve may also have entered '-' directives
    for (int set = 0; set < cx->setsOfCurves(); set++)
      for ( long i = 0; i < cx->nCurves(set); i++ )
      {
	if ( cx->toDoCurve(set,i,skip) )
	  cx->plotPointsOfCurve(set,i,*this);
      }
    plotBufferedData();
  }
  void extinction(int former)
  {
    if (cx->approve() && cx->keepHistory())
    { vnl_vector<double> pt = cx->phenotype(former);
      if (pt[0] != 0 || pt[1] != 0)
      {
	if (!os.is_open())
	  os.open(cx->historyFile().c_str());
	os << pt[0] << ' ' << pt[1] << endl;
	os.flush();
      }
    }
    updateDisplayWithout(former);
  }

  virtual void writeGnuplotHeaders() {
    string title = cx->title();
    if (title.length())
      gnuplot << "set title '" << title << "'\n";
    gnuplot << "set nokey\n";
    gnuplot << "set parametric\n";
    //gnuplot << "set pointsize 10\n";
    gnuplot << "delay=0\n";
  }
};

#endif // CONSTRDISPLAY_H
