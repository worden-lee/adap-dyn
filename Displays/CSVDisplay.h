/* -*- C++ -*- */
// #include "Integrator.h"
// #include <math.h> 
// #include "OutputController.h"
// #include "GnuplotDisplay.h"
// #include "CVIntegrator.h"
// #include "Site.h"
#include <string>
#include <fstream>
#include <iostream>

class CSVDisplay : public Display
{
//  friend class CSVController;
private:
  ofstream logfile;
  bool onBlankLine;
  //  string myOutfile;

public:
  CSVDisplay () {
    // included for compatibility, but a file is not open!
    // use openFile if this must be used
  }

  // preferred constructor form - pass in the filename to use
  CSVDisplay(string filename) {
    logfile.open(filename.c_str());
    onBlankLine = true;
    //    myOutfile = filename;
  }

  ~CSVDisplay () {
    logfile.close();
  }

  // open a file, if there's not one loaded already
  void openFile(string filename) {
    if (logfile) {
      cerr << "CSVDisplay::openFile called, but a file is already open!\n";
    } else {
      logfile.open(filename.c_str());
      onBlankLine = true;
      //      myOutfile = filename;
    }
  }

  // write a whole line (i.e., header data)
  CSVDisplay &writeLine(const string line) {
    if (!logfile) {
       cerr << "CSVDisplay::writeLine called, but no logfile is open!\n";
    }
    logfile << line;
    // check if the line included a newline, and write as needed
    string::const_reverse_iterator lastchar = line.rbegin();
    if (*lastchar != '\n') logfile << endl;
    onBlankLine = true;
    return *this;
  }

  // Write a new datum -- should do comma insertion automagically
  // this should automatically handle other numericals
  CSVDisplay &operator<<(const double& data) { 
    if (&data == 0 || data == HUGE) {
      // HUGE: the international code for "it's broken!"
      *this << string("NA");
      return *this;
    }
    if (!logfile) {
      cerr << "CSVDisplay::<< called, but no logfile is open!\n";
    }
    if (onBlankLine) {
      onBlankLine = false;
    } else {
      logfile << ",";
    }
    logfile << data;
    return *this;
  }

  CSVDisplay &operator<<(const string data) {
    if (!logfile) {
      cerr << "CSVDisplay::<< called, but no logfile is open!\n";
    }
    if (onBlankLine) {
      onBlankLine = false;
    } else {
      logfile << ",";
    }
    logfile << data;
    return *this;
  }

  // start a new row of data
  // writes a new line only if we're not on a fresh line
  void newRow() {
    logfile << endl;
    onBlankLine = true;
  }

public:
  // included for compatibility, but I don't think these should be
  // called directly...
  void updateDisplay(void) {} 
  void recordCommunity(void) {}

  void flush() {
    logfile.flush();
  }
};

#if 0 // need some work to make this file compatible in climate/
//template<class CSVDisplay>
class CSVController : public DisplayController<CSVDisplay>
{
protected:
  int* speciationsForEachLineage; // for an array
  int equilibriaReached;
  
public:
  // We're stretching the DisplayController spec a little, because
  // this class will control many CSVDisplay classes: one that is
  // continually updated with general data as the simulation runs, and
  // many that contain snapshots of data for each live species.
  // 
  // So *display only points to the community logfile CSVDisplay. We
  // initialize the snapshot writers as needed.
  //
  // Of course there is no displayPeriod needed, because the CSVs are
  // not displayed, only written to disk.
  CSVController(double filePeriod, Site *s=0)
    : DisplayController<CSVDisplay>(HUGE,filePeriod,s)
    {
      // create the community logfile
      createDisplay();
      // vars this controller keeps track of
#define MAXLINEAGES 1000
      speciationsForEachLineage = new int[MAXLINEAGES];
      equilibriaReached = 0;
      recordCounter = time();
    }

  string csvdir(){
    return outdir()+"/csv";
  }

  void createDisplay(){
    mkdir(csvdir().c_str(), S_IRWXU|S_IRWXG|S_IRWXO);
    string communityOutfile = csvdir() + "/community-log.csv";
    display = new CSVDisplay(communityOutfile);
    display->writeLine("temp,eco_t,evo_t,n_equil,n_species,n_resources,dT/dmaxT");
  }

  // The default is fine for most of these notifications
  // but for speciations and equilibria we keep a running count  
  void speciation(double t, const Index& parent, const Index& daughter)
  { 
    int lin = ((LCommunity*)site->community)->lineage(parent);
    if ( lin < 0 || lin > MAXLINEAGES )
      site->outputcontroller->log("CSVController - speciation in invalid lineage\n");
    else
      ++speciationsForEachLineage[lin];
    //    update(time());     
  }

  void equilibrium(double t) { 
    equilibriaReached++;
    update(time());
  }

  void step(double t) {}

  // more junk for completeness
  void flush(void) {
    display->flush();
  }
  void finish(void) {
    delete display;
  }
  void recordCommunity(void) {
    update(time());
  }
  void updateDisplay(void) {}

  void update(double t)
  {
//     if (!display)
//       createDisplay();
    if (recordEvery > 0 && recordCounter+recordEvery <= t && t != 0)
      {
	updateCommunityState();
	recordSnapshot();
	recordCounter = t;
      }
  }

  // update the community logfile
  void updateCommunityState()
  {
    LCommunity &comm = *(LCommunity*)site->community;
    //CVIntegrator &integ = *(CVIntegrator*)site->integrator;
    // record temperature
    *display << site->integrator->state(comm.temperatureIndex());
    // record ecological time
    *display << site->integrator->time();
    // record evolutionary time
    *display << site->integrator->evolutionaryTime();
    // # equilibria
    *display << equilibriaReached;
    // # species
    *display << comm.nSpecies();
    // # resources
    *display << comm.nResources(); 
    // the relevant community graph? Skip this for now
    // record dT/dmaxT - the effect of changing max temp
    *display << ((LSite*)site)->calculations->dT_dMaxT;
    display->newRow();
    display->flush();
    
  }

  // record a snapshot: a line of variable for each live species
  void recordSnapshot(void) 
  {
    // quick refs
    LCommunity &comm = *(LCommunity*)site->community;
    //CVIntegrator &integ = *(CVIntegrator*)site->integrator;
    Calculations &calc = *((LSite*)site)->calculations; 
    // globals
    double t = site->integrator->time();
    double t_ev = site->integrator->evolutionaryTime();
    double T = site->integrator->state(comm.temperatureIndex());
    // Open a new file for the snapshot
    string snapbase = 
      csvdir() + "/snapshot-" + stringf("%g", t_ev) + '-' + stringf("%g", t);
    string snapshotfile = snapbase + ".csv";
    // TODO: if the file already exists, do something about it?
    // for now, just overwrite it...
    CSVDisplay* snapshot = new CSVDisplay(snapshotfile);
    // write a header line
    string header = "lineage,species,lineage_mutations,dT_dtau,dRstar_dT,"
      "pop,T0,Rstar,stability,dT_without,extinct_without";
    snapshot->writeLine(header);

    calc._recache();
    for (int i=0; i<comm.nSpecies(); i++) {
      //      if (liveSpeciesIndexPtrs[i] != 0) {
      //	Index currentSpeciesIndex = *liveSpeciesIndexPtrs[i];
      if (!comm.alive[i]) continue;
      species &currentSpecies = comm.sp[i];
      Index iIndex(i,comm.speciesIndexing());
      int currentLineage = comm.lineage(iIndex);
      //int speciesIndexKey = iIndex.key();
      int speciesUniqueIndexKey = comm.speciesUniqueIndexing().index(iIndex);
      // record which # lineage
      *snapshot << currentLineage;
      // record which species
      //	*snapshot << iIndex.key();
      *snapshot << speciesUniqueIndexKey;
      // record # mutations so far
      *snapshot << speciationsForEachLineage[currentLineage];
      // record which interval this lineage started in (do we care?)
      // -- for now, I think I'll just skip this bit
      // -- it'd be nice if one of those cues told us when a new lin started
      // TODO: These next ones might be returning reciprocals..
      if(lparameters.doPowerRelations())
      {
	// record dtau_i/dT
	*snapshot << calc.dT_dtau[i];
	// record dR*_i/dT
	*snapshot << calc.dRstar_dT[i];
      }
      // record pop(i)
      *snapshot << site->integrator->state(iIndex);
      // record T0_i
      *snapshot << currentSpecies.T0;
      // record R*_i
      *snapshot << calc.Rstar(T, currentSpecies.T0, currentSpecies.gamma);
      if(lparameters.doPowerRelations())
      {
	// record stability contribution calculated
	*snapshot << calc.stabilityContributions[i];
	// effect of this species disappearing
	*snapshot << calc.deltaTWithoutEach[i];
	*snapshot << calc.extinctionsWithoutEach[i];
      }

      snapshot->newRow();
    }
    // Close the file
    snapshot->flush();
    delete snapshot;

    // now write a text file about the community
    string commfile = snapbase + ".txt";
    ofstream co(commfile.c_str());
    co << "Community at evolutionary time " << t_ev
       << " / ecological time " << t << "\n\n";

    co << "Temperature: " << T << '\n';
    co << "Dead temperature: "
       << lparameters.maxTemperature()*comm._resources[0].heating << '\n';
    co << "Dominant eigenvalue: "
       << static_cast<LSite*>(site)->calculations->baselineStability << '\n'
       << comm.speciesCount() << " species\n"
       << '\n';

    for (int i=0; i<comm.nSpecies(); i++) 
      if (comm.alive[i])
      {
	species &spi = comm.sp[i];
	Index iIndex(i,comm.speciesIndexing());
	int iLineage = comm.lineage(iIndex);
    
	co << "Species " << site->outputcontroller->basename(iIndex) 
	   << " (" << site->outputcontroller->basename(spi.source())
	   << "->" << site->outputcontroller->basename(spi.product())
	   << ")\n"
	   << "lineage: " << iLineage << '\n'
	   << "population: " << site->integrator->state(iIndex) << '\n'
	   << "tau: " << spi.T0 << '\n'
	   << "R*: " << calc.Rstar(T, spi.T0, spi.gamma) << '\n';
        if (lparameters.doPowerRelations())
          co << "dR*/dT: " << calc.dRstar_dT[i] << '\n'
             << "dT/dtau: " << calc.dT_dtau[i] << '\n'
             << "stability contribution: "
             << calc.stabilityContributions[i] - calc.baselineStability << '\n'
             << "if removed: deltaT = " << calc.deltaTWithoutEach[i]
             << ", " << 1+calc.extinctionsWithoutEach[i]
             << " extinctions\n"
             << '\n';
	
        if (lparameters.doPowerRelations())
          for (int j=0; j<comm.nSpecies(); j++) 
            if (comm.alive[j])
            {
              //species &jSpecies = comm.sp[j];
              Index jIndex(j,comm.speciesIndexing());
              
              co << "Effect on " << site->outputcontroller->basename(jIndex)
                 << ":\n"
                 << "    dR*[" << comm.speciesUniqueIndexing().index(jIndex)
                 << "]/dtau[" << comm.speciesUniqueIndexing().index(iIndex)
                 << "] = " << calc.dRstar_dT[j]*calc.dT_dtau[i] << '\n';
              if (j!=i)
                co << "    if removed: "
                   << (calc.causesExtinction[i][j]? "goes extinct":"survives")
                   << '\n';
              co << '\n';
            }
      }
  }
};
#endif
