/* -*- C++ -*- 
 */
#ifndef PARAMETERS_H
#define PARAMETERS_H
#include "Indexing.h"
#include "util.h"
#include <time.h>
#include <math.h>

#ifdef __unix__
#define UNIX
#define DIR_SEP '/'
#else
#define WINDOWS
#define DIR_SEP '\\'
#endif

class Parameters;
// global parameter object
extern Parameters &parameters;

typedef map<string, string> dictionary_t;

class Parameters : protected dictionary_t
{
protected:
  inline void set(string key, string val)
  { operator[](key) = val; }
  inline const string *get(string key) const
  { const_iterator it=find(key); \
    return (it==end()) ? 0 : &it->second; \
  }
  // this is not needed by the program but useful in debugger
  inline const string *get(const char *key) const
  { return get(string(key));
  }

// Declaring parameters at compile time:
// they are all represented internally by associating a
// value string with a key string in the dictionary,
// accessed via set() and get().
// But we also provide typed set/get functions for
// the rest of the program to use.
// So for instance,
// 
// DECLARE_PARAM(double,abortOnFullMoon)
// 
// would expand to
// 
// inline void setabortOnFullMoon(double val)
// { set("abortOnFullMoon",double_to_string(val)); }
// inline double abortOnFullMoon()
// { const string *str = get("abortOnFullMoon");
//   return str? string_to_double(*str) : double();
// }
//
// the various string_to_* and *_to_string functions are
// declared in util.h.  If some other type is needed those
// may have to be added to.
  
// RTYPE is helpful in case of enum types, apparently
#define DECLARE_PARAM_3(RTYPE,TYPE,NAME) \
  inline void set##NAME(TYPE val) \
  { set(#NAME, TYPE##_to_string(val)); } \
  inline RTYPE NAME() \
  { const string *str = get(#NAME);\
    return str? string_to_##RTYPE(*str) : RTYPE(); \
  }

#define DECLARE_PARAM(TYPE,NAME) \
  DECLARE_PARAM_3(TYPE,TYPE,NAME)

// here's some nasty stuff for making the get/set
// business work with enum types.  We should
// probably just disallow enums and use string
// values instead.

// brought over from the climate project, 9/13/2006

// macros here to generate string_to_* and *_to_string
// functions for enum types
//
// say you have enum todo_t { CONTINUE, PANIC };
// then ENUM_DECLARATIONS(todo_t) gives you
//  map<string,todo_t> s_todo_t_map;
//  map<todo_t,string> todo_t_s_map;
//  inline todo_t string_to_todo_t(const string &x)
//  { return s_todo_t_map[x]; }
//  inline string todo_t_to_string(en_t x)
//  { return todo_t_s_map[x]; }
//
// and
//  DECLARE_ENUM_VALUE(todo_t,CONTINUE)
//  DECLARE_ENUM_VALUE(todo_t,PANIC)
// gives you
//  s_todo_t_map["CONTINUE"] = CONTINUE;
//  todo_t_s_map[CONTINUE] = "CONTINUE";
//  s_todo_t_map["PANIC"] = PANIC;
//  todo_t_s_map[PANIC] = "PANIC";

// this one has to go in the class declaration,
// once for each enum type
#define ENUM_DECLARATIONS(en_t) \
  map<string,en_t> s_##en_t##_map; \
  map<en_t,string> en_t##_s_map; \
  inline en_t string_to_##en_t(const string &x) \
  { return s_##en_t##_map[x]; } \
  inline string en_t##_to_string(en_t x) \
  { return en_t##_s_map[x]; }

// and this one has to go in the LParameters constructor,
// once for each possible value of each enum type.
#define DECLARE_ENUM_VALUE(en_t,val) \
  s_##en_t##_map[#val] = val; \
  en_t##_s_map[val] = #val;

public:
  // extinction
  DECLARE_PARAM(bool, doExtinction);
  DECLARE_PARAM(double, extinctionThreshold);
  DECLARE_PARAM(double, explosionCeiling);
  DECLARE_PARAM(bool, doExtinctionNotificationTheOldWay);
  // speciation
  DECLARE_PARAM(bool, doSpeciation);
  DECLARE_PARAM(double, introductionSize);
  DECLARE_PARAM(bool, speciateAtEquilibrium);
  DECLARE_PARAM(double, equilibriumThreshold);
  DECLARE_PARAM(double, equilibriumTime);
  DECLARE_PARAM(double, maxWaitForEquilibrium);
  DECLARE_PARAM(double, speciationRate);
  DECLARE_PARAM(bool, useInvasionProbability);
  DECLARE_PARAM(double, invasionProbabilityFactor);
  DECLARE_PARAM(bool, usePopulationInInvasionProbability);
  //DECLARE_PARAM(long, giveUpAfter);
//   DECLARE_PARAM(long, parentTriesBeforeESS);
//   DECLARE_PARAM(long, speciationMutationTriesPerParent);
  DECLARE_PARAM(long, mutationTriesBeforeESS);
  DECLARE_PARAM(bool, quitAfterReachingESS);
  DECLARE_PARAM(bool, useBoundedMutation);
  DECLARE_PARAM(bool, useClamping);
  DECLARE_PARAM(bool, preserveSigns);
  // invasion (assembly model)
  DECLARE_PARAM(bool, doImmigration);
  DECLARE_PARAM(double, immigrationRate);
  // space
  DECLARE_PARAM(int, rank);
  DECLARE_PARAM(int, size);
  DECLARE_PARAM(int, nProcessors);
  DECLARE_PARAM(int, rowLength);
  DECLARE_PARAM(int, sitesPerProcessor);
  DECLARE_PARAM(int, totalGridSize);
  DECLARE_PARAM(int, nRows);
  //  int    rowsPerProcessor;
  // diffusion
  DECLARE_PARAM(double, diffusionConstant);
  DECLARE_PARAM(double, diffusionTimeStep);
  // output
  DECLARE_PARAM(double, outputTimeStep);
  DECLARE_PARAM(string, outputDirectory);
  DECLARE_PARAM(double, reopenOutputFilesInterval);
  DECLARE_PARAM(bool, reallyWrite);
  DECLARE_PARAM(bool, outputTiming);
  DECLARE_PARAM(bool, outputValues); //, refers, to, graphs of populations etc.
  DECLARE_PARAM(bool, outputMPIMessages);
  DECLARE_PARAM(bool, outputODEMessages);
  DECLARE_PARAM(bool, outputInternalMessages);
  DECLARE_PARAM(bool, outputInitMessage);
  DECLARE_PARAM(bool, outputSummaries);
  DECLARE_PARAM(double, minCountForAvg);
  DECLARE_PARAM(bool, outputSpeciations);
  DECLARE_PARAM(bool, outputImmigrations);
  DECLARE_PARAM(bool, outputExtinctions);
  DECLARE_PARAM(bool, outputCommunities);
  DECLARE_PARAM(bool, logToCout);
  DECLARE_PARAM(bool, disableDisplaying);
  DECLARE_PARAM(string, ghostviewCommand);
//  DECLARE_PARAM(bool, runGnuplot);
//  DECLARE_PARAM(bool, useGnuplotLogFile);
//  DECLARE_PARAM(string, gnuplotLogFilename);
  // whatever
  DECLARE_PARAM(long, randSeed);
  DECLARE_PARAM(double, runLength);
  DECLARE_PARAM(bool, useAuxDataFiles);
  DECLARE_PARAM(double, integrationGrain); // used by CVIntegrator
  DECLARE_PARAM(string, defaultSettingsFile);
  
  // parameter values get set in here
  //Parameters();
  // this is some kind of unholy hybrid of old style and new style. fix it.
  Parameters *inheritFrom;
  Parameters(Parameters*upstream=0);

  virtual ~Parameters() {}

  // get some settings from a file
  virtual void loadDefaults(void)
  { parseSettingsFile(defaultSettingsFile()); } 
  void parseSettingsFile(string filename);
  void parseSettings(istream &);
  virtual string cleanLine(string);
  virtual void parseLine(string);

  // write them back out to a file
  virtual void writeAllSettings(ostream&);

  virtual void rankAndSizeAre(int r, int s)
  {
    setrank(r); setnProcessors(s); setsize(s);
    settotalGridSize( nProcessors() * sitesPerProcessor() );
    setnRows((int)ceil(((double)totalGridSize()) / rowLength()));
    //      rowsPerProcessor = sitesPerProcessor / rowLength; 
    Indexing::stagger(r,s);
  }

  // this is called after all settings are read in
  virtual void afterSetting(void)
  { if (randSeed() == 0.0)
      setrandSeed(time(0));
  }

  // this is called after all the objects are set up and initialized
  virtual void finishInitialize(void)
  { }
};

#endif //PARAMETERS_H
