/*
 * ODE parallel main program
 *
 * i.e. the main() function is run for each thread
 *
 * 4/28/2000  Lee Worden
 *
 * ...does this still have any relation to mpi? 8/31/2006
*/
#include "Simulation.h"
#include "Parameters.h"
#include <string>

void handleArgs(int argc, char **argv)
{
  for ( ++argv; (*argv); ++argv )
  {
    string arg(*argv);
    if (arg=="-f")
    {
      //cout << "argument: -f " << argv[1] << endl;
      string filename(*++argv);
      parameters.parseSettingsFile(filename);
    }
    else if (arg[0] == '-' && arg[1] == '-')
    { // this currently can only be --param=val
      string line(arg,2,string::npos);
      string::size_type eq = line.find('=');
      if (eq != string::npos)
        line[eq] = ' ';
      //cout << line << endl;
      parameters.parseLine(line);
    }
    else
    {
      cout << "unknown argument: " << arg << endl;
    }
  }
}

int main(int argc, char **argv)
{
  //parameters.loadDefaults();
  handleArgs(argc,argv);
  parameters.afterSetting();
  simulation.setup();
  simulation.doSimulation(parameters.runLength());
  simulation.finish();
  exit(0);
}
