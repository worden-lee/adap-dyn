#include "OutputController.h"
#include "Node.h"
#include "util.h"
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdarg.h>

void NodeOutputController::writeTag(ostream &o)
{
  if ( parameters.totalGridSize() > 1 )
  {
    o << "[" << parameters.rank() << "] ";
  }
}

void NodeOutputController::logTag(void)
{
  writeTag(logfile);
  if ( parameters.logToCout() )
    writeTag(cout);
}

void NodeOutputController::logWithAllCommunities(char*fmt,...)
{
  va_list va;
  va_start(va,fmt);
  logTag();
  logRawV(fmt,va);
  va_end(va);
  if ( parameters.outputCommunities() )
  { va_start(va,fmt);
    const char *v = vfstring(fmt,va);
    va_end(va);
    communities << v;
    for(int i=0; i<((Node*)site)->nSites; i++)
      ((Node*)site)->sites[i]->outputcontroller->communities << v;
  }
}

void NodeOutputController::logRawV(const char*fmt,va_list va)
{ const char *v = vfstring(fmt,va);
  logfile << v;
  if ( parameters.logToCout() )
    cout << v;
}

void NodeOutputController::initDir(void)
{
  openMemberFiles();
  
  log("init node, rank = %i, size = %i\n",
      parameters.rank(), parameters.size());
  logfile.flush();
}

void NodeOutputController::reopenMemberFiles(void)
{
  int rank = parameters.rank();

  // "" is a fine directory
//   if ( strlen(dir) == 0 )
//   {
//     cerr << "outputcontroller not correctly initialized\n";
//     cerr.flush();
//   }  
  
  {  // try to create the output directory
    string path( parameters.outputDirectory() + DIR_SEP + dir);
    (void)mkdir(path.c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
    
    struct stat st;
    if (stat(path.c_str(),&st))
      cerr << "[" << rank << "] problem with directory '" << path << "'\n";
    else if ( !S_ISDIR(st.st_mode) )
      cerr << "[" << rank << "] '" << path << "' is not a directory\n";
  }

  string path( parameters.outputDirectory() + DIR_SEP + dir +
	       DIR_SEP + "log." + int_to_string(rank) );
  logfile.open(path.c_str());
  if ( !logfile.is_open() )
  {
    cerr << "[" << rank << "] Failed to open log file (" << path << ")\n";
    cerr.flush();
  }

  if ( parameters.outputCommunities() )
  {
    path = parameters.outputDirectory() + DIR_SEP + dir +
      DIR_SEP + "community." + int_to_string(rank);

    communities.open(path.c_str());
    if ( !communities.is_open() )
    {
      cerr << "[" << rank << "] Failed to open community file ("
	   << path << ")\n";
      cerr.flush();
    }
  }
}
