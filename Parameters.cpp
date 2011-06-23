#include "Parameters.h"
#include <fstream>
#include <list>
#include <functional>
#include <sys/stat.h>

Parameters::Parameters(Parameters*upstream)  : inheritFrom(upstream)
{
  if (!inheritFrom)
  { // if settings/default.settings is found, use it.
    // if someone invokes the program from a different directory, they
    // have to provide the settings using -f.
    const char *def_filename = "settings/default.settings";
    struct stat statbuf;
    int status = stat(def_filename, &statbuf);
    if (status == 0)
    { setdefaultSettingsFile(def_filename);
      loadDefaults();
    }
  }  
  // finishConstruct(); 
}

void Parameters::parseSettingsFile(string filename)
{
  cout << "parseSettingsFile(\"" << filename << "\")\n";
  ifstream settings(filename.c_str());
  if (settings.is_open())
    parseSettings(settings);
  else
  { cout << "couldn't open " << filename << endl;
    exit(-1);
  }
}

// each line here must be
// "key val"
// with just one space between the two;
// or, "#comment"; or, "include filename";
// otherwise the line will be ignored or misinterpreted.
// leading and trailing spaces are ignored.
void Parameters::parseSettings(istream &file)
{
  string line;
  while (getline(file,line))
  {
    parseLine(line);
  }
}

string Parameters::cleanLine(string line)
{
  // get rid of anything starting with #
  string::size_type pound = line.find('#');
  // or //
  string::size_type slashes = line.find("//");
  if (pound == string::npos || slashes < pound)
    pound = slashes;
  if (pound != string::npos)
    line.erase(pound);

  // get rid of leading + trailing space
  string::iterator e = line.end();
  while (e != line.begin() && isspace(e[-1]))
    --e;
  string::iterator b = line.begin();
  while (b != e && isspace(*b))
    ++b;
  //line.erase(e,line.end());
  line.replace(line.begin(),line.end(),b,e);

  return line;
}

void Parameters::parseLine(string line)
{
  line = cleanLine(line);

  // check for include directive
  string directive("include ");
  if (line.substr(0,directive.size())==directive)
  {
    parseSettingsFile(line.substr(directive.size(),string::npos));
    return;
  }

  // what's left has to have a space as separator
  string::size_type space = line.find(' ');
  if (space == string::npos)
    return;
  string key(line,0,space);
  string val(line,space+1);
  cout << key << "|" << val << endl;
  set(key,val);
}

void Parameters::writeAllSettings(ostream&outs)
{
  list<string> keys;
  // get all keys from the dictionary
  transform(begin(),end(), inserter(keys,keys.begin()),
	    _Select1st< pair<string,string> >());
  // output key/value pairs using sorted order
  for (list<string>::iterator ki = keys.begin(); ki != keys.end(); ++ki)
    if (const string *getki = get(*ki))
      outs << *ki << ' ' << *getki << '\n';
}
