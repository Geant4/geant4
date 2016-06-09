//$Id: genconf.cpp,v 1.35 2008/10/15 21:51:24 marcocle Exp $	//

#ifdef _WIN32
  // Disable a warning in Boost program_options headers:
  // inconsistent linkage in program_options/variables_map.hpp
  #pragma warning ( disable : 4273 )
  #define popen _popen
  #define pclose _pclose 
  #define fileno _fileno 
  #include <stdlib.h>
#endif

// Include files----------------------------------------------------------------
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include "LibSymbolInfo.h"

using namespace std;

namespace windef {
  void usage(){
    cerr << "Usage: genwindef [-l <dllname>] [-o <output-file> | exports.def]  <obj or lib filenames>" << endl;
    exit(1);
  }
}

//--- Command main program-----------------------------------------------------
int main ( int argc, char** argv )
//-----------------------------------------------------------------------------
{
  string outfile("exports.def");
  string library("UnknownLib");
  string objfiles;
  bool debug(false);

  int arg;
  if (argc < 3) windef::usage();
  arg = 1;
  while (argv[arg][0] == '-') {
    if (strcmp(argv[arg], "--") == 0) {
      windef::usage();
    } 
    else if (strcmp(argv[arg], "-l") == 0) {
      arg++; 
      if (arg == argc) windef::usage();
      library = argv[arg];
    } 
    else if (strcmp(argv[arg], "-o") == 0) {
      arg++; 
      if (arg == argc) windef::usage();
      outfile = argv[arg];
    } 
    arg++;
  }
  if (arg == argc) windef::usage();
  for (arg; arg < argc; arg++) {
     objfiles += argv[arg];
     if( arg+1 < argc) objfiles += " ";
  }

  CLibSymbolInfo libsymbols;
  ofstream out(outfile.c_str());
  if(out.fail()) {
    cerr << "windef: Error opening file " << outfile << endl;
    return 1;
  }
  out << "LIBRARY " << library << endl;
  out << "EXPORTS" << endl;

  libsymbols.DumpSymbols(const_cast<char*>(objfiles.c_str()), out);

  out.close();


  return 0;
}



