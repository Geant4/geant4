#include <strstream.h>
#include "genericRead.hh"
#include <ctype.h>
#include <stdlib.h>

genericRead::genericRead(char* file) : Name(file)
{
  try {
    in = new ifstream(file);
  }
  catch ( ... ) {
    cerr << "Cannot open file " << file << endl;
    throw;
  }
}

void genericRead::read() 
{
  const int maxLineLength = 120;
  char c;
  int noLines = 1,N=0;
  while ( *in ) {
    while ( in->get(c) && isspace(c) )
      if ( c=='\n' ) ++noLines;
    if ( !(*in) ) {
      cerr << "Read in complete...\n";
      break;
    }
    if ( c == '#' ) {
      char Line[maxLineLength];
      in->get(Line,maxLineLength,'\n');
      in->get(c);
      ++noLines;
    } 
    else{
      in->putback(c);
      try {
	readIn(*in);
	++N;
      }
      catch (const String& s) { 
	cerr << Name << ": Read Error in Line " << noLines << ":\n";
	cerr << s << endl;
	exit(1);
      }
    }
  }
}

