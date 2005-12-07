#ifndef _NTSTFileRead_
#define _NTSTFileRead_ 1

#include <strstream.h>
#include "globals.hh"

#include <fstream.h>

class NTSTFileRead{

public:
  NTSTFileRead(const char* FileName, G4bool echo=false);
  ~NTSTFileRead();
  char* ReadLine();
  istrstream &StreamLine();
  
private:
  char _Line[255];
  int _LineLength;
  ifstream* _Istr;
  G4bool _echo;
  
  istrstream *stuff;
};

#endif
