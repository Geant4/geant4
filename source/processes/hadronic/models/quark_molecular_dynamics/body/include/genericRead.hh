#ifndef __GENREAD__
#define __GENREAD__

#include "g4std/iostream"
#include "g4std/fstream"
#include "String.hh"

class genericRead
{
  G4std::istream* in;
  String Name;
  virtual void readIn(G4std::istream&) = 0;
public:
  genericRead(G4std::istream& i) : in(&i) {}
  genericRead(char* file);
  virtual ~genericRead() { delete in; }
  void read();
};

template<class t>
class FileRead : private genericRead
{
  virtual void readIn(G4std::istream& in) { t*x = new t(in); }
public:
  FileRead(const String& fileName) 
    : genericRead(fileName) { genericRead::read(); }
};

#endif
