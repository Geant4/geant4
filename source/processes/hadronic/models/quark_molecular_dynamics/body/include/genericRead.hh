#ifndef __GENREAD__
#define __GENREAD__

#include <iostream.h>
#include <fstream.h>
#include "String.hh"

class genericRead
{
  istream* in;
  String Name;
  virtual void readIn(istream&) = 0;
public:
  genericRead(istream& i) : in(&i) {}
  genericRead(char* file);
  virtual ~genericRead() { delete in; }
  void read();
};

template<class t>
class FileRead : private genericRead
{
  virtual void readIn(istream& in) { t*x = new t(in); }
public:
  FileRead(const String& fileName) 
    : genericRead(fileName) { genericRead::read(); }
};

#endif
