#ifndef __GENREAD__
#define __GENREAD__

#include <iostream>
#include <fstream>
#include "String.hh"

class genericRead
{
  std::istream* in;
  String Name;
  virtual void readIn(std::istream&) = 0;
public:
  genericRead(std::istream& i) : in(&i) {}
  genericRead(char* file);
  virtual ~genericRead() { delete in; }
  void read();
};

template<class t>
class FileRead : private genericRead
{
  virtual void readIn(std::istream& in) { t*x = new t(in); }
public:
  FileRead(const String& fileName) 
    : genericRead(fileName) { genericRead::read(); }
};

#endif
