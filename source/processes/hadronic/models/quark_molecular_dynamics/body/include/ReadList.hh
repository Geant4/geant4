#ifndef __READLIST__
#define __READLIST__

#ifdef IS_GCC
#pragma interface
#endif

#include "globals.hh"
#include "g4std/iostream"
#include "g4std/vector"

//
//  Reads members of a class t from file.
//  '#' can be used for comments
//
//  Requirements to template class t:
//
//  - constructor  t(G4std::istream& )   for reading in properties
//  - operator<<   for being printed out
//  - operator==   to identify entry in the vector
//

template<class t>
class ReadList : public vector<t>
{
  friend G4std::ostream& operator<<(G4std::ostream&,ReadList<t>&);
protected:
  int N;
  G4std::istream* file;
public:
  ReadList(char* fileName);
  void readIn();
  virtual void bookIn(const t& item) { insert(end(),item); }
  virtual ~ReadList();
  int noSpecies() const { return N; }
  int getIndex(const t& h);
#ifdef IS_GCC
  ReadList<t>& operator=(const ReadList<t>&) { return *this; }
#endif
};

template<class t>
class ReadList_P : public vector<t*>
{
  friend G4std::ostream& operator<<(G4std::ostream&,ReadList_P<t>&);
protected:
  int N;
  G4std::istream* file;
public:
  ReadList_P(char* fileName);
  void readIn();
  virtual void bookIn(t* item) { insert(end(),item); }
  virtual ~ReadList_P();
  int noSpecies() const { return N; }
  int getIndex(const t& h);
#ifdef IS_GCC
  ReadList_P<t>& operator=(const ReadList_P<t>&) { return *this; }
#endif
};

#ifdef XLC_QNOTEMPINC
#ifndef IS_GCC
#include "ReadList.tcc"
#endif
#endif

#endif
