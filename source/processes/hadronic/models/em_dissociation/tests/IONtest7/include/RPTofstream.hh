#ifndef RPTofstream_h
#define RPTofstream_h 1
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// MODULE:		RPTofstream.hh
//
// Version:		0.a1
// Date:
// Author:		P R Truscott
// Organisation:	QinetiQ Ltd
// Customer:		ESA/ESTEC, NOORDWIJK
// Contract:		12115/JG/NL Work Order No. 3
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
//
// DESCRIPTION
// -----------
//
//
// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//
// PUBLIC MEMBER FUNCTIONS
// -----------------------
//
////////////////////////////////////////////////////////////////////////////////
//
#include "globals.hh"
#include "G4UnitsTable.hh"
#include <strstream>
#include "fstream.h"
#include <string>
////////////////////////////////////////////////////////////////////////////////
//
class RPTofstream : public ofstream
{
  public:
    RPTofstream () : ofstream() {};
  //    RPTofstream (const char *name, int mode=ios::out, int prot=0664)
  //      : ofstream(name, mode, prot) {};
    RPTofstream (const char *name)
        : ofstream(name) {};
    ~RPTofstream () {};

  public:
    void outG4String (G4String, G4int);
    void outG4BestUnit (G4BestUnit, G4int);
//    friend RPTofstream& operator<< (RPTofstream &, G4string &);
//    friend RPTofstream& operator<< (RPTofstream &, string &);
//    friend RPTofstream& operator<< (RPTofstream &, char );
};
////////////////////////////////////////////////////////////////////////////////
//
inline void RPTofstream::outG4String (G4String q, G4int w)
{
  *this << q.substr(0,w);
  for (G4int i = q.length(); i < w; i++) *this << " ";
}
////////////////////////////////////////////////////////////////////////////////
//
inline void RPTofstream::outG4BestUnit (G4BestUnit q, G4int w)
{
  char tmpChar[80] = {'0'};
  std::ostrstream os(tmpChar,80);
  os <<q;

  this->outG4String(G4String(tmpChar),w);
}
////////////////////////////////////////////////////////////////////////////////
//
/*//inline RPTofstream& operator << (RPTofstream &s, G4String &q)
inline RPTofstream& operator << (RPTofstream &s, string &q)
//inline RPTofstream& operator << (RPTofstream &s, char c)
{
//  string q = string(&c);
  G4int w = s.width();
  G4cout <<"GOT HERE!!!!" <<G4endl;
  ofstream *ofs = (ofstream*) &s;
  *ofs << q.substr(0,w);
  for (G4int i = q.length(); i < w; i++) *ofs <<" ";
  return s;
}*/
////////////////////////////////////////////////////////////////////////////////
#endif
