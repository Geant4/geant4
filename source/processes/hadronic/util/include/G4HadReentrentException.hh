#ifndef G4HadReentrentException_h
#define G4HadReentrentException_h

#include <exception>
#include <iostream>
#include "G4HadronicException.hh"

class G4HadReentrentException : public G4HadronicException
{
  public:
  G4HadReentrentException(G4String in, G4int at, G4String mess)
  : G4HadronicException(in, at, mess) {}
  virtual ~G4HadReentrentException() throw () {}
  
  void Report(std::ostream & aS)
  {
    aS<< "G4HadReentrentException:" <<std::endl;
    G4HadronicException::Report(aS);
  }
};

#endif
