#ifndef G4HadronicException_h
#define G4HadronicException_h

#include <exception>
#include <iostream>
#include "globals.hh"

class G4HadronicException : public std::exception
{
  public:
  G4HadronicException(G4String in, G4int at, G4String mess)
  {
    theMessage = mess;
    theName = in;
    theLine = at;
  }
  virtual ~G4HadronicException() throw () {}
  
  void Report(std::ostream & aS)
  {
    aS<< "In " <<theName<<", line "<<theLine<<": "<<theMessage<<std::endl;
  }
  
  private:
  G4String theMessage;
  G4String theName;
  G4int theLine;
};

#endif
