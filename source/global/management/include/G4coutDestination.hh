#ifndef G4COUTDESTINATION_HH
#define G4COUTDESTINATION_HH

#include "globals.hh"

class G4coutDestination
{
public:
  virtual int ReceiveG4cout(G4String){return 0;}
  virtual int ReceiveG4cerr(G4String){return 0;}
};
#endif

