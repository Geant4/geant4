///////////////////////////////////////////////////////////////////////////////
// File: CCalG4Hit.hh
// Description: G4 Hit class for Calorimeters (Ecal, Hcal, ...)
///////////////////////////////////////////////////////////////////////////////

#ifndef CCalG4Hit_h
#define CCalG4Hit_h 1

#include "G4VHit.hh"

#include "CCalHit.hh"

class CCalG4Hit : public G4VHit, public CCalHit {

public:

  CCalG4Hit(){}
  ~CCalG4Hit(){}
  int  operator==(const CCalG4Hit &right) {return 0;}

  void Draw(){}
  void Print(){CCalHit::print();}
};
#endif
