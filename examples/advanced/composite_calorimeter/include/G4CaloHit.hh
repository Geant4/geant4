///////////////////////////////////////////////////////////////////////////////
// File: G4CaloHit.hh
// Date: 11.1998 
// Modifications:
///////////////////////////////////////////////////////////////////////////////

#ifndef G4CaloHit_h
#define G4CaloHit_h 1

#include "G4VHit.hh"

#include "CMSCaloHit.hh"

class G4CaloHit : public G4VHit, public CMSCaloHit {

public:

  G4CaloHit(){}
  ~G4CaloHit(){}
  int  operator==(const G4CaloHit &right){return 0;}

  void Draw(){}
  void Print(){CMSCaloHit::print();}
};
#endif
