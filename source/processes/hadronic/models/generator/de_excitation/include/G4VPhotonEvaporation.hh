// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class file 
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group
//      CERN, Geneva, Switzerland
//
//      File name:     G4VPhotonEvaporation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------

#ifndef G4VPHOTONEVAPORATION_HH
#define G4VPHOTONEVAPORATION_HH

#include "globals.hh"
#include "G4Fragment.hh"

class G4VPhotonEvaporation 
{
public:

  G4VPhotonEvaporation() {};
  virtual ~G4VPhotonEvaporation() {};
  
  G4bool operator==(const G4VPhotonEvaporation &right) const;
  G4bool operator!=(const G4VPhotonEvaporation &right) const;
  
  virtual G4FragmentVector* BreakItUp(const G4Fragment &theNucleus) = 0;
  
private:  

  G4VPhotonEvaporation(const G4VPhotonEvaporation &right);
    const G4VPhotonEvaporation& operator=(const G4VPhotonEvaporation &right);

};

#endif
