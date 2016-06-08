// This code implementation is the intellectual property of
// the GEANT4 collaboration.
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
//      File name:     G4DiscreteGammaDeexcitation
//
//      Author:        Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------
//
#ifndef G4DiscreteGammaDeexcitation_hh
#define G4DiscreteGammaDeexcitation_hh 

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"

#include "G4DiscreteGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4NuclearLevelManager.hh"

class G4DiscreteGammaDeexcitation : public G4VGammaDeexcitation
{
public:

  // Constructor
  G4DiscreteGammaDeexcitation();

  // Destructor
  ~G4DiscreteGammaDeexcitation();

  // Functions

public:

  virtual G4VGammaTransition* CreateTransition();

  virtual G4bool CanDoTransition() const;

private:
  G4int _nucleusZ;
  G4int _nucleusA;
  G4double _tolerance;
  G4NuclearLevelManager _levelManager;
};


#endif
