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
//      File name:     G4ContinuumGammaDeexcitation
//
//      Authors:       Carlo Dallapiccola (dallapiccola@umdhep.umd.edu)
//                     Maria Grazia Pia (pia@genova.infn.it)
// 
//      Creation date: 23 October 1998
//
//      Modifications: 
//      
// -------------------------------------------------------------------
//
//
// Class G4ContinuumGammaDeexcitation.hh
//
#ifndef G4ContinuumGammaDeexcitation_hh
#define G4ContinuumGammaDeexcitation_hh 

#include "G4VGammaDeexcitation.hh"

#include "globals.hh"

#include "G4ContinuumGammaTransition.hh"
#include "G4Fragment.hh"
#include "G4NuclearLevelManager.hh"

class G4ContinuumGammaDeexcitation : public G4VGammaDeexcitation
{
public:

  // Constructor
  G4ContinuumGammaDeexcitation();

  // Destructor
  ~G4ContinuumGammaDeexcitation();

  // Functions

public:

  virtual G4VGammaTransition* CreateTransition();

  virtual G4bool CanDoTransition() const;

private:

  G4int _nucleusZ;
  G4int _nucleusA;  
  G4NuclearLevelManager _levelManager;


};

#endif
