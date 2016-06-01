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
//      File name:     G4ContinuumGammaTransition
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
//  Header file for G4ContinuumGammaTransition
//

#ifndef G4ContinuumGammaTransition_hh
#define G4ContinuumGammaTransition_hh

#include "globals.hh"
#include "G4VGammaTransition.hh"
#include "G4NuclearLevelManager.hh"
#include "G4VLevelDensityParameter.hh"

class G4ContinuumGammaTransition : public G4VGammaTransition
{
public:

  // Constructor
  G4ContinuumGammaTransition(const G4NuclearLevelManager& levelManager,
			     G4int Z, G4int A, G4double excitation, G4int verbose);

  // Destructor
  ~G4ContinuumGammaTransition();

  // Functions

  virtual G4double GammaEnergy();
  virtual G4double GetEnergyTo() const;
  virtual void SetEnergyFrom(const G4double energy);

private:

  G4double E1Pdf(G4double energy);

  G4int _A;
  G4int _Z;
  G4double _eMin;
  G4double _eMax;
  G4double _maxLevelE;
  G4double _minLevelE;
  G4double _excitation;
  G4double _eGamma;
  G4NuclearLevelManager _levelManager;

};

#endif
