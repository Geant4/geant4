// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LowEnergyPolarizedCompton.hh,v 1.1 2001-05-23 16:39:40 flongo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ------------------------------------------------------------
//      GEANT 4 class header file
//      CERN Geneva Switzerland
//
//      For information related to this code contact:
//      CERN, IT Division, ASD group

// --------- G4LowEnergyPolarizedCompton class -----
//
//           by G.Depaola & F.Longo (21 may 2001)
//
// ************************************************************
// ------------------------------------------------------------

// class description
//
// inherit from G4LowEnergyCompton
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4LowEnergyPolarizedCompton_h
#define G4LowEnergyPolarizedCompton_h 1

#include "G4LowEnergyCompton.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4LowEnergyPolarizedCompton : public G4LowEnergyCompton
{  
 public:  // with description
  
  G4LowEnergyPolarizedCompton(const G4String& processName = "polarLowEnCompt");
  
  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);
  
 private:
  
  G4ThreeVector SetNewPolarization(G4double,G4double, G4double, G4double);
  
  G4double SetPhi(G4double, G4double);
  
  void SystemOfRefChange(G4ThreeVector&, G4ThreeVector&, G4ThreeVector&, 
                         G4ThreeVector&);

  // hide assignment operator as private 
  G4LowEnergyPolarizedCompton& operator=(const G4LowEnergyPolarizedCompton &right);
  G4LowEnergyPolarizedCompton(const G4LowEnergyPolarizedCompton& );

};

#endif
 








