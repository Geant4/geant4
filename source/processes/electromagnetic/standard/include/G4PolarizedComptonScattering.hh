//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PolarizedComptonScattering.hh,v 1.8 2005/05/04 16:16:12 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// --------- G4PolarizedComptonScattering physics process ----------------------
//                   by Vicente Lara, March 1998
//
// -----------------------------------------------------------------------------

// class description
//
// inherit from G4ComptonScattering
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef G4PolarizedComptonScattering_h
#define G4PolarizedComptonScattering_h 1

#include "G4ComptonScattering52.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4PolarizedComptonScattering : public G4ComptonScattering52
{
 public:  // with description

  G4PolarizedComptonScattering(const G4String& processName = "polarCompt");

  G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep);

 private:
  
  G4ThreeVector SetNewPolarization(G4double, G4double,
  				   G4double, G4double,
  				   G4ThreeVector& );

  G4double SetPhi(G4double, G4double, G4double, G4double);


  void SystemOfRefChange(G4ThreeVector&, G4ThreeVector&, G4ThreeVector&, 
                         G4ThreeVector&);

  // hide assignment operator as private 
  G4PolarizedComptonScattering& operator=(const G4PolarizedComptonScattering &right);
  G4PolarizedComptonScattering(const G4PolarizedComptonScattering& );

};
  
#endif
 
