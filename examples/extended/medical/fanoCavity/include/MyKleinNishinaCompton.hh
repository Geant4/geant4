//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file medical/fanoCavity/include/MyKleinNishinaCompton.hh
/// \brief Definition of the MyKleinNishinaCompton class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef MyKleinNishinaCompton_h
#define MyKleinNishinaCompton_h 1

#include "G4KleinNishinaCompton.hh"

class DetectorConstruction;
class MyKleinNishinaMessenger;
class G4ParticleChangeForGamma;

class MyKleinNishinaCompton : public G4KleinNishinaCompton
{

public:

  MyKleinNishinaCompton(DetectorConstruction*,
                        const G4ParticleDefinition* p = 0, 
                        const G4String& nam = "myKlein-Nishina");

 ~MyKleinNishinaCompton();
                                      
  virtual G4double CrossSectionPerVolume(
                                const G4Material*,
                                const G4ParticleDefinition*,
                                      G4double kinEnergy, 
                                      G4double cut,
                                      G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
                                 const G4MaterialCutsCouple*,
                                 const G4DynamicParticle*,
                                 G4double tmin,
                                 G4double maxEnergy);
                                      
  void SetCSFactor(G4double factor) {fCrossSectionFactor = factor;};

protected:

  DetectorConstruction*    fDetector;
  MyKleinNishinaMessenger* fMessenger;
  G4double                 fCrossSectionFactor;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
