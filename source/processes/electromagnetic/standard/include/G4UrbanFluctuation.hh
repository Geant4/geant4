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
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4UrbanFluctuation
//
// Author:        V.Ivanchenko make a class with the Laszlo Urban model
//
// Creation date: 14.02.2022
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss fluctuations made by L.Urban for 
// Geant4 10.X series for updated design of 11.X by V.Ivanchenko

// -------------------------------------------------------------------
//

#ifndef G4UrbanFluctuation_h
#define G4UrbanFluctuation_h 1

#include "G4UniversalFluctuation.hh"

class G4UrbanFluctuation : public G4UniversalFluctuation
{

public:

  explicit G4UrbanFluctuation(const G4String& nam = "UrbanFluc");

  ~G4UrbanFluctuation() override;

  // hide assignment operator
  G4UrbanFluctuation & operator =
  (const G4UrbanFluctuation &right) = delete;
  G4UrbanFluctuation(const G4UrbanFluctuation&) = delete;

protected:

  G4double SampleGlandz(CLHEP::HepRandomEngine* rndm, 
                        const G4Material*, const G4double tcut) override;

private:

  // material properties
  const G4Material* lastMaterial = nullptr;
  G4double f1Fluct = 0.0;
  G4double f2Fluct = 0.0;
  G4double e1Fluct = 0.0;
  G4double e2Fluct = 0.0;
  G4double e1LogFluct = 0.0;
  G4double e2LogFluct = 0.0;
  G4double esmall = 0.0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
