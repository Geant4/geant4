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
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4DynamicParticleFluctuation
//
// Author:        V.Ivanchenko
//
// Creation date: 23.08.2024
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss fluctuations using G4DynamicParticle
// information and without any access to G4ParticleDefinition. 

// -------------------------------------------------------------------
//

#ifndef G4DynamicParticleFluctuation_h
#define G4DynamicParticleFluctuation_h 1

#include "G4UniversalFluctuation.hh"
#include "G4DynamicParticle.hh"

class G4DynamicParticleFluctuation : public G4UniversalFluctuation
{

public:

  G4DynamicParticleFluctuation(const G4String& nam = "dynPartFluc");

  ~G4DynamicParticleFluctuation() override = default;

  G4double SampleFluctuations(const G4MaterialCutsCouple*,
			      const G4DynamicParticle*,
                              const G4double, const G4double,
			      const G4double, const G4double) override;

  G4double Dispersion(const G4Material*,
		      const G4DynamicParticle*,
                      const G4double, const G4double,
                      const G4double) override;

  // hide assignment operator
  G4DynamicParticleFluctuation & operator=
  (const G4DynamicParticleFluctuation &right) = delete;
  G4DynamicParticleFluctuation(const G4DynamicParticleFluctuation&) = delete;

private:

  void InitialiseLocal(const G4DynamicParticle*);
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

