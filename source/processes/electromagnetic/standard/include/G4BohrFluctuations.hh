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
// $Id: G4BohrFluctuations.hh 96934 2016-05-18 09:10:41Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4BohrFluctuations
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications: 
//
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
//
// Class Description:
//
// Implementation of Gaussion energy loss Fluctuations

// -------------------------------------------------------------------
//

#ifndef G4BohrFluctuations_h
#define G4BohrFluctuations_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4Material.hh"
#include "G4DynamicParticle.hh"

class G4BohrFluctuations : public G4VEmFluctuationModel
{

public:

  explicit G4BohrFluctuations(const G4String& nam = "BohrFluc");

  virtual ~G4BohrFluctuations();

  virtual G4double SampleFluctuations(const G4MaterialCutsCouple*,
                                      const G4DynamicParticle*,
                                      G4double, G4double, G4double) override;

  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
                              G4double, G4double) override;

  virtual void InitialiseMe(const G4ParticleDefinition*) override;

private:

  // hide assignment operator
  G4BohrFluctuations & operator=(const  G4BohrFluctuations &right) = delete;
  G4BohrFluctuations(const  G4BohrFluctuations&) = delete;

  const G4ParticleDefinition* particle;

  G4double particleMass;
  G4double chargeSquare;

  G4double minNumberInteractionsBohr;
  G4double minFraction;
  G4double xmin;
  G4double minLoss;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

