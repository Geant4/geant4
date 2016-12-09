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
// $Id: G4IonFluctuations.hh 100399 2016-10-20 07:38:12Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4IonFluctuations
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 02.04.2003
//
// Modifications: 
//
// 16-10-03 Changed interface to Initialisation (V.Ivanchenko)
// 01-06-08 Added initialisation of effective charge prestep (V.Ivanchenko)
//
// Class Description:
//
// Implementation of ion energy loss fluctuations

// -------------------------------------------------------------------
//

#ifndef G4IonFluctuations_h
#define G4IonFluctuations_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4UniversalFluctuation.hh"

class G4Pow;

class G4IonFluctuations : public G4VEmFluctuationModel
{

public:

  explicit G4IonFluctuations(const G4String& nam = "IonFluc");

  virtual ~G4IonFluctuations();

  // Sample fluctuations
  virtual G4double SampleFluctuations(const G4MaterialCutsCouple*,
                                      const G4DynamicParticle*,
                                      G4double tmax,
                                      G4double length,
                                      G4double meanLoss) override;

  // Compute dispertion 
  virtual G4double Dispersion(const G4Material*,
                              const G4DynamicParticle*,
                              G4double tmax,
                              G4double length) override;

  // Initialisation prerun
  virtual void InitialiseMe(const G4ParticleDefinition*) override;

  // Initialisation prestep
  virtual void SetParticleAndCharge(const G4ParticleDefinition*, 
                                    G4double q2) override;

private:

  G4double Factor(const G4Material*, G4double Zeff);
  G4double RelativisticFactor(const G4Material*, G4double Zeff);

  // hide assignment operator
  G4IonFluctuations & operator=(const  G4IonFluctuations &right) = delete;
  G4IonFluctuations(const  G4IonFluctuations&) = delete;

  G4UniversalFluctuation      uniFluct;
  const G4ParticleDefinition* particle;

  G4Pow*   g4calc; 

  G4double particleMass;
  G4double charge;
  G4double chargeSquare;
  G4double effChargeSquare;

  // data members to speed up the fluctuation calculation
  G4double parameter;
  G4double theBohrBeta2;
  G4double minFraction;
  G4double xmin;
  G4double minLoss;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

