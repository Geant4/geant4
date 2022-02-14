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

class G4Pow;
class G4UniversalFluctuation;

class G4IonFluctuations : public G4VEmFluctuationModel
{

public:

  explicit G4IonFluctuations(const G4String& nam = "IonFluc");

  ~G4IonFluctuations() override;

  // Sample fluctuations
  G4double SampleFluctuations(const G4MaterialCutsCouple*,
			      const G4DynamicParticle*,
			      const G4double tcut,
			      const G4double tmax,
			      const G4double length,
			      const G4double meanLoss) override;

  // Compute dispertion 
  G4double Dispersion(const G4Material*,
		      const G4DynamicParticle*,
                      const G4double tcut,
		      const G4double tmax,
		      const G4double length) override;

  // Initialisation prerun
  void InitialiseMe(const G4ParticleDefinition*) override;

  // Initialisation prestep
  void SetParticleAndCharge(const G4ParticleDefinition*, 
			    G4double q2) override;

  // hide assignment operator
  G4IonFluctuations & operator=(const  G4IonFluctuations &right) = delete;
  G4IonFluctuations(const  G4IonFluctuations&) = delete;

private:

  G4double Factor(const G4Material*, G4double Zeff);
  G4double RelativisticFactor(const G4Material*, G4double Zeff);

  const G4ParticleDefinition* particle = nullptr;
  G4UniversalFluctuation* uniFluct;
  G4Pow* g4calc; 

  G4double particleMass;
  G4double charge = 1.0;
  G4double chargeSquare = 1.0;
  G4double effChargeSquare = 1.0;

  // data members to speed up the fluctuation calculation
  G4double parameter;
  G4double theBohrBeta2;
  G4double minFraction = 0.2;
  G4double xmin = 0.2;
  G4double minLoss;
  // cash
  G4double kineticEnergy = 0.0;
  G4double beta2 = 0.0;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

