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
// File name:     G4AtimaFluctuations
//
// Author:        Jose Luis Rodriguez Sanchez on base of ATIMA code
//
// Creation date: 16.01.2018
//
// Modifications:
//
//
// Class Description:
//
// Implementation of ion energy loss fluctuations
//
// -------------------------------------------------------------------
//

#ifndef G4AtimaFluctuations_h
#define G4AtimaFluctuations_h 1


#include "G4VEmFluctuationModel.hh"
#include "G4ParticleDefinition.hh"
#include "G4UniversalFluctuation.hh"

class G4Pow;

class G4AtimaFluctuations : public G4VEmFluctuationModel
{

public:

  explicit G4AtimaFluctuations(const G4String& nam = "IonFlucAtima");

  virtual ~G4AtimaFluctuations();

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

  G4double EnergyTable_interpolate(const G4double* table,G4double xval, const G4double* y);

  // hide assignment operator
  G4AtimaFluctuations & operator=(const  G4AtimaFluctuations &right) = delete;
  G4AtimaFluctuations(const  G4AtimaFluctuations&) = delete;

  G4UniversalFluctuation      uniFluct;
  const G4ParticleDefinition* particle;

  G4Pow*   g4calc; 

  G4double particleMass;
  G4double charge;
  G4double chargeSquare;
  G4double effChargeSquare;

  G4double MLN10;
  G4double atomic_mass_unit;
  G4double dedx_constant;
  G4double electron_mass;
  G4double fine_structure;
  G4double domega2dx_constant;

  static G4double stepE;
  static G4double tableE[200];
  static const G4double ls_X_coefficients_a[110][200];
  static const G4double ls_X_coefficients_ahi[110][200];
  static const G4double element_atomic_weights[110];

  // data members to speed up the fluctuation calculation
  G4double minLoss;
  // cash
  G4double kineticEnergy;
  G4double beta2;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

