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
// File name:     G4BraggIonModel
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 13.10.2004
//
// Modifications:
// 09-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 11-05-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 25-04-06 Add stopping data from ASTAR (V.Ivanchenko)
// 12-08-08 Added methods GetParticleCharge, GetChargeSquareRatio, 
//          CorrectionsAlongStep needed for ions(V.Ivanchenko)

//
// Class Description:
//
// Implementation of energy loss and delta-electron production
// by heavy slow charged particles using ICRU'49, NIST, and ICRU90 
// evaluated data for alpha and protons

// -------------------------------------------------------------------
//

#ifndef G4BraggIonModel_h
#define G4BraggIonModel_h 1

#include "G4BraggModel.hh"

class G4ASTARStopping;

class G4BraggIonModel : public G4BraggModel
{

public:

  explicit G4BraggIonModel(const G4ParticleDefinition* p = nullptr,
			   const G4String& nam = "BraggIon");

  ~G4BraggIonModel() override;

  void Initialise(const G4ParticleDefinition*, 
		  const G4DataVector&) override;

  G4double ComputeCrossSectionPerAtom(
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double Z, G4double A,
				 G4double cutEnergy,
				 G4double maxEnergy) override;

  G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy,
				 G4double maxEnergy) override;

  G4double ComputeDEDXPerVolume(const G4Material*,
                                const G4ParticleDefinition*,
                                G4double kineticEnergy,
                                G4double cutEnergy) override;

  // Compute ion charge not applied to alpha
  G4double GetChargeSquareRatio(const G4ParticleDefinition*,
				const G4Material*,
				G4double kineticEnergy) override;

  // add correction to energy loss and ompute non-ionizing energy loss
  void CorrectionsAlongStep(const G4MaterialCutsCouple*,
			    const G4DynamicParticle*,
			    const G4double& length,
			    G4double& eloss) override;

  // hide assignment operator
  G4BraggIonModel & operator=(const  G4BraggIonModel &right) = delete;
  G4BraggIonModel(const  G4BraggIonModel&) = delete;

private:

  G4double HeEffChargeSquare(const G4double z, 
                             const G4double kinEnergyInMeV) const;

  G4int HasMaterialForHe(const G4Material* material) const;

  G4double HeStoppingPower(const G4double kinEnergy) const;

  G4double HeElectronicStoppingPower(const G4int z, const G4double kinEnergy) const;

  G4double HeDEDX(const G4Material* material, const G4double kinEnergy);

  static G4ASTARStopping* fASTAR;

  G4double heChargeSquare = 4.0;
  G4double HeMass;
  G4double massFactor;

  G4int iASTAR = -1;    // index in ASTAR
  G4bool isAlpha = false;
  G4bool isFirstAlpha = false;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
