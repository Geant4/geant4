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
// File name:     G4eeToTwoGammaModel
//
// Author:        Vladimir Ivanchenko on base of Michel Maire code
//
// Creation date: 02.08.2004
//
// Modifications:
// 08-04-05 Major optimisation of internal interfaces (V.Ivanchenko)
// 18-04-05 Compute CrossSectionPerVolume (V.Ivanchenko)
// 06-02-06 ComputeCrossSectionPerElectron, ComputeCrossSectionPerAtom (mma)
// 20-10-06 Add theGamma as a member (V.Ivanchenko)
//
//
// Class Description:
//
// Implementation of e+ annihilation into 2 gamma

// -------------------------------------------------------------------
//

#ifndef G4eeToTwoGammaModel_h
#define G4eeToTwoGammaModel_h 1

#include "G4VEmModel.hh"

class G4ParticleChangeForGamma;

class G4eeToTwoGammaModel : public G4VEmModel
{

public:

  explicit G4eeToTwoGammaModel(const G4ParticleDefinition* p = nullptr,
                               const G4String& nam = "eplus2gg");

  ~G4eeToTwoGammaModel() override;

  void Initialise(const G4ParticleDefinition*, const G4DataVector&) override;
  
  virtual G4double ComputeCrossSectionPerElectron(G4double kinEnergy); 
  
  G4double ComputeCrossSectionPerAtom(
                                 const G4ParticleDefinition*,
                                 G4double kinEnergy, 
                                 G4double Z, 
                                 G4double A = 0., 
                                 G4double cutEnergy = 0.,
                                 G4double maxEnergy = DBL_MAX) override;

  G4double CrossSectionPerVolume(const G4Material*,
				 const G4ParticleDefinition*,
				 G4double kineticEnergy,
				 G4double cutEnergy = 0.0,
				 G4double maxEnergy = DBL_MAX) override;

  void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  // hide assignment operator
  G4eeToTwoGammaModel & operator=(const  G4eeToTwoGammaModel &right) = delete;
  G4eeToTwoGammaModel(const  G4eeToTwoGammaModel&) = delete;

private:

  G4double pi_rcl2;
  const G4ParticleDefinition* theGamma;
  G4ParticleChangeForGamma* fParticleChange;

  static G4bool fSampleAtomicPDF;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
