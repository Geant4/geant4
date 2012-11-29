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
// $Id$
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremsstrahlungModel
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 07.01.2002
//
// Modifications:
//
// 23-12-02 Change interface in order to move to cut per region (V.Ivanchenko)
// 24-01-03 Make models region aware (V.Ivanchenko)
// 13-02-03 Add name (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 07-02-06  public function ComputeCrossSectionPerAtom() (mma)
//
//
// Class Description:
//
// Implementation of energy loss for gamma emission by electrons and
// positrons

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlungModel_h
#define G4eBremsstrahlungModel_h 1

#include "G4VEmModel.hh"

class G4Element;
class G4ParticleChangeForLoss;

class G4eBremsstrahlungModel : public G4VEmModel
{

public:

  G4eBremsstrahlungModel(const G4ParticleDefinition* p = 0, 
			 const G4String& nam = "eBrem");

  virtual ~G4eBremsstrahlungModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeDEDXPerVolume(const G4Material*,
					const G4ParticleDefinition*,
					G4double kineticEnergy,
					G4double cutEnergy);
					
  virtual G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
                                                  G4double tkin, 
                                                  G4double Z,   G4double,
                                                  G4double cut,
						  G4double maxE = DBL_MAX);
  
  virtual G4double CrossSectionPerVolume(const G4Material*,
					 const G4ParticleDefinition*,
					 G4double kineticEnergy,
					 G4double cutEnergy,
					 G4double maxEnergy);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

protected:

  const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple);

private:

  void SetParticle(const G4ParticleDefinition* p);

  G4double ComputeBremLoss(G4double Z, G4double tkin, G4double cut);

  G4double PositronCorrFactorLoss(G4double Z, G4double tkin, G4double cut);

  G4double PositronCorrFactorSigma(G4double Z, G4double tkin, G4double cut);

  G4DataVector* ComputePartialSumSigma(const G4Material* material,
                                             G4double tkin, G4double cut);

  G4double SupressionFunction(const G4Material* material, G4double tkin,
                                    G4double gammaEnergy);

  inline G4double ScreenFunction1(G4double ScreenVariable);

  inline G4double ScreenFunction2(G4double ScreenVariable);

  // hide assignment operator
  G4eBremsstrahlungModel & operator=(const  G4eBremsstrahlungModel &right);
  G4eBremsstrahlungModel(const  G4eBremsstrahlungModel&);

protected:

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  G4bool   isElectron;

private:

  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double probsup;
  G4double MigdalConstant;
  G4double LPMconstant;
  G4bool   isInitialised;

  std::vector<G4DataVector*> partialSumSigma;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eBremsstrahlungModel::ScreenFunction1(G4double ScreenVariable)

// compute the value of the screening function 3*PHI1 - PHI2

{
  G4double screenVal;

  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*std::log(ScreenVariable+0.952);
  else
    screenVal = 42.392 - ScreenVariable* (7.796 - 1.961*ScreenVariable);

  return screenVal;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
G4double G4eBremsstrahlungModel::ScreenFunction2(G4double ScreenVariable)

// compute the value of the screening function 1.5*PHI1 - 0.5*PHI2

{
  G4double screenVal;

  if (ScreenVariable > 1.)
    screenVal = 42.24 - 8.368*std::log(ScreenVariable+0.952);
  else
    screenVal = 41.734 - ScreenVariable* (6.484 - 1.250*ScreenVariable);

  return screenVal;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
