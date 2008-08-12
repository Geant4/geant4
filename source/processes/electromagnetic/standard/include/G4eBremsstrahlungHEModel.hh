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
// $Id: G4eBremsstrahlungHEModel.hh,v 1.2 2008-08-12 17:50:23 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eBremsstrahlungHEModel
//                extention of standard G4eBremsstrahlungModel
//
// Author:        Andreas Schaelicke 
//
// Creation date: 28.03.2008
//
// Modifications:
//
//
// Class Description:
//
// Implementation of energy loss for gamma emission by electrons and
// positrons including an improved version of the LPM effect

// -------------------------------------------------------------------
//

#ifndef G4eBremsstrahlungHEModel_h
#define G4eBremsstrahlungHEModel_h 1

#include "G4VEmModel.hh"

class G4Element;
class G4ParticleChangeForLoss;

class G4eBremsstrahlungHEModel : public G4VEmModel
{

public:

  G4eBremsstrahlungHEModel(const G4ParticleDefinition* p = 0, 
			   const G4String& nam = "eBremHE");

  virtual ~G4eBremsstrahlungHEModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double MinEnergyCut(const G4ParticleDefinition*, 
			const G4MaterialCutsCouple*);

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

  inline void SetLPMflag(G4bool val);
  inline G4bool LPMflag() const;

  inline void SetEnergyThreshold(G4double val);
  inline G4double EnergyThreshold() const;

  G4double SupressionFunction(const G4Material* material, G4double tkin,
                                    G4double gammaEnergy);

protected:

  inline G4double MaxSecondaryEnergy(const G4ParticleDefinition*,
				     G4double kineticEnergy);

  const G4Element* SelectRandomAtom(const G4MaterialCutsCouple* couple);

private:

  void SetParticle(const G4ParticleDefinition* p);

  G4double ComputeBremLoss(G4double Z, G4double tkin, G4double cut);

  G4double PositronCorrFactorLoss(G4double Z, G4double tkin, G4double cut);

  G4double PositronCorrFactorSigma(G4double Z, G4double tkin, G4double cut);

  G4DataVector* ComputePartialSumSigma(const G4Material* material,
                                             G4double tkin, G4double cut);

  inline G4double ScreenFunction1(G4double ScreenVariable);

  inline G4double ScreenFunction2(G4double ScreenVariable);

  // hide assignment operator
  G4eBremsstrahlungHEModel & operator=(const  G4eBremsstrahlungHEModel &right);
  G4eBremsstrahlungHEModel(const  G4eBremsstrahlungHEModel&);

protected:

  const G4ParticleDefinition* particle;
  G4ParticleDefinition*       theGamma;
  G4ParticleChangeForLoss*    fParticleChange;

  G4double minThreshold;
  G4bool   isElectron;

private:

  G4double highKinEnergy;
  G4double lowKinEnergy;
  G4double probsup;
  G4double MigdalConstant;
  G4double LPMconstant;
  G4double highEnergyTh;
  G4bool   theLPMflag;
  G4bool   isInitialised;

  std::vector<G4DataVector*> partialSumSigma;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4eBremsstrahlungHEModel::ScreenFunction1(G4double ScreenVariable)

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
G4double G4eBremsstrahlungHEModel::ScreenFunction2(G4double ScreenVariable)

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

inline 
G4double G4eBremsstrahlungHEModel::MaxSecondaryEnergy(
                                 const G4ParticleDefinition*,
    				       G4double kineticEnergy)
{
  return kineticEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eBremsstrahlungHEModel::SetLPMflag(G4bool val) 
{
  G4cout<<"G4eBremsstrahlungHEModel::SetLPMflag("<<val<<")\n";
  theLPMflag = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4bool G4eBremsstrahlungHEModel::LPMflag() const 
{
  return theLPMflag;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eBremsstrahlungHEModel::SetEnergyThreshold(G4double val) 
{
  highEnergyTh = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
G4double G4eBremsstrahlungHEModel::EnergyThreshold() const 
{
  return highEnergyTh;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
