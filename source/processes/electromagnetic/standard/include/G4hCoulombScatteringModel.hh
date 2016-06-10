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
// $Id: G4hCoulombScatteringModel.hh 79067 2014-02-14 09:48:52Z gcosmo $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4hCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 08.06.2012 from G4eCoulombScatteringModel
//
// Modifications:
//
// Class Description:
//
// Implementation of Coulomb Scattering of pointlike charge particle 
// on Atomic Nucleus for interval of scattering anles in Lab system 
// thetaMin - ThetaMax, nucleus recoil is neglected.
//   The model based on analysis of J.M.Fernandez-Varea et al. 
// NIM B73(1993)447 originated from G.Wentzel Z.Phys. 40(1927)590 with 
// screening parameter from H.A.Bethe Phys. Rev. 89 (1953) 1256.
// 

// -------------------------------------------------------------------
//

#ifndef G4hCoulombScatteringModel_h
#define G4hCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "globals.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4WentzelVIRelXSection.hh"

class G4ParticleChangeForGamma;
class G4ParticleDefinition;
class G4IonTable;
class G4NistManager;

class G4hCoulombScatteringModel : public G4VEmModel
{

public:

  G4hCoulombScatteringModel(G4bool combined = true);
 
  virtual ~G4hCoulombScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
				G4double kinEnergy, 
				G4double Z, 
				G4double A, 
				G4double cut,
				G4double emax);

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy);

  // defines low energy limit of the model
  inline void SetLowEnergyThreshold(G4double val);

  // user definition of low-energy threshold of recoil
  inline void SetRecoilThreshold(G4double eth);

protected:

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  inline void SetupParticle(const G4ParticleDefinition*);

private:

  // hide assignment operator
  G4hCoulombScatteringModel & operator=(const G4hCoulombScatteringModel &right);
  G4hCoulombScatteringModel(const  G4hCoulombScatteringModel&);

  //protected:
 
  G4IonTable*               theIonTable;
  G4ParticleChangeForGamma* fParticleChange;
  G4WentzelVIRelXSection*   wokvi;
  G4NistManager*            fNistManager;

  const std::vector<G4double>* pCuts;

  const G4MaterialCutsCouple* currentCouple;
  const G4Material*           currentMaterial;
  G4int                       currentMaterialIndex;

  G4double                  cosThetaMin;
  G4double                  cosThetaMax;
  G4double                  cosTetMinNuc;
  G4double                  cosTetMaxNuc;
  G4double                  recoilThreshold;
  G4double                  elecRatio;
  G4double                  mass;

  // projectile
  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* theProton;

  G4double                  lowEnergyThreshold;

  G4bool                    isCombined;  

  //private:

  G4bool                    isInitialised;             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
void G4hCoulombScatteringModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4hCoulombScatteringModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    wokvi->SetupParticle(p);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4hCoulombScatteringModel::SetLowEnergyThreshold(G4double val)
{
  lowEnergyThreshold = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4hCoulombScatteringModel::SetRecoilThreshold(G4double eth)
{
  recoilThreshold = eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
