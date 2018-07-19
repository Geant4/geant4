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
// $Id: G4hCoulombScatteringModel.hh 104307 2017-05-24 09:01:45Z gcosmo $
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
// Implementation of Coulomb Scattering of a charge particle 
// on Atomic Nucleus for interval of scattering anles in Lab system 
// thetaMin - ThetaMax.
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

  explicit G4hCoulombScatteringModel(G4bool combined = true);
 
  virtual ~G4hCoulombScatteringModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual void InitialiseLocal(const G4ParticleDefinition*, 
                               G4VEmModel* masterModel) override;

  virtual G4double ComputeCrossSectionPerAtom(
                                const G4ParticleDefinition*,
				G4double kinEnergy, 
				G4double Z, 
				G4double A, 
				G4double cut,
				G4double emax) override;

  virtual void SampleSecondaries(std::vector<G4DynamicParticle*>*,
				 const G4MaterialCutsCouple*,
				 const G4DynamicParticle*,
				 G4double tmin,
				 G4double maxEnergy) override;

  virtual G4double MinPrimaryEnergy(const G4Material*,
				    const G4ParticleDefinition*,
				    G4double) final;

  // defines low energy limit of the model
  inline void SetLowEnergyThreshold(G4double val);

  // user definition of low-energy threshold of recoil
  inline void SetRecoilThreshold(G4double eth);

  // defines low energy limit on energy transfer to atomic electron
  inline void SetFixedCut(G4double);

  // low energy limit on energy transfer to atomic electron
  inline G4double GetFixedCut() const;

protected:

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  inline void SetupParticle(const G4ParticleDefinition*);

private:

  // hide assignment operator
  G4hCoulombScatteringModel & operator=
  (const G4hCoulombScatteringModel &right) = delete;
  G4hCoulombScatteringModel(const  G4hCoulombScatteringModel&) = delete;
 
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
  G4double                  recoilThreshold;
  G4double                  elecRatio;
  G4double                  mass;

  G4double                  fixedCut;

  // projectile
  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* theProton;

  G4bool                    isCombined;  
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

inline void G4hCoulombScatteringModel::SetRecoilThreshold(G4double eth)
{
  recoilThreshold = eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4hCoulombScatteringModel::SetFixedCut(G4double val)
{
  fixedCut = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline G4double G4hCoulombScatteringModel::GetFixedCut() const
{
  return fixedCut;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
