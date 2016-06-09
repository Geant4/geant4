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
// $Id: G4eCoulombScatteringModel.hh,v 1.49.2.1 2010/04/06 09:14:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-03-patch-02 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
//
// File name:     G4eCoulombScatteringModel
//
// Author:        Vladimir Ivanchenko 
//
// Creation date: 19.02.2006
//
// Modifications:
// 01.08.06 V.Ivanchenko extend upper limit of table to TeV and review the
//          logic of building - only elements from G4ElementTable
// 08.08.06 V.Ivanchenko build internal table in ekin scale, introduce faclim
// 19.08.06 V.Ivanchenko add inline function ScreeningParameter and
//                       make some members protected
// 09.10.07 V.Ivanchenko reorganized methods, add cut dependence in scattering off e- 
// 09.06.08 V.Ivanchenko add SelectIsotope and sampling of the recoil ion 
// 17.06.09 C.Consoalndi modified SetupTarget method - remove kinFactor
//				       
//
// Class Description:
//
// Implementation of eCoulombScattering of pointlike charge particle 
// on Atomic Nucleus for interval of scattering anles in Lab system 
// thetaMin - ThetaMax, nucleus recoil is neglected.
//   The model based on analysis of J.M.Fernandez-Varea et al. 
// NIM B73(1993)447 originated from G.Wentzel Z.Phys. 40(1927)590 with 
// screening parameter from H.A.Bethe Phys. Rev. 89 (1953) 1256.
// 

// -------------------------------------------------------------------
//

#ifndef G4eCoulombScatteringModel_h
#define G4eCoulombScatteringModel_h 1

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "globals.hh"
#include "G4NistManager.hh"
#include <vector>

class G4ParticleChangeForGamma;
class G4ParticleDefinition;

class G4eCoulombScatteringModel : public G4VEmModel
{

public:

  G4eCoulombScatteringModel(const G4String& nam = "eCoulombScattering");
 
  virtual ~G4eCoulombScatteringModel();

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
  inline void SetLowEnergyLimit(G4double val);

  // obsolete method
  inline void SetRecoilThreshold(G4double eth);

protected:

  G4double CrossSectionPerAtom();

  G4double SampleCosineTheta();

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void SetupKinematic(G4double kinEnergy, G4double cut);
  
  inline void SetupTarget(G4double Z, G4double kinEnergy); 

private:

  void ComputeMaxElectronScattering(G4double cut);

  // hide assignment operator
  G4eCoulombScatteringModel & operator=(const G4eCoulombScatteringModel &right);
  G4eCoulombScatteringModel(const  G4eCoulombScatteringModel&);

protected:
 
  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  G4ParticleTable*          theParticleTable; 
  G4ParticleChangeForGamma* fParticleChange;
  G4NistManager*            fNistManager;

  const std::vector<G4double>* pCuts;

  const G4MaterialCutsCouple* currentCouple;
  const G4Material*           currentMaterial;
  const G4Element*            currentElement;
  G4int                       currentMaterialIndex;

  G4double                  coeff;
  G4double                  cosThetaMin;
  G4double                  cosThetaMax;
  G4double                  cosTetMinNuc;
  G4double                  cosTetMaxNuc;
  G4double                  cosTetMaxNuc2;
  G4double                  cosTetMaxElec;
  G4double                  cosTetMaxElec2;
  G4double                  q2Limit;
  G4double                  recoilThreshold;
  G4double                  elecXSection;
  G4double                  nucXSection;
  G4double                  ecut;

  // projectile
  const G4ParticleDefinition* particle;

  G4double                  chargeSquare;
  G4double                  spin;
  G4double                  mass;
  G4double                  tkin;
  G4double                  mom2;
  G4double                  invbeta2;
  G4double                  etag;
  G4double                  lowEnergyLimit;

  // target
  G4double                  targetZ;
  G4double                  targetMass;
  G4double                  screenZ;
  G4double                  formfactA;
  G4int                     idxelm;
  G4int                     iz;

private:

  G4double                  alpha2;
  G4double                  faclim;

  static G4double ScreenRSquare[100];
  static G4double FormFactor[100];

  G4bool                    isInitialised;             
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline
void G4eCoulombScatteringModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4eCoulombScatteringModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    spin = particle->GetPDGSpin();
    G4double q = particle->GetPDGCharge()/eplus;
    chargeSquare = q*q;
    tkin = 0.0;
    etag = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eCoulombScatteringModel::SetupKinematic(G4double ekin, 
						      G4double cut)
{
  if(ekin != tkin || ecut != cut) {
    tkin = ekin;
    mom2 = tkin*(tkin + 2.0*mass);
    invbeta2 = 1.0 +  mass*mass/mom2;
    cosTetMinNuc = cosThetaMin;
    cosTetMaxNuc = cosThetaMax;
    if(mass < MeV && cosThetaMin < 1.0 && ekin <= 10.*cut) {
      cosTetMinNuc = ekin*(cosThetaMin + 1.0)/(10.*cut) - 1.0;
    }
    ComputeMaxElectronScattering(cut);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4eCoulombScatteringModel::SetupTarget(G4double Z, G4double e)
{
  if(Z != targetZ || e != etag) {
    etag    = e; 
    targetZ = Z;
    iz= G4int(Z);
    if(iz > 99) iz = 99;
    targetMass = fNistManager->GetAtomicMassAmu(iz)*amu_c2;
    screenZ = ScreenRSquare[iz]/mom2;
    screenZ *=(1.13 + std::min(1.0,3.76*Z*Z*invbeta2*alpha2));
    if(mass > MeV) { screenZ *= 2.0; }
    formfactA = FormFactor[iz]*mom2;
    cosTetMaxNuc2 = cosTetMaxNuc;
    if(1 == iz && particle == theProton && cosTetMaxNuc2 < 0.0) {
      cosTetMaxNuc2 = 0.0;
    }
    cosTetMaxElec2 = cosTetMaxElec;
  } 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eCoulombScatteringModel::SetLowEnergyLimit(G4double val)
{
  lowEnergyLimit = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4eCoulombScatteringModel::SetRecoilThreshold(G4double eth)
{
  recoilThreshold = eth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif
