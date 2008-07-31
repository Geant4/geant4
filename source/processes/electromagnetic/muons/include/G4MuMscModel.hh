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
// $Id: G4MuMscModel.hh,v 1.15 2008-07-31 13:11:57 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4MuMscModel
//
// Author:        V.Ivanchenko on base of L.Urban model
//
// Creation date: 25.10.2007
//
// Modifications:
//
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model
// J.M. Fernandez-Varea et al., NIM B73 (1993) 447;
// G.Wentzel, Z. Phys. 40 (1927) 590.

// -------------------------------------------------------------------
//

#ifndef G4MuMscModel_h
#define G4MuMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VMscModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4NistManager.hh"

class G4LossTableManager;
class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4ParticleDefinition;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuMscModel : public G4VMscModel
{

public:

  G4MuMscModel(const G4String& nam = "MuMscUni");

  virtual ~G4MuMscModel();

  void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition*,
				      G4double KineticEnergy,
				      G4double AtomicNumber,
				      G4double AtomicWeight=0., 
				      G4double cut = DBL_MAX,
				      G4double emax= DBL_MAX);

  void SampleScattering(const G4DynamicParticle*, G4double safety);

  void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double,
			 G4double);

  G4double ComputeTruePathLengthLimit(const G4Track& track,
				      G4PhysicsTable* theLambdaTable,
				      G4double currentMinimalStep);

  G4double ComputeGeomPathLength(G4double truePathLength);

  G4double ComputeTrueStepLength(G4double geomStepLength);

private:

  G4double ComputeXSectionPerVolume();

  void ComputeMaxElectronScattering(G4double cut);

  inline G4double GetLambda(G4double kinEnergy);

  inline void SetupParticle(const G4ParticleDefinition*);

  inline void SetupKinematic(G4double kinEnergy, G4double cut);
  
  inline void SetupTarget(G4double Z, G4double kinEnergy);

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  //  hide assignment operator
  G4MuMscModel & operator=(const  G4MuMscModel &right);
  G4MuMscModel(const  G4MuMscModel&);

  const G4ParticleDefinition* theProton;
  const G4ParticleDefinition* theElectron;
  const G4ParticleDefinition* thePositron;

  G4ParticleChangeForMSC*   fParticleChange;

  G4SafetyHelper*           safetyHelper;
  G4PhysicsTable*           theLambdaTable;
  G4PhysicsTable*           theLambda2Table;
  G4LossTableManager*       theManager;
  const G4DataVector*       currentCuts;

  G4NistManager*            fNistManager;

  G4double numlimit;
  G4double tlimitminfix;
  G4double invsqrt12;
  G4double lowBinEnergy;
  G4double highBinEnergy;

  // cash
  G4double preKinEnergy;
  G4double xSection;
  G4double ecut;
  G4double lambda0;
  G4double tPathLength;
  G4double zPathLength;
  G4double lambdaeff;
  G4double currentRange; 
  G4double par1;
  G4double par2;
  G4double par3;

  G4double xsece1;
  G4double xsece2;
  G4double xsecn2;
  G4double zcorr;
  G4double xsecn[40];
  G4double xsece[40];

  G4int    nbins;
  G4int    nwarnings;
  G4int    nwarnlimit;

  G4int    currentMaterialIndex;

  const G4MaterialCutsCouple* currentCouple;
  const G4Material* currentMaterial;

  // single scattering parameters
  G4double coeff;
  G4double constn;
  G4double cosThetaMin;
  G4double cosThetaMax;
  G4double cosTetMaxNuc;
  G4double cosTetMaxElec;
  G4double q2Limit;
  G4double alpha2;
  G4double a0;

  // projectile
  const G4ParticleDefinition* particle;

  G4double chargeSquare;
  G4double spin;
  G4double mass;
  G4double tkin;
  G4double mom2;
  G4double invbeta2;
  G4double etag;

  // target
  G4double targetZ;
  G4double screenZ;
  G4double formfactA;
  G4double FF[100];

  // flags
  G4bool   isInitialized;
  G4bool   inside;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4MuMscModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
    currentMaterial = cup->GetMaterial();
    currentMaterialIndex = currentCouple->GetIndex(); 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4MuMscModel::GetLambda(G4double e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = CrossSection(currentCouple,particle,e,
		     (*currentCuts)[currentMaterialIndex]);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline 
void G4MuMscModel::SetupParticle(const G4ParticleDefinition* p)
{
  // Initialise mass and charge
  if(p != particle) {
    particle = p;
    mass = particle->GetPDGMass();
    spin = particle->GetPDGSpin();
    G4double q = particle->GetPDGCharge()/eplus;
    chargeSquare = q*q;
    tkin = 0.0;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

inline void G4MuMscModel::SetupKinematic(G4double ekin, G4double cut)
{
  if(ekin != tkin || ecut != cut) {
    tkin  = ekin;
    mom2  = tkin*(tkin + 2.0*mass);
    invbeta2 = 1.0 +  mass*mass/mom2;
    cosTetMaxNuc = cosThetaMax;
    if(ekin <= 10.*cut && mass < MeV) {
      cosTetMaxNuc = ekin*(cosThetaMax + 1.0)/(10.*cut) - 1.0;
    }
    cosTetMaxNuc = std::max(cosTetMaxNuc, 1.0 - 0.5*q2Limit/mom2);
    ComputeMaxElectronScattering(cut);
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
  
inline void G4MuMscModel::SetupTarget(G4double Z, G4double e)
{
  if(Z != targetZ || e != etag) {
    etag    = e; 
    targetZ = Z;
    G4int iz = G4int(Z);
    if(iz > 99) iz = 99;
    G4double x = fNistManager->GetZ13(iz);
    screenZ = a0*x*x*(1.13 + 3.76*invbeta2*Z*Z*chargeSquare*alpha2)/mom2;
    //    screenZ = a0*x*x*(1.13 + 3.76*Z*Z*chargeSquare*alpha2)/mom2;
    formfactA = FF[iz];
    if(formfactA == 0.0) {
      x = fNistManager->GetA27(iz); 
      formfactA = constn*x*x;
      FF[iz] = formfactA;
    }
    formfactA *= mom2;
  } 
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

