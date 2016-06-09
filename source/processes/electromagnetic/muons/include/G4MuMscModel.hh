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
// $Id: G4MuMscModel.hh,v 1.4 2007/11/09 19:48:09 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
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

// -------------------------------------------------------------------
//

#ifndef G4MuMscModel_h
#define G4MuMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4eCoulombScatteringModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"
#include "G4MaterialCutsCouple.hh"

class G4LossTableManager;
class G4ParticleChangeForMSC;
class G4SafetyHelper;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4MuMscModel : public G4eCoulombScatteringModel
{

public:

  G4MuMscModel(G4double frange = 0.2,
	       G4double thetaMax = 0.04,
	       G4double tMax = TeV*TeV, 
	       const G4String& nam = "MuMscUni");

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

  inline void SetStepLimitType(G4MscStepLimitType);

  inline void SetLateralDisplasmentFlag(G4bool val);

  inline G4double GetLambda(G4double kinEnergy);

  inline G4double GetLambda2(G4double kinEnergy);

  //  inline void SetThetaLimit(G4double);

  inline void SetRangeFactor(G4double);

private:

  void BuildTables();

  G4double ComputeLambda2(G4double kinEnergy, G4double cut);

  inline void DefineMaterial(const G4MaterialCutsCouple*);

  //  hide assignment operator
  G4MuMscModel & operator=(const  G4MuMscModel &right);
  G4MuMscModel(const  G4MuMscModel&);

  G4ParticleChangeForMSC*     fParticleChange;

  G4SafetyHelper*             safetyHelper;
  G4PhysicsTable*             theLambdaTable;
  G4PhysicsTable*             theLambda2Table;
  G4LossTableManager*         theManager;
  const G4DataVector*         currentCuts;

  G4double dtrl;
  G4double facrange;
  G4double thetaLimit;
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

  G4int    currentMaterialIndex;

  G4int    nbins;
  G4int    nwarnings;
  G4int    nwarnlimit;

  const G4MaterialCutsCouple* currentCouple;

  G4MscStepLimitType steppingAlgorithm;

  G4bool   samplez;
  G4bool   latDisplasment;
  G4bool   isInitialized;
  G4bool   buildTables;
  G4bool   newrun;
  G4bool   inside;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4MuMscModel::SetLateralDisplasmentFlag(G4bool val) 
{ 
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
inline
void G4MuMscModel::SetThetaLimit(G4double val) 
{ 
  thetaLimit = val;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4MuMscModel::SetRangeFactor(G4double val) 
{ 
  facrange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4MuMscModel::SetStepLimitType(G4MscStepLimitType val) 
{ 
  steppingAlgorithm = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4MuMscModel::DefineMaterial(const G4MaterialCutsCouple* cup) 
{ 
  if(cup != currentCouple) {
    currentCouple = cup;
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
    x = CrossSection(currentCouple,particle,e);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4MuMscModel::GetLambda2(G4double e)
{
  G4double x;
  if(theLambda2Table) {
    G4bool b;
    x = ((*theLambda2Table)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = ComputeLambda2(e, (*currentCuts)[currentMaterialIndex]);
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

