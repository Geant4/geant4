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
// $Id: G4UrbanMscModel90.hh,v 1.1 2007/12/07 17:35:52 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4UrbanMscModel90
//
// Author:        V.Ivanchenko clone Laszlo Urban model
//
// Creation date: 07.12.2007
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

#ifndef G4UrbanMscModel90_h
#define G4UrbanMscModel90_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"

class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UrbanMscModel90 : public G4VEmModel
{

public:

  G4UrbanMscModel90(G4double facrange, G4double dtrl, G4double lambdalimit, 
		  G4double facgeom,G4double skin, 
		  G4bool samplez, G4MscStepLimitType stepAlg, 
		  const G4String& nam = "UrbanMscUni");

  virtual ~G4UrbanMscModel90();

  virtual void Initialise(const G4ParticleDefinition*, const G4DataVector&);

  G4double ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
				      G4double KineticEnergy,
				      G4double AtomicNumber,
				      G4double AtomicWeight=0., 
				      G4double cut =0.,
				      G4double emax=DBL_MAX);

  void SampleSecondaries(std::vector<G4DynamicParticle*>*, 
			 const G4MaterialCutsCouple*,
			 const G4DynamicParticle*,
			 G4double,
			 G4double);

  void SampleScattering(const G4DynamicParticle*,
			G4double safety);

  G4double ComputeTruePathLengthLimit(const G4Track& track,
				      G4PhysicsTable* theLambdaTable,
				      G4double currentMinimalStep);

  G4double ComputeGeomPathLength(G4double truePathLength);

  G4double ComputeTrueStepLength(G4double geomStepLength);

  G4double ComputeTheta0(G4double truePathLength,
                         G4double KineticEnergy);

  void SetStepLimitType(G4MscStepLimitType);

  void SetLateralDisplasmentFlag(G4bool val);

  void SetRangeFactor(G4double);

  void SetGeomFactor(G4double);

  void SetSkin(G4double);

private:

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  G4double SampleDisplacement();

  G4double LatCorrelation();

  G4double GetLambda(G4double kinEnergy);

  void GeomLimit(const G4Track& track);

  void SetParticle(const G4ParticleDefinition* p);

  //  hide assignment operator
  G4UrbanMscModel90 & operator=(const  G4UrbanMscModel90 &right);
  G4UrbanMscModel90(const  G4UrbanMscModel90&);

  const G4ParticleDefinition* particle;
  G4ParticleChangeForMSC*     fParticleChange;

  G4SafetyHelper*             safetyHelper;
  G4PhysicsTable*             theLambdaTable;
  const G4MaterialCutsCouple* couple;
  G4LossTableManager*         theManager;


  G4double mass;
  G4double charge;

  G4double masslimite,masslimitmu;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double dtrl;

  G4double lambdalimit;
  G4double facrange;
  G4double frscaling1,frscaling2;
  G4double tlimit;
  G4double tlimitmin;
  G4double tlimitminfix;

  G4double nstepmax;
  G4double geombig;
  G4double geommin;
  G4double geomlimit;
  G4double facgeom;
  G4double skin;
  G4double skindepth;
  G4double smallstep;

  G4double presafety;
  G4double facsafety;

  G4double lambda0;
  G4double lambdaeff;
  G4double tPathLength;
  G4double zPathLength;
  G4double par1,par2,par3 ;

  G4double stepmin ;

  G4double currentKinEnergy;
  G4double currentRange; 
  G4double currentRadLength;

  G4double Zeff;

  G4int    currentMaterialIndex;

  G4MscStepLimitType steppingAlgorithm;

  G4bool   samplez;
  G4bool   latDisplasment;
  G4bool   isInitialized;

  G4bool   inside;
  G4bool   insideskin;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetLateralDisplasmentFlag(G4bool val) 
{ 
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetSkin(G4double val) 
{ 
  skin = val;
  stepmin       = tlimitminfix;
  skindepth     = skin*stepmin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetRangeFactor(G4double val) 
{ 
  facrange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetGeomFactor(G4double val) 
{ 
  facgeom = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetStepLimitType(G4MscStepLimitType val) 
{ 
  steppingAlgorithm = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4UrbanMscModel90::GetLambda(G4double e)
{
  G4double x;
  if(theLambdaTable) {
    G4bool b;
    x = ((*theLambdaTable)[currentMaterialIndex])->GetValue(e, b);
  } else {
    x = CrossSection(couple,particle,e);
  }
  if(x > DBL_MIN) x = 1./x;
  else            x = DBL_MAX;
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel90::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/eplus;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

