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
// $Id: G4UrbanMscModel.hh,v 1.31 2007/10/29 08:42:43 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-01 $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4UrbanMscModel
//
// Author:        Laszlo Urban
//
// Creation date: 01.03.2001
//
// Modifications:
//
// 27-03-03 Move model part from G4MultipleScattering (V.Ivanchenko)
// 27-03-03 Rename (V.Ivanchenko)
//
// 05-08-03 angle distribution has been modified (L.Urban)
// 26-11-03 new data member currentRange (L.Urban)
// 01-03-04 changes in data members + signature changed in SampleCosineTheta
// 11-03-04 changes in data members (L.Urban)
// 23-04-04 changes in data members and in signature of SampleCosineTheta
//          (L.Urban)
// 17-08-04 name of data member facxsi changed to factail (L.Urban)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 15-04-05 optimize internal interface - add SampleSecondaries method (V.Ivanchenko)
// 11-08-05 computation of lateral correlation added (L.Urban)
// 02-10-05 nuclear size correction computation removed, the correction
//          included in the (theoretical) tabulated values (L.Urban)
// 16-02-06 data members b and xsi have been removed (L.Urban)
// 17-02-06 Save table of transport cross sections not mfp (V.Ivanchenko)
// 07-03-06 Create G4UrbanMscModel and move there step limit calculation (V.Ivanchenko)
// 10-05-06 SetMscStepLimitation at initialisation (V.Ivantchenko)
// 11-05-06 name of data member safety changed to presafety, some new data
//          members added (frscaling1,frscaling2,tlimitminfix,nstepmax)
//         (L.Urban)
// 13-10-06 data member factail removed, data member tkinlimit changed
//          to lambdalimit,
//          new data members tgeom,tnow,skin,skindepth,Zeff,geomlimit
//          G4double GeomLimit(const G4Track& track) changed to
//              void GeomLimit(const G4Track& track) (L.Urban)
// 20-10-06 parameter theta0 now computed in the (public)
//          function ComputeTheta0,
//          single scattering modified allowing not small
//          angles as well (L.Urban)
// 31-01-07 code cleaning (L.Urban)
// 06-02-07 Move SetMscStepLimitation method into the source (VI)
// 15-02-07 new data member : smallstep (L.Urban)
// 10-04-07 remove navigator, smallstep, tnow (V.Ivanchenko)
//
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4UrbanMscModel_h
#define G4UrbanMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4VEmModel.hh"
#include "G4PhysicsTable.hh"
#include "G4MscStepLimitType.hh"

class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UrbanMscModel : public G4VEmModel
{

public:

  G4UrbanMscModel(G4double facrange, G4double dtrl, G4double lambdalimit, 
		  G4double facgeom,G4double skin, 
		  G4bool samplez, G4MscStepLimitType stepAlg, 
		  const G4String& nam = "UrbanMscUni");

  virtual ~G4UrbanMscModel();

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
  G4UrbanMscModel & operator=(const  G4UrbanMscModel &right);
  G4UrbanMscModel(const  G4UrbanMscModel&);

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
void G4UrbanMscModel::SetLateralDisplasmentFlag(G4bool val) 
{ 
  latDisplasment = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetSkin(G4double val) 
{ 
  skin = val;
  stepmin       = tlimitminfix;
  skindepth     = skin*stepmin;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetRangeFactor(G4double val) 
{ 
  facrange = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetGeomFactor(G4double val) 
{ 
  facgeom = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetStepLimitType(G4MscStepLimitType val) 
{ 
  steppingAlgorithm = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4UrbanMscModel::GetLambda(G4double e)
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
void G4UrbanMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/eplus;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

