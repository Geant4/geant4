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
// $Id:   $
// GEANT4 tag $Name:  $
//
// -------------------------------------------------------------------
//
//
// GEANT4 Class header file
//
//
// File name:     G4UrbanAdjointMscModel
//
// Author:        Laszlo Urban
//
// Creation date: 19.02.2013
//
// Created from G4UrbanAdjointMscModel96
//
// New parametrization for theta0
// Correction for very small step length
//
// Class Description:
//
// Implementation of the model of multiple scattering based on
// H.W.Lewis Phys Rev 78 (1950) 526 and L.Urban model

// -------------------------------------------------------------------
//

#ifndef G4UrbanAdjointMscModel_h
#define G4UrbanAdjointMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VMscModel.hh"
#include "G4MscStepLimitType.hh"
#include "G4Log.hh"
#include "G4Exp.hh"
#include "G4Electron.hh"


class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UrbanAdjointMscModel : public G4VMscModel
{

public:

  explicit G4UrbanAdjointMscModel(const G4String& nam = "UrbanMsc");

  virtual ~G4UrbanAdjointMscModel();

  virtual void Initialise(const G4ParticleDefinition*, 
			  const G4DataVector&) override;

  virtual void StartTracking(G4Track*) override;

  virtual G4double 
  ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
			     G4double KineticEnergy,
			     G4double AtomicNumber,
			     G4double AtomicWeight=0., 
			     G4double cut =0.,
			     G4double emax=DBL_MAX) override;

  virtual G4ThreeVector& SampleScattering(const G4ThreeVector&, 
					  G4double safety) override;

  virtual G4double 
  ComputeTruePathLengthLimit(const G4Track& track,
			     G4double& currentMinimalStep) override;

  virtual G4double ComputeGeomPathLength(G4double truePathLength) override;

  virtual G4double ComputeTrueStepLength(G4double geomStepLength) override;

  G4double ComputeTheta0(G4double truePathLength, G4double KineticEnergy);

  inline void SetNewDisplacementFlag(G4bool);

private:

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  void SampleDisplacement(G4double sinTheta, G4double phi);

  void SampleDisplacementNew(G4double sinTheta, G4double phi);

  inline void SetParticle(const G4ParticleDefinition*);

  inline void UpdateCache();

  inline G4double Randomizetlimit();
  
  inline G4double SimpleScattering(G4double xmeanth, G4double x2meanth);

  //  hide assignment operator
  G4UrbanAdjointMscModel & operator=(const  G4UrbanAdjointMscModel &right) = delete;
  G4UrbanAdjointMscModel(const  G4UrbanAdjointMscModel&) = delete;

  CLHEP::HepRandomEngine*     rndmEngineMod;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* positron;
  G4ParticleChangeForMSC*     fParticleChange;

  const G4MaterialCutsCouple* couple;
  G4LossTableManager*         theManager;

  G4double mass;
  G4double charge,ChargeSquare;
  G4double masslimite,lambdalimit,fr;

  G4double taubig;
  G4double tausmall;
  G4double taulim;
  G4double currentTau;
  G4double tlimit;
  G4double tlimitmin;
  G4double tlimitminfix,tlimitminfix2;
  G4double tgeom;

  G4double geombig;
  G4double geommin;
  G4double geomlimit;
  G4double skindepth;
  G4double smallstep;

  G4double presafety;

  G4double lambda0;
  G4double lambdaeff;
  G4double tPathLength;
  G4double zPathLength;
  G4double par1,par2,par3;

  G4double stepmin;

  G4double currentKinEnergy;
  G4double currentRange; 
  G4double rangeinit;
  G4double currentRadLength;

  G4int    currentMaterialIndex;

  G4double Zold;
  G4double Zeff,Z2,Z23,lnZ;
  G4double coeffth1,coeffth2;
  G4double coeffc1,coeffc2,coeffc3,coeffc4;

  G4bool   firstStep;
  G4bool   insideskin;

  G4bool   latDisplasmentbackup ;
  G4bool   displacementFlag;

  G4double rangecut;
  G4double drr,finalr;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UrbanAdjointMscModel::SetNewDisplacementFlag(G4bool val)
{
  displacementFlag = val;
}

inline
void G4UrbanAdjointMscModel::SetParticle(const G4ParticleDefinition* p)
{ const G4ParticleDefinition* p1 =p;

  if (p->GetParticleName() =="adj_e-") p1= G4Electron::Electron();

  if (p1 != particle) {
    particle = p1;
    mass = p1->GetPDGMass();
    charge = p1->GetPDGCharge()/CLHEP::eplus;
    ChargeSquare = charge*charge;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4UrbanAdjointMscModel::Randomizetlimit()
{
  G4double temptlimit = tlimit;
  if(tlimit > tlimitmin)
  {
    G4double delta = tlimit-tlimitmin;
    do {
         temptlimit = G4RandGauss::shoot(rndmEngineMod,tlimit,0.1*delta);
         // Loop checking, 10-Apr-2016, Laszlo Urban         
       } while ((temptlimit < tlimit-delta) ||
                (temptlimit > tlimit+delta));
  }
  else { temptlimit = tlimitmin; }

  return temptlimit;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UrbanAdjointMscModel::UpdateCache()
{
  lnZ = G4Log(Zeff);
  // correction in theta0 formula
  G4double w = G4Exp(lnZ/6.);
  G4double facz = 0.990395+w*(-0.168386+w*0.093286) ;
  coeffth1 = facz*(1. - 8.7780e-2/Zeff);
  coeffth2 = facz*(4.0780e-2 + 1.7315e-4*Zeff);

  // tail parameters
  G4double Z13 = w*w;
  coeffc1  = 2.3785    - Z13*(4.1981e-1 - Z13*6.3100e-2);
  coeffc2  = 4.7526e-1 + Z13*(1.7694    - Z13*3.3885e-1);
  coeffc3  = 2.3683e-1 - Z13*(1.8111    - Z13*3.2774e-1);
  coeffc4  = 1.7888e-2 + Z13*(1.9659e-2 - Z13*2.6664e-3);

  Z2   = Zeff*Zeff;
  Z23  = Z13*Z13;               
                                              
  Zold = Zeff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4UrbanAdjointMscModel::SimpleScattering(G4double xmeanth, G4double x2meanth)
{
  // 'large angle scattering'
  // 2 model functions with correct xmean and x2mean
  G4double a = (2.*xmeanth+9.*x2meanth-3.)/(2.*xmeanth-3.*x2meanth+1.);
  G4double prob = (a+2.)*xmeanth/a;

  // sampling
  G4double cth = 1.;
  if(rndmEngineMod->flat() < prob) {
    cth = -1.+2.*G4Exp(G4Log(rndmEngineMod->flat())/(a+1.));
  } else {
    cth = -1.+2.*rndmEngineMod->flat();
  }
  return cth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

