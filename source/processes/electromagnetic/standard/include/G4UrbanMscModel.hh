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
// Creation date: 19.02.2013
//
// Created from G4UrbanMscModel96 
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

#ifndef G4UrbanMscModel_h
#define G4UrbanMscModel_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <CLHEP/Units/SystemOfUnits.h>

#include "G4VMscModel.hh"
#include "G4MscStepLimitType.hh"
#include "G4Log.hh"
#include "G4Exp.hh"

class G4ParticleChangeForMSC;
class G4SafetyHelper;
class G4LossTableManager;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4UrbanMscModel : public G4VMscModel
{

public:

  explicit G4UrbanMscModel(const G4String& nam = "UrbanMsc");

  ~G4UrbanMscModel() override;

  void Initialise(const G4ParticleDefinition*, 
		  const G4DataVector&) override;

  void StartTracking(G4Track*) override;

  G4double 
  ComputeCrossSectionPerAtom(const G4ParticleDefinition* particle,
			     G4double KineticEnergy,
			     G4double AtomicNumber,
			     G4double AtomicWeight=0., 
			     G4double cut =0.,
			     G4double emax=DBL_MAX) override;

  G4ThreeVector& SampleScattering(const G4ThreeVector&, 
				  G4double safety) override;

  G4double ComputeTruePathLengthLimit(const G4Track& track,
			              G4double& currentMinimalStep) override;

  G4double ComputeGeomPathLength(G4double truePathLength) override;

  G4double ComputeTrueStepLength(G4double geomStepLength) override;

  G4double ComputeTheta0(G4double truePathLength, G4double KineticEnergy);

private:

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  void SampleDisplacement(G4double sinTheta, G4double phi);

  void SampleDisplacementNew(G4double sinTheta, G4double phi);

  inline void SetParticle(const G4ParticleDefinition*);

  inline void UpdateCache();

  inline G4double Randomizetlimit();
  
  inline G4double SimpleScattering(G4double xmeanth, G4double x2meanth);

  //  hide assignment operator
  G4UrbanMscModel & operator=(const  G4UrbanMscModel &right) = delete;
  G4UrbanMscModel(const  G4UrbanMscModel&) = delete;

  CLHEP::HepRandomEngine*     rndmEngineMod;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* positron;
  G4ParticleChangeForMSC*     fParticleChange;

  const G4MaterialCutsCouple* couple;
  G4LossTableManager*         theManager;

  G4double mass;
  G4double charge,ChargeSquare;
  G4double masslimite,fr;

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
  G4double stepmina,stepminb;

  G4double currentKinEnergy;
  G4double currentLogKinEnergy;
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

  G4bool   latDisplasmentbackup;
  G4bool   dispAlg96;

  G4double rangecut;
  G4double drr,finalr;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
void G4UrbanMscModel::SetParticle(const G4ParticleDefinition* p)
{
  if (p != particle) {
    particle = p;
    mass = p->GetPDGMass();
    charge = p->GetPDGCharge()/CLHEP::eplus;
    ChargeSquare = charge*charge;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4UrbanMscModel::Randomizetlimit()
{
  G4double res = tlimitmin;
  if(tlimit > tlimitmin)
  {
    res = G4RandGauss::shoot(rndmEngineMod,tlimit,0.1*(tlimit-tlimitmin));
    res = std::max(res, tlimitmin);
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void G4UrbanMscModel::UpdateCache()                                   
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

  stepmina = 15.99/(1.+0.119*Zeff);
  stepminb =  4.39/(1.+0.079*Zeff);
                                              
  Zold = Zeff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4double G4UrbanMscModel::SimpleScattering(G4double xmeanth, G4double x2meanth)
{
  // 'large angle scattering'
  // 2 model functions with correct xmean and x2mean
  G4double a = (2.*xmeanth+9.*x2meanth-3.)/(2.*xmeanth-3.*x2meanth+1.);
  G4double prob = (a+2.)*xmeanth/a;

  // sampling
  G4double rdm = rndmEngineMod->flat();
  G4double cth = (rndmEngineMod->flat() < prob)
    ? -1.+2.*G4Exp(G4Log(rdm)/(a+1.)) : -1.+2.*rdm;
  return cth;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

