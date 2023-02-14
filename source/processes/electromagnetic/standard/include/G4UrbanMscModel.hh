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

  //  hide assignment operator
  G4UrbanMscModel & operator=(const  G4UrbanMscModel &right) = delete;
  G4UrbanMscModel(const  G4UrbanMscModel&) = delete;

private:

  G4double SampleCosineTheta(G4double trueStepLength, G4double KineticEnergy);

  void SampleDisplacement(G4double sinTheta, G4double phi);

  void SampleDisplacementNew(G4double sinTheta, G4double phi);

  void InitialiseModelCache();

  inline void SetParticle(const G4ParticleDefinition*);

  inline G4double Randomizetlimit();
  
  inline G4double SimpleScattering();

  inline G4double ComputeStepmin();

  inline G4double ComputeTlimitmin();

  CLHEP::HepRandomEngine*     rndmEngineMod;

  const G4ParticleDefinition* particle;
  const G4ParticleDefinition* positron;
  G4ParticleChangeForMSC*     fParticleChange;
  const G4MaterialCutsCouple* couple;

  G4double mass;
  G4double charge,chargeSquare;
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

  G4double currentKinEnergy;
  G4double currentLogKinEnergy;
  G4double currentRange; 
  G4double rangeinit;
  G4double currentRadLength;

  G4double drr,finalr;

  G4double tlow;
  G4double invmev;
  G4double xmeanth = 0.0;
  G4double x2meanth = 1./3.;
  G4double rndmarray[2];

  struct mscData {
    G4double Z23, sqrtZ;
    G4double coeffth1, coeffth2;
    G4double coeffc1, coeffc2, coeffc3, coeffc4;
    G4double stepmina, stepminb;
    G4double doverra, doverrb;
    G4double posa, posb, posc, posd, pose;
  };
  static std::vector<mscData*> msc;

  // index of G4MaterialCutsCouple
  G4int idx;

  G4bool firstStep;
  G4bool insideskin;

  G4bool latDisplasmentbackup;
  G4bool dispAlg96;
  G4bool fPosiCorrection = true;
  G4bool isFirstInstance = false;
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
    chargeSquare = charge*charge;
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

inline G4double G4UrbanMscModel::SimpleScattering()
{
  // 'large angle scattering'
  // 2 model functions with correct xmean and x2mean
  G4double a = (2.*xmeanth+9.*x2meanth-3.)/(2.*xmeanth-3.*x2meanth+1.);
  G4double prob = (a+2.)*xmeanth/a;

  // sampling
  rndmEngineMod->flatArray(2, rndmarray);
  return (rndmarray[1] < prob) ? 
    -1.+2.*G4Exp(G4Log(rndmarray[0])/(a+1.)) : -1.+2.*rndmarray[0];
}

inline G4double G4UrbanMscModel::ComputeStepmin()
{
  // define stepmin using estimation of the ratio 
  // of lambda_elastic/lambda_transport
  G4double rat = currentKinEnergy*invmev;
  return lambda0*1.e-3/(2.e-3+rat*(msc[idx]->stepmina+msc[idx]->stepminb*rat));
}

inline G4double G4UrbanMscModel::ComputeTlimitmin()
{
  G4double x = (particle == positron) ? 
    0.7*msc[idx]->sqrtZ*stepmin : 0.87*msc[idx]->Z23*stepmin;
  if(currentKinEnergy < tlow) { x *= 0.5*(1.+currentKinEnergy/tlow); }
  return std::max(x, tlimitminfix);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

