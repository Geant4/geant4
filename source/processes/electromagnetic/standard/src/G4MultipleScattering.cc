//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4MultipleScattering.cc,v 1.27 2005/04/15 14:41:13 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// -----------------------------------------------------------------------------
// 16/05/01 value of cparm changed , L.Urban
// 18/05/01 V.Ivanchenko Clean up against Linux ANSI compilation
// 07/08/01 new methods Store/Retrieve PhysicsTable (mma)
// 23-08-01 new angle and z distribution,energy dependence reduced,
//          Store,Retrieve methods commented out temporarily, L.Urban
// 27-08-01 in BuildPhysicsTable:aParticleType.GetParticleName()=="mu+" (mma)
// 28-08-01 GetContinuousStepLimit and AlongStepDoIt moved from .icc file (mma)
// 03-09-01 value of data member factlim changed, L.Urban
// 10-09-01 small change in GetContinuousStepLimit, L.Urban
// 11-09-01 G4MultipleScatteringx put as default G4MultipleScattering
//          store/retrieve physics table reactivated (mma)
// 13-09-01 corr. in ComputeTransportCrossSection, L.Urban
// 14-09-01 protection in GetContinuousStepLimit, L.Urban
// 17-09-01 migration of Materials to pure STL (mma)
// 27-09-01 value of data member factlim changed, L.Urban
// 31-10-01 big fixed in PostStepDoIt,L.Urban
// 17-04-02 NEW angle distribution + boundary algorithm modified, L.Urban
// 22-04-02 boundary algorithm modified -> important improvement in timing (L.Urban)
// 24-04-02 some minor changes in boundary algorithm, L.Urban
// 06-05-02 bug fixed in GetContinuousStepLimit, L.Urban
// 24-05-02 changes in angle distribution and boundary algorithm, L.Urban
// 11-06-02 bug fixed in ComputeTransportCrossSection, L.Urban
// 12-08-02 bug fixed in PostStepDoIt (lateral displacement), L.Urban
// 15-08-02 new angle distribution, L.Urban
// 26-09-02 angle distribution + boundary algorithm modified, L.Urban
// 15-10-02 temporary fix for proton scattering
// 30-10-02 modified angle distribution,mods in boundary algorithm,
//          changes in data members, L.Urban
// 11-12-02 precision problem in ComputeTransportCrossSection
//          for small Tkin/for heavy particles cured from L.Urban
// 20-01-03 Migrade to cut per region (V.Ivanchenko)
// 05-02-03 changes in data members, new sampling for geom.
//          path length, step dependence reduced with new
//          method
// 28-03-03 Move to model design (V.Ivanchenko)
// 08-08-03 STD substitute standard  (V.Ivanchenko)
// 23-04-04 value of data member dtrl changed from 0.15 to 0.05 (L.Urban)
// 17-08-04 name of facxsi changed to factail (L.Urban)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 07-02-05 correction in order to have a working Setsamplez function (L.Urban)
// 15-04-05 optimize internal interface (V.Ivanchenko)
// -----------------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4MultipleScattering.hh"
#include "G4MscModel.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4MultipleScattering::G4MultipleScattering(const G4String& processName)
  : G4VMultipleScattering(processName),
    totBins(120),
    facrange(0.199),
    dtrl(0.05),
    NuclCorrPar (0.0615),
    FactPar(0.40),
    factail(1.0),
    cf(1.001),
    stepnolastmsc(-1000000),
    nsmallstep(5),
    samplez(true),
    boundary(true),
    isInitialized(false)
{
  lowKineticEnergy = 0.1*keV;
  highKineticEnergy= 100.*TeV;

  tlimit           = 1.e10*mm;
  tlimitmin        = 1.e-7*mm;

  SetBinning(totBins);
  SetMinKinEnergy(lowKineticEnergy);
  SetMaxKinEnergy(highKineticEnergy);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4MultipleScattering::~G4MultipleScattering()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering::InitialiseProcess(const G4ParticleDefinition* particle)
{
  if(isInitialized) return;

  if (particle->GetParticleType() == "nucleus") {
    boundary = false;
    SetLateralDisplasmentFlag(false);
    SetBuildLambdaTable(false);
    Setsamplez(false) ;
  } else {
    SetLateralDisplasmentFlag(true);
    SetBuildLambdaTable(true);
  }
  G4MscModel* em = new G4MscModel(dtrl,NuclCorrPar,FactPar,factail,samplez);
  em->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
  em->SetLowEnergyLimit(lowKineticEnergy);
  em->SetHighEnergyLimit(highKineticEnergy);
  AddEmModel(1, em);
  isInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4MultipleScattering::TruePathLengthLimit(const G4Track&  track,
                                                            G4double& lambda,
                                                            G4double  currentMinimalStep)
{

  G4double tPathLength = currentMinimalStep;

  // special treatment near boundaries ?
  if (boundary) {

    G4int stepno = track.GetCurrentStepNumber() ;
    // first step
    if (stepno == 1) {
      stepnolastmsc = -1000000 ;
      tlimit = 1.e10;
    } else if (stepno > 1) {

      if (track.GetStep()->GetPreStepPoint()->GetStepStatus() == fGeomBoundary) {

        stepnolastmsc = stepno;
        //  if : diff.treatment for small/not small Z
	G4double range = CurrentRange();
        if (range > lambda) tlimit = facrange*range;
        else                tlimit = facrange*lambda;

        if(tlimit < tlimitmin) tlimit = tlimitmin;
        if(tPathLength > tlimit) tPathLength = tlimit;

      } else if (stepno > stepnolastmsc && stepno - stepnolastmsc < nsmallstep
              && tPathLength > tlimit) {
        tlimit *= cf;
        tPathLength = tlimit;
      }
    }
  }

  return tPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4MultipleScattering::PrintInfo()
{
  if(boundary) {
    G4cout << "      Boundary algorithm is active with facrange= "
           << facrange
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

