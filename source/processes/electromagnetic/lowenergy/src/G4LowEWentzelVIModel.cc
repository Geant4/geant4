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
// $Id: G4LowEWentzelVIModel.cc 74528 2013-10-12 17:24:24Z vnivanch $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:   G4LowEWentzelVIModel
//
// Author:      V.Ivanchenko 
//
// Creation date: 11.02.2014 from G4WentzelVIModel
//
// Modifications:
//
// Class Description:
//

// -------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4LowEWentzelVIModel.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

using namespace std;

G4LowEWentzelVIModel::G4LowEWentzelVIModel() :
  G4WentzelVIModel(false,"LowEnWentzelVI")
{
  SetSingleScatteringFactor(0.5);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LowEWentzelVIModel::~G4LowEWentzelVIModel()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4LowEWentzelVIModel::ComputeTruePathLengthLimit(
                             const G4Track& track,
			     G4double& currentMinimalStep)
{
  G4double tlimit = currentMinimalStep;
  const G4DynamicParticle* dp = track.GetDynamicParticle();
  G4StepPoint* sp = track.GetStep()->GetPreStepPoint();
  G4StepStatus stepStatus = sp->GetStepStatus();
  singleScatteringMode = false;
  //G4cout << "G4LowEWentzelVIModel::ComputeTruePathLengthLimit stepStatus= " 
  //	 << stepStatus << "  " << track.GetDefinition()->GetParticleName() 
  //	 << G4endl;

  // initialisation for each step, lambda may be computed from scratch
  preKinEnergy  = dp->GetKineticEnergy();
  DefineMaterial(track.GetMaterialCutsCouple());
  lambdaeff = GetTransportMeanFreePath(particle,preKinEnergy);
  currentRange = GetRange(particle,preKinEnergy,currentCouple);
  cosTetMaxNuc = wokvi->SetupKinematic(preKinEnergy, currentMaterial);
  /*
  G4cout << "lambdaeff= " << lambdaeff << " Range= " << currentRange
  	 << " tlimit= " << tlimit << " 1-cost= " << 1 - cosTetMaxNuc << G4endl;
  */
  // extra check for abnormal situation
  // this check needed to run MSC with eIoni and eBrem inactivated
  tlimit = std::min(tlimit, currentRange);

  // stop here if small range particle
  if(tlimit < tlimitminfix) { 
    return ConvertTrueToGeom(tlimit, currentMinimalStep); 
  }

  // pre step
  G4double presafety = sp->GetSafety();
  // far from geometry boundary
  if(currentRange < presafety) {
    return ConvertTrueToGeom(tlimit, currentMinimalStep);
  }

  // compute presafety again if presafety <= 0 and no boundary
  // i.e. when it is needed for optimization purposes
  if(stepStatus != fGeomBoundary && presafety < tlimitminfix) {
    presafety = ComputeSafety(sp->GetPosition(), tlimit); 
    if(currentRange < presafety) {
      return ConvertTrueToGeom(tlimit, currentMinimalStep);
    }
  }
  /*   
  G4cout << "e(MeV)= " << preKinEnergy/MeV
	 << "  " << particle->GetParticleName() 
	 << " CurLimit(mm)= " << tlimit/mm <<" safety(mm)= " << presafety/mm
	 << " R(mm)= " <<currentRange/mm
	 << " L0(mm^-1)= " << lambdaeff*mm 
	 <<G4endl;
  */
  // natural limit for high energy
  G4double rlimit = std::max(facrange*currentRange, lambdaeff);
  //G4double rlimit = std::max(facrange*currentRange, 
  //			     0.7*(1.0 - cosTetMaxNuc)*lambdaeff);

  // low-energy e-
  rlimit = std::max(rlimit, facsafety*presafety);
   
  // cut correction
  //G4double rcut = currentCouple->GetProductionCuts()->GetProductionCut(1);
  //G4cout << "rcut= " << rcut << " rlimit= " << rlimit << " presafety= " 
  // << presafety << " 1-cosThetaMax= " <<1-cosThetaMax 
  //<< " 1-cosTetMaxNuc= " << 1-cosTetMaxNuc << G4endl;
  //if(rcut > rlimit) { rlimit = std::min(rlimit, rcut*sqrt(rlimit/rcut)); }
 
  tlimit = std::min(tlimit, rlimit);
  tlimit = std::max(tlimit, tlimitminfix);

  // step limit in infinite media
  tlimit = std::min(tlimit, 50*currentMaterial->GetRadlen()/facgeom);

  //compute geomlimit and force few steps within a volume
  if (steppingAlgorithm == fUseDistanceToBoundary 
      && stepStatus == fGeomBoundary) {

    G4double geomlimit = ComputeGeomLimit(track, presafety, currentRange);
    tlimit = std::min(tlimit, geomlimit/facgeom);
  } 
  /*  
  G4cout << particle->GetParticleName() << " E(MeV)= " << preKinEnergy
	 << " L0= " << lambdaeff << " R= " << currentRange
	 << " tlimit= " << tlimit  
  	 << " currentMinimalStep= " << currentMinimalStep << G4endl;
  */
  return ConvertTrueToGeom(tlimit, currentMinimalStep);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
