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
// $Id: G4VMultipleScattering.cc,v 1.21 2003/11/11 14:05:10 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VMultipleScattering
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 25.03.2003
//
// Modifications:
//
// 13.04.03 Change printout (V.Ivanchenko)
// 04-06-03 Fix compilation warnings (V.Ivanchenko)
// 16-07-03 Use G4VMscModel interface (V.Ivanchenko)
// 03-11-03 Fix initialisation problem in RetrievePhysicsTable (V.Ivanchenko)
// 04-11-03 Update PrintInfoDefinition (V.Ivanchenko)
//
//
// Class Description:
//
// It is the generic process of multiple scattering it includes common
// part of calculations for all charged particles

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VMultipleScattering.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::G4VMultipleScattering(const G4String& name, G4ProcessType type):
                 G4VContinuousDiscreteProcess(name, type),
  navigator(0),
  theLambdaTable(0),
  currentCouple(0),
  nBins(120),
  boundary(false),
  latDisplasment(true),
  buildLambdaTable(true)
{
  minKinEnergy = 100.0*eV;
  maxKinEnergy = 100.0*TeV;
  SetVerboseLevel(0);
  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::~G4VMultipleScattering()
{
  (G4LossTableManager::Instance())->DeRegister(this);
  delete modelManager;
  if (theLambdaTable) {
    theLambdaTable->clearAndDestroy();
    delete theLambdaTable;
  }
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  currentParticle = &part;
  currentCouple = 0;
  if(0 < verboseLevel) {
    G4cout << "========================================================" << G4endl;
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }

  InitialiseProcess(part);

  if(latDisplasment) navigator = G4TransportationManager::GetTransportationManager()
			       ->GetNavigatorForTracking();
  const G4DataVector* theCuts = modelManager->Initialise(&part, 0, 10.0, verboseLevel);

  if (buildLambdaTable) {

    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();
    theLambdaTable = new G4PhysicsTable(numOfCouples);

    for (size_t i=0; i<numOfCouples; i++) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = PhysicsVector(couple);
      modelManager->FillLambdaVector(aVector, couple, false);

      // Insert vector for this material into the table
      theLambdaTable->insert(aVector) ;
    }

    if(0 < verboseLevel) {
      G4cout << "Lambda table is built for "
             << part.GetParticleName()
             << G4endl;
    }
    if(2 < verboseLevel) G4cout << *theLambdaTable << G4endl;
    if(5 < verboseLevel) G4cout << theCuts << G4endl;

  }

  G4String num = part.GetParticleName();
  if (verboseLevel>0 || num == "e-" || num == "mu+" || num == "proton"
                     || num == "pi-") PrintInfoDefinition();

  if(0 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::AddEmModel(G4int order, G4VEmModel* p,
                                 const G4Region* region)
{
  G4VEmFluctuationModel* fm = 0;
  modelManager->AddEmModel(order, p, fm, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VMultipleScattering::PostStepDoIt(const G4Track& track,
                                                       const G4Step& step)
{
  fParticleChange.Initialize(track);
  G4double kineticEnergy = track.GetKineticEnergy();
  G4double truestep = step.GetStepLength();

  if (kineticEnergy > 0.0) {
    G4double cth  = currentModel->SampleCosineTheta(truestep);
    G4double sth  = sqrt(1.-cth*cth);
    G4double phi  = twopi*G4UniformRand();
    G4double dirx = sth*cos(phi);
    G4double diry = sth*sin(phi);

    G4ThreeVector oldDirection = track.GetMomentumDirection();
    G4ThreeVector newDirection(dirx,diry,cth);
    newDirection.rotateUz(oldDirection);
    fParticleChange.SetMomentumChange(newDirection);


  /*
  if(0 < verboseLevel) {
    const G4ParticleDefinition* pd = dynParticle->GetDefinition();
    G4cout << "G4VMultipleScattering::PostStepDoIt: Sample secondary; E= " << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit(pd)
           << ", " <<  currentModel->HighEnergyLimit(pd) << ")"
           << G4endl;
  }
  */

    //    G4cout << "PostStep: sth= " << sth << " trueLength= " << truestep << " tLast= " << truePathLength << G4endl;

    if (latDisplasment) {

      G4double safety = step.GetPostStepPoint()->GetSafety();
      if ( safety > 0.0) {
        G4double r = currentModel->SampleDisplacement();
        if (r > safety) r = safety;

	//    G4cout << "r= " << r << " safety= " << safety << G4endl;

        // sample direction of lateral displacement
        G4double phi  = twopi*G4UniformRand();
        G4double dirx = cos(phi);
        G4double diry = sin(phi);

        G4ThreeVector latDirection(dirx,diry,0.0);
        latDirection.rotateUz(oldDirection);

        // compute new endpoint of the Step
        G4ThreeVector newPosition = (step.GetPostStepPoint())->GetPosition()
					    + r*latDirection;

        navigator->LocateGlobalPointWithinVolume(newPosition);

        fParticleChange.SetPositionChange(newPosition);
      }
    }
  }

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::PrintInfoDefinition()
{
  G4cout << G4endl << GetProcessName() << ":  Model variant of multiple scattering " 
         << "for " << currentParticle->GetParticleName()
         << G4endl;
  if (theLambdaTable) {
    G4cout << "      Lambda tables from "
           << G4BestUnit(MinKinEnergy(),"Energy")
	   << " to "
           << G4BestUnit(MaxKinEnergy(),"Energy")
           << " in " << nBins << " bins."
           << G4endl;
  }
  if (1 < verboseLevel) {
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VMultipleScattering::PhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4int nbins = 3;
  //G4int nbins =  nDEDXBins;
  if( couple->IsUsed() ) nbins = nBins;
  //  G4double xmax = maxKinEnergy*exp( log(maxKinEnergy/minKinEnergy) / ((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nbins);
  return v;
}

G4bool G4VMultipleScattering::StorePhysicsTable(G4ParticleDefinition* part,
			 	     const G4String& directory,
				           G4bool ascii)
{
  G4bool res = true;
  if ( theLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    G4bool yes = theLambdaTable->StorePhysicsTable(name,ascii);
    if( !yes ) res = false;
  }

  if ( res ) {
    G4cout << "Physics table are stored for " << part->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  } else {
    G4cout << "Fail to store Physics Table for " << part->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VMultipleScattering::RetrievePhysicsTable(G4ParticleDefinition* part,
					     const G4String& directory,
			  	                   G4bool ascii)
{
  currentParticle = part;
  if(0 < verboseLevel) {
    G4cout << "========================================================" << G4endl;
    G4cout << "G4VMultipleScattering::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
	   << GetProcessName() << G4endl;
  }

  InitialiseProcess(*part);
  if(latDisplasment) navigator = G4TransportationManager::GetTransportationManager()
                               ->GetNavigatorForTracking();

  const G4DataVector* theCuts = modelManager->Initialise(part, 0, 10.0, verboseLevel);
  if(5 < verboseLevel) G4cout << theCuts << G4endl;
  if(!buildLambdaTable) return true;

  G4String num = part->GetParticleName();
  const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4String filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  theLambdaTable = new G4PhysicsTable(numOfCouples);
  G4bool res = theLambdaTable->RetrievePhysicsTable(filename,ascii);
  if ( res ) {
    if (-1 < verboseLevel) {
        G4cout << "Lambda table for " << num << " is retrieved from <"
               << filename << ">"
               << G4endl;
    }
  } else {
    theLambdaTable->clearAndDestroy();
    theLambdaTable = 0;
    if (-1 < verboseLevel) {
        G4cout << "Lambda table for " << num << " in file <"
               << filename << "> is not exist"
               << G4endl;
    }
  }

  if (verboseLevel>0 || num == "e-" || num == "mu+" || num == "proton")
           PrintInfoDefinition();
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

