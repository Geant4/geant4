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
// $Id: G4VEmProcess.cc,v 1.3 2003/10/14 07:50:34 vnivanch Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4VEmProcess
//
// Author:        Vladimir Ivanchenko on base of Laszlo Urban code
//
// Creation date: 01.10.2003
//
// Modifications:
//
//
// Class Description:
//
// It is the unified process for e+ annililation at rest and in fly.

// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "G4VEmProcess.hh"
#include "G4LossTableManager.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmModel.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4ProcessManager.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::G4VEmProcess(const G4String& name, G4ProcessType type):
                      G4VDiscreteProcess(name, type),
  theLambdaTable(0),
  particle(0),
  secondaryParticle(0),
  currentCouple(0),
  nLambdaBins(90),
  integral(true),
  meanFreePath(true)
{

  minKinEnergy         = 0.1*keV;
  maxKinEnergy         = 100.0*GeV;

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::~G4VEmProcess()
{
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  theLambdaTable = 0;
  modelManager->Clear();
  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::Initialise()
{
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  theLambdaTable = 0;
  modelManager->Clear();
  theCuts = modelManager->Initialise(particle,secondaryParticle,2.,verboseLevel);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if( !particle ) particle = &part;
  currentCouple = 0;
  preStepLambda = 0.0;
  if(0 < verboseLevel) {
    G4cout << "G4VEmProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }

  G4bool cutsWasModified = false;
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  for (size_t j=0; j<numOfCouples; j++){
    if (theCoupleTable->GetMaterialCutsCouple(j)->IsRecalcNeeded()) {
      cutsWasModified = true;
      break;
    }
  }
  if( !cutsWasModified ) return;

  Initialise();
  theLambdaTable = BuildLambdaTable();
  PrintInfoDefinition();

  if(0 < verboseLevel && theCuts) {
    G4cout << "G4VEmProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsTable* G4VEmProcess::BuildLambdaTable()
{

  if(0 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  G4PhysicsTable* theTable = new G4PhysicsTable(numOfCouples);

  for(size_t i=0; i<numOfCouples; i++) {

    // create physics vector and fill it
    const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
    G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
    modelManager->FillLambdaVector(aVector, couple);

    // Insert vector for this material into the table
    theTable->insert(aVector) ;
  }

  if(0 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << *theTable << G4endl;
    }
  }

  return theTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetParticle(const G4ParticleDefinition* p)
{
  particle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetSecondaryParticle(const G4ParticleDefinition* p)
{
  secondaryParticle = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::AddEmModel(G4int order, G4VEmModel* p,
                                      G4VEmFluctuationModel*,
                                const G4Region* region)
{
  modelManager->AddEmModel(order, p, 0, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::UpdateEmModel(const G4String& nam, G4double emin,
                                                              G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEmProcess::PostStepDoIt(const G4Track& track,
                                                  const G4Step& step)
{
  aParticleChange.Initialize(track);
  G4double finalT = track.GetKineticEnergy();

  // Integral approach
  if (integral) {
    G4bool b;
    G4double postStepLambda =
             (((*theLambdaTable)[currentMaterialIndex])->GetValue(finalT,b));

    if(preStepLambda*G4UniformRand() > postStepLambda)
      return G4VDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = SelectModel(finalT);
  G4double tcut = (*theCuts)[currentMaterialIndex];
  const G4DynamicParticle* dynParticle = track.GetDynamicParticle();

  /*
  if(0 < verboseLevel) {
    const G4ParticleDefinition* pd = dynParticle->GetDefinition();
    G4cout << "G4VEmProcess::PostStepDoIt: Sample secondary; E= "
           << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit(pd)
           << ", " <<  currentModel->HighEnergyLimit(pd) << ")"
           << G4endl;
  }
  */

  SecondariesPostStep(currentModel,currentCouple,dynParticle,tcut,finalT);
  return &aParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PrintInfoDefinition()
{
  G4cout << G4endl << GetProcessName() << ":  " << G4endl
         << "      Lambda tables from threshold to "
         << G4BestUnit(maxKinEnergy,"Energy")
         << " in " << nLambdaBins << " bins."
         << G4endl;
  /*
      G4cout << "DEDXTable address= " << theDEDXTable << G4endl;
      if(theDEDXTable) G4cout << (*theDEDXTable) << G4endl;
      G4cout << "RangeTable address= " << theRangeTable << G4endl;
      if(theRangeTable) G4cout << (*theRangeTable) << G4endl;
      G4cout << "InverseRangeTable address= " << theInverseRangeTable << G4endl;
      if(theInverseRangeTable) G4cout << (*theInverseRangeTable) << G4endl;
  */
  if(0 < verboseLevel) {
    G4cout << "Tables are built for " << particle->GetParticleName()
           << " IntegralFlag= " <<  integral
           << G4endl;

    if(2 < verboseLevel) {
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEmProcess::LambdaPhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4double cut  = (*theCuts)[couple->GetIndex()];
  G4int nbins = 3;
  if( couple->IsUsed() ) nbins = nLambdaBins;
  G4double tmin = std::max(MinPrimaryEnergy(particle, couple->GetMaterial(), cut),
                               minKinEnergy);
  if(tmin >= maxKinEnergy) tmin = 0.5*maxKinEnergy;
  //  G4double xmax = maxKinEnergy*exp(log(maxKinEnergy/tmin)/((G4double)(nbins-1)) );
  G4PhysicsVector* v = new G4PhysicsLogVector(tmin, maxKinEnergy, nbins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MicroscopicCrossSection(G4double kineticEnergy,
                                                 const G4MaterialCutsCouple* couple)
{
  // Cross section per atom is calculated
  DefineMaterial(couple);
  G4double cross = 0.0;
  G4bool b;
  if(theLambdaTable) {
    cross = (((*theLambdaTable)[currentMaterialIndex])->
                           GetValue(kineticEnergy, b));

    cross /= currentMaterial->GetTotNbOfAtomsPerVolume();
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MeanFreePath(const G4Track& track,
                                              G4double s,
                                              G4ForceCondition* cond)
{
  return GetMeanFreePath(track, s, cond);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEmProcess::StorePhysicsTable(G4ParticleDefinition* part,
			 	     const G4String& directory,
				           G4bool ascii)
{
  G4bool yes = true;

  if ( theLambdaTable ) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    yes = theLambdaTable->StorePhysicsTable(name,ascii);
  }

  if ( yes ) {
    G4cout << "Physics tables are stored for " << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  } else {
    G4cout << "Fail to store Physics Tables for " << particle->GetParticleName()
           << " and process " << GetProcessName()
	   << " in the directory <" << directory
	   << "> " << G4endl;
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEmProcess::RetrievePhysicsTable(G4ParticleDefinition* part,
			  	        const G4String& directory,
			  	              G4bool ascii)
{
  currentCouple = 0;
  preStepLambda = 0.0;
  if(0 < verboseLevel) {
    G4cout << "G4VEmProcess::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
	   << GetProcessName() << G4endl;
  }
  G4bool yes = true;

  const G4String particleName = part->GetParticleName();
  if( !particle ) particle = part;

  Initialise();

  G4String filename;
  const G4ProductionCutsTable* theCoupleTable=
           G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();

  filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  theLambdaTable = new G4PhysicsTable(numOfCouples);
  yes = theLambdaTable->RetrievePhysicsTable(filename,ascii);
  if ( yes ) {
      if (-1 < verboseLevel) {
        G4cout << "Lambda table for " << particleName << " is retrieved from <"
               << filename << ">"
               << G4endl;
      }
      PrintInfoDefinition();
  } else {
      theLambdaTable->clearAndDestroy();
      theLambdaTable = 0;
      if (-1 < verboseLevel) {
        G4cout << "Lambda table for " << particleName << " in file <"
               << filename << "> is not exist"
               << G4endl;
      }
  }

  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
