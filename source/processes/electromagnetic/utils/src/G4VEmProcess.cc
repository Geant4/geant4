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
// $Id: G4VEmProcess.cc,v 1.9 2004-08-08 12:09:44 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 30-06-04 make it to be pure discrete process (V.Ivanchenko)
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
#include "G4DataVector.hh"
#include "G4PhysicsTable.hh"
#include "G4PhysicsVector.hh"
#include "G4PhysicsLogVector.hh"
#include "G4VParticleChange.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::G4VEmProcess(const G4String& name, G4ProcessType type):
                      G4VDiscreteProcess(name, type),
  theLambdaTable(0),
  theEnergyOfCrossSectionMax(0),
  theCrossSectionMax(0),
  particle(0),
  secondaryParticle(0),
  nLambdaBins(90),
  lambdaFactor(0.1),
  currentCouple(0),
  integral(false),
  meanFreePath(true),
  aboveCSmax(true)
{

  minKinEnergy    = 0.1*keV;
  maxKinEnergy    = 100.0*GeV;
  mfpKinEnergy    = DBL_MAX;

  pParticleChange = &fParticleChange;

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::~G4VEmProcess()
{
  Clear();
  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::Initialise()
{
  Clear();
  const std::vector<G4double>* theCutsGamma =
        modelManager->Initialise(particle,secondaryParticle,2.,verboseLevel);
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  theCutsGamma    = theCoupleTable->GetEnergyCutsVector(0);
  theCutsElectron = theCoupleTable->GetEnergyCutsVector(1);
  theCutsPositron = theCoupleTable->GetEnergyCutsVector(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::Clear()
{
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  if(theEnergyOfCrossSectionMax) delete [] theEnergyOfCrossSectionMax;
  if(theCrossSectionMax) delete [] theCrossSectionMax;
  theLambdaTable = 0;
  theEnergyOfCrossSectionMax = 0;
  theCrossSectionMax = 0;
  modelManager->Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  currentCouple = 0;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;

  if( !particle ) particle = &part;

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
  FindLambdaMax();
  PrintInfoDefinition();

  if(0 < verboseLevel) {
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

void G4VEmProcess::AddEmModel(G4int order, G4VEmModel* p, const G4Region* region)
{
  modelManager->AddEmModel(order, p, 0, region);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::UpdateEmModel(const G4String& nam, G4double emin, G4double emax)
{
  modelManager->UpdateEmModel(nam, emin, emax);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4VEmProcess::PostStepDoIt(const G4Track& track,
                                              const G4Step& step)
{
  fParticleChange.Initialize(track);

  // Do not make anything if particle is stopped, the annihilation then
  // should be performed by the AtRestDoIt!
  if (track.GetTrackStatus() == fStopButAlive) return &fParticleChange;

  G4double finalT = track.GetKineticEnergy();

  // Integral approach
  if (integral) {
    G4double lx = GetLambda(finalT);
    if(preStepLambda<lx && 0 < verboseLevel) {
      G4cout << "WARING: for " << particle->GetParticleName() 
             << " and " << GetProcessName() 
             << " E(MeV)= " << finalT/MeV
             << " preLambda= " << preStepLambda << " < " << lx << " (postLambda) "
	     << G4endl;  
    }

    if(preStepLambda*G4UniformRand() > lx)
      return G4VDiscreteProcess::PostStepDoIt(track,step);
  }

  G4VEmModel* currentModel = SelectModel(finalT);
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

  std::vector<G4DynamicParticle*>* newp = SecondariesPostStep(currentModel,
                                                              currentCouple,
		                                              dynParticle);
  if (newp) {
    G4int num = newp->size();
    if(num > 0) {
      fParticleChange.SetNumberOfSecondaries(num);
      G4double gcut = (*theCutsGamma)[currentMaterialIndex];
      G4double ecut = (*theCutsElectron)[currentMaterialIndex];
      G4double pcut = (*theCutsPositron)[currentMaterialIndex];
      G4double edep = 0.0;
      for (G4int i=0; i<num; i++) {
        G4DynamicParticle* dp = (*newp)[i];
        const G4ParticleDefinition* p = dp->GetDefinition();
        G4double e = dp->GetKineticEnergy();
	G4bool good = true;
	if (p == G4Gamma::Gamma()) {
	   if (e < gcut) {
	     good = false;
	     edep += e;
	   }
	} else if (p == G4Electron::Electron()) {
	   if (e < ecut) {
	     good = false;
	     edep += e;
	   }
	} else if (p == G4Positron::Positron()) {
	   if (e < pcut) {
	     good = false;
	     edep += e + electron_mass_c2;
	   }
        }
        if (good) aParticleChange.AddSecondary(dp);
        else      delete dp;
      }
      fParticleChange.SetLocalEnergyDeposit(edep);
    }
    delete newp;
  }

  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PrintInfoDefinition()
{
  G4cout << G4endl << GetProcessName() << ":  " << G4endl
         << "      Lambda tables from threshold to "
         << G4BestUnit(maxKinEnergy,"Energy")
         << " in " << nLambdaBins << " bins."
         << G4endl;

  if(0 < verboseLevel) {
    G4cout << "Tables are built for " << particle->GetParticleName()
           << G4endl;

    if(2 < verboseLevel) {
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
    }
  }
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
  mfpKinEnergy  = DBL_MAX;

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
      FindLambdaMax();
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

void G4VEmProcess::FindLambdaMax()
{
  if(1 < verboseLevel) {
    G4cout << "### FindLambdaMax: " << particle->GetParticleName() 
           << " and process " << GetProcessName() << G4endl; 
  }
  size_t n = theLambdaTable->length();
  G4PhysicsVector* pv = (*theLambdaTable)[0];
  size_t nb = pv->GetVectorLength();
  G4double emax = pv->GetLowEdgeEnergy(nb);
  G4double e, s, smax = 0.0;
  theEnergyOfCrossSectionMax = new G4double [n];
  theCrossSectionMax = new G4double [n];
  G4bool b;

  for (size_t i=0; i<n; i++) {
    pv = (*theLambdaTable)[i];
    smax = 0.0;
    for (size_t j=0; j<nb; j++) {
      e = pv->GetLowEdgeEnergy(j);
      s = pv->GetValue(e,b);
      if(s > smax) {
	smax = s;
	emax = e;
      }
    }
    theEnergyOfCrossSectionMax[i] = emax;
    theCrossSectionMax[i] = smax;
    if(2 < verboseLevel) {
      G4cout << "For " << particle->GetParticleName() 
	     << " Max CS at i= " << i << " emax(MeV)= " << emax/MeV
	     << " lambda= " << smax << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetLambdaBinning(G4int nbins)
{
  nLambdaBins = nbins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4int G4VEmProcess::LambdaBinning() const
{
  return nLambdaBins;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetMinKinEnergy(G4double e)
{
  minKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MinKinEnergy() const
{
  return minKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetMaxKinEnergy(G4double e)
{
  maxKinEnergy = e;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::MaxKinEnergy() const
{
  return maxKinEnergy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::ActivateFluorescence(G4bool, const G4Region*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::ActivateAugerElectronProduction(G4bool, const G4Region*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

const G4PhysicsTable* G4VEmProcess::LambdaTable() const
{
  return theLambdaTable;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetIntegral(G4bool val)
{
  integral = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEmProcess::IsIntegral() const
{
  return integral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
