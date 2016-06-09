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
// $Id: G4VEmProcess.cc,v 1.30 2005/05/30 08:22:38 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-01 $
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
// 30-09-08 optimise integral option (V.Ivanchenko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivanchenko)
// 11-03-05 Shift verbose level by 1, add applyCuts and killPrimary flags (V.Ivanchenko)
// 14-03-05 Update logic PostStepDoIt (V.Ivanchenko)
// 08-04-05 Major optimisation of internal interfaces (V.Ivantchenko)
// 18-04-05 Use G4ParticleChangeForGamma (V.Ivantchenko)
//
//
// Class Description:
//

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
#include "G4PhysicsTableHelper.hh"
#include "G4NistManager.hh"

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
  aboveCSmax(true),
  buildLambdaTable(true),
  applyCuts(false),
  startFromNull(true),
  nRegions(0)
{
  SetVerboseLevel(1);
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 100.0*GeV;
  theGamma     = G4Gamma::Gamma();
  theElectron  = G4Electron::Electron();
  thePositron  = G4Positron::Positron();

  pParticleChange = &fParticleChange;

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmProcess::~G4VEmProcess()
{
  Clear();
  if(theLambdaTable) theLambdaTable->clearAndDestroy();
  delete modelManager;
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(!particle) particle = &part;
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::PreparePhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
	   << " local particle " << particle->GetParticleName() 
           << G4endl;
  }

  if(particle == &part) {
    Clear();
    InitialiseProcess(particle);
    theCutsGamma =
        modelManager->Initialise(particle,secondaryParticle,2.,verboseLevel);
    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    theCutsGamma    = theCoupleTable->GetEnergyCutsVector(idxG4GammaCut);
    theCutsElectron = theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut);
    theCutsPositron = theCoupleTable->GetEnergyCutsVector(idxG4PositronCut);
    if(buildLambdaTable)
      theLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::Clear()
{
  if(theEnergyOfCrossSectionMax) delete [] theEnergyOfCrossSectionMax;
  if(theCrossSectionMax) delete [] theCrossSectionMax;
  theEnergyOfCrossSectionMax = 0;
  theCrossSectionMax = 0;
  modelManager->Clear();
  currentCouple = 0;
  preStepLambda = 0.0;
  mfpKinEnergy  = DBL_MAX;
  preStepMFP    = DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
	   << " buildLambdaTable= " << buildLambdaTable
           << G4endl;
  }

  if(buildLambdaTable) {
    BuildLambdaTable();
    FindLambdaMax();
  }
  if(0 < verboseLevel) PrintInfoDefinition();

  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::BuildLambdaTable()
{
  if(1 < verboseLevel) {
    G4cout << "G4VEnergyLossSTD::BuildLambdaTable() for process "
           << GetProcessName() << " and particle "
           << particle->GetParticleName()
           << G4endl;
  }

  // Access to materials
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  for(size_t i=0; i<numOfCouples; i++) {

    if (theLambdaTable->GetFlag(i)) {

      // create physics vector and fill it
      const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
      G4PhysicsVector* aVector = LambdaPhysicsVector(couple);
      modelManager->FillLambdaVector(aVector, couple, startFromNull);
      G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
    }
  }

  if(1 < verboseLevel) {
    G4cout << "Lambda table is built for "
           << particle->GetParticleName()
           << G4endl;
    if(2 < verboseLevel) {
      G4cout << *theLambdaTable << G4endl;
    }
  }
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
  if(p) p->SetParticleChange(pParticleChange);
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
  fParticleChange.InitializeForPostStep(track);

  // Do not make anything if particle is stopped, the annihilation then
  // should be performed by the AtRestDoIt!
  if (track.GetTrackStatus() == fStopButAlive) return &fParticleChange;

  G4double finalT = track.GetKineticEnergy();

  // Integral approach
  if (integral) {
    G4double lx = GetLambda(finalT, currentCouple);
    if(preStepLambda<lx && 1 < verboseLevel) {
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

  /*  
  if(0 < verboseLevel) {
    G4cout << "G4VEmProcess::PostStepDoIt: Sample secondary; E= "
           << finalT/MeV
           << " MeV; model= (" << currentModel->LowEnergyLimit()
           << ", " <<  currentModel->HighEnergyLimit() << ")"
           << G4endl;
  }
  */

  
  std::vector<G4DynamicParticle*>* newp = 
    SecondariesPostStep(currentModel,currentCouple,track.GetDynamicParticle());

  if (newp) {
    G4int num = newp->size();
    fParticleChange.SetNumberOfSecondaries(num);
    G4double edep = fParticleChange.GetLocalEnergyDeposit();
     
    for (G4int i=0; i<num; i++) {
      G4DynamicParticle* dp = (*newp)[i];
      const G4ParticleDefinition* p = dp->GetDefinition();
      G4double e = dp->GetKineticEnergy();
      G4bool good = true;
      if (p == theGamma) {
	if (e < (*theCutsGamma)[currentMaterialIndex]) good = false;

      } else if (p == theElectron) {
	if (e < (*theCutsElectron)[currentMaterialIndex]) good = false;

      } else if (p == thePositron) {
	if (e < (*theCutsPositron)[currentMaterialIndex]) {
	  good = false;
	  e += 2.0*electron_mass_c2;
	}
      }

      if (good || !applyCuts) {
	fParticleChange.AddSecondary(dp);
	
      } else {
	delete dp;
	edep += e;
      }
    } 
    fParticleChange.ProposeLocalEnergyDeposit(edep);
    delete newp;
  }

  return G4VDiscreteProcess::PostStepDoIt(track,step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::PrintInfoDefinition()
{
  if(verboseLevel > 0) {
    G4cout << G4endl << GetProcessName() << ": " ;
    PrintInfo();
  }

  if (!buildLambdaTable)  return;
  
  if(verboseLevel > 0) {
    G4cout << "      tables are built for  "
           << particle->GetParticleName()
           << G4endl
           << "      Lambda tables from "
           << G4BestUnit(minKinEnergy,"Energy") 
           << " to "
           << G4BestUnit(maxKinEnergy,"Energy")
           << " in " << nLambdaBins << " bins."
           << G4endl;
  }

  if(verboseLevel > 1) {
    G4cout << "Tables are built for " << particle->GetParticleName()
           << G4endl;

  if(verboseLevel > 2) {
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
  } else {
    G4VEmModel* model = SelectModel(kineticEnergy);
    cross = model->CrossSectionPerVolume(currentMaterial,particle,kineticEnergy);
  }

  return cross;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VEmProcess::ComputeCrossSectionPerAtom(G4double kineticEnergy, G4double Z)
{
  G4VEmModel* model = SelectModel(kineticEnergy);
  G4double A = G4NistManager::Instance()->FindOrBuildElement(G4int(Z), false)->GetN();
  G4double cross = model->ComputeCrossSectionPerAtom(particle,kineticEnergy,Z,A);
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

G4bool G4VEmProcess::StorePhysicsTable(const G4ParticleDefinition* part,
			 	       const G4String& directory,
				             G4bool ascii)
{
  G4bool yes = true;

  if ( theLambdaTable && part == particle) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    yes = theLambdaTable->StorePhysicsTable(name,ascii);

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
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VEmProcess::RetrievePhysicsTable(const G4ParticleDefinition* part,
			  	          const G4String& directory,
			  	                G4bool ascii)
{
  if(1 < verboseLevel) {
    G4cout << "G4VEmProcess::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
	   << GetProcessName() << G4endl;
  }
  G4bool yes = true;

  if(!buildLambdaTable || particle != part) return yes;

  const G4String particleName = part->GetParticleName();
  G4String filename;

  filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  yes = G4PhysicsTableHelper::RetrievePhysicsTable(theLambdaTable,filename,ascii);
  if ( yes ) {
    if (0 < verboseLevel) {
      G4cout << "Lambda table for " << particleName << " is Retrieved from <"
             << filename << ">"
             << G4endl;
    }
  } else {
    if (1 < verboseLevel) {
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
    G4cout << "### G4VEmProcess::FindLambdaMax: " << particle->GetParticleName() 
           << " and process " << GetProcessName() << G4endl; 
  }
  size_t n = theLambdaTable->length();
  G4PhysicsVector* pv = (*theLambdaTable)[0];
  G4double e, s, emax, smax;
  theEnergyOfCrossSectionMax = new G4double [n];
  theCrossSectionMax = new G4double [n];
  G4bool b;

  for (size_t i=0; i<n; i++) {
    pv = (*theLambdaTable)[i];
    emax = DBL_MAX;
    smax = 0.0;
    if(pv) {
      size_t nb = pv->GetVectorLength();
      emax = pv->GetLowEdgeEnergy(nb);
      smax = 0.0;
      for (size_t j=0; j<nb; j++) {
	e = pv->GetLowEdgeEnergy(j);
	s = pv->GetValue(e,b);
	if(s > smax) {
	  smax = s;
	  emax = e;
	}
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

void G4VEmProcess::ActivateDeexcitation(G4bool, const G4Region*)
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
  if(integral) buildLambdaTable = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4VEmProcess::IsIntegral() const
{
  return integral;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VEmProcess::LambdaPhysicsVector(const G4MaterialCutsCouple*)
{
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nLambdaBins);
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetBuildTableFlag(G4bool val)
{
  buildLambdaTable = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetStartFromNullFlag(G4bool val)
{
  startFromNull = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VEmProcess::SetApplyCuts(G4bool val)
{
  applyCuts = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
