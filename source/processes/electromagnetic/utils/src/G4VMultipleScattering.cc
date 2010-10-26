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
// $Id: G4VMultipleScattering.cc,v 1.86 2010-10-26 11:30:46 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
// 01-03-04 SampleCosineTheta signature changed
// 22-04-04 SampleCosineTheta signature changed back to original
// 27-08-04 Add InitialiseForRun method (V.Ivanchneko)
// 08-11-04 Migration to new interface of Store/Retrieve tables (V.Ivantchenko)
// 11-03-05 Shift verbose level by 1 (V.Ivantchenko)
// 15-04-05 optimize internal interface (V.Ivanchenko)
// 15-04-05 remove boundary flag (V.Ivanchenko)
// 27-10-05 introduce virtual function MscStepLimitation() (V.Ivanchenko)
// 12-04-07 Add verbosity at destruction (V.Ivanchenko)
// 27-10-07 Virtual functions moved to source (V.Ivanchenko)
// 11-03-08 Set skin value does not effect step limit type (V.Ivanchenko)
// 24-06-09 Removed hidden bin in G4PhysicsVector (V.Ivanchenko)
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
#include "G4PhysicsTableHelper.hh"
#include "G4GenericIon.hh"
#include "G4Electron.hh"
#include "G4EmConfigurator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::G4VMultipleScattering(const G4String& name, 
					     G4ProcessType type):
  G4VContinuousDiscreteProcess(name, type),
  buildLambdaTable(true),
  theLambdaTable(0),
  firstParticle(0),
  stepLimit(fUseSafety),
  skin(1.0),
  facrange(0.04),
  facgeom(2.5),
  latDisplasment(true),
  isIon(false),
  currentParticle(0),
  currentCouple(0)
{
  SetVerboseLevel(1);
  SetProcessSubType(fMultipleScattering);

  // Size of tables assuming spline
  minKinEnergy = 0.1*keV;
  maxKinEnergy = 10.0*TeV;
  nBins        = 77;

  // default limit on polar angle
  polarAngleLimit = 0.0;

  pParticleChange = &fParticleChange;

  modelManager = new G4EmModelManager();
  (G4LossTableManager::Instance())->Register(this);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::~G4VMultipleScattering()
{
  if(1 < verboseLevel) {
    G4cout << "G4VMultipleScattering destruct " << GetProcessName() 
	   << G4endl;
  }
  delete modelManager;
  if (theLambdaTable) {
    theLambdaTable->clearAndDestroy();
    delete theLambdaTable;
  }
  (G4LossTableManager::Instance())->DeRegister(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::AddEmModel(G4int order, G4VEmModel* p,
				       const G4Region* region)
{
  G4VEmFluctuationModel* fm = 0;
  modelManager->AddEmModel(order, p, fm, region);
  if(p) { p->SetParticleChange(pParticleChange); }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::SetModel(G4VMscModel* p, G4int index)
{
  G4int n = mscModels.size();
  if(index >= n) { for(G4int i=n; i<=index; ++i) {mscModels.push_back(0);} }
  mscModels[index] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMscModel* G4VMultipleScattering::Model(G4int index)
{
  G4VMscModel* p = 0;
  if(index >= 0 && index <  G4int(mscModels.size())) { p = mscModels[index]; }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VEmModel* 
G4VMultipleScattering::GetModelByIndex(G4int idx, G4bool ver) const
{
  return modelManager->GetModel(idx, ver);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if (!firstParticle) {
    currentCouple = 0;
    if(part.GetParticleType() == "nucleus" && 
       part.GetParticleSubType() == "generic") {
      firstParticle = G4GenericIon::GenericIon();
      isIon = true;
    } else {
      firstParticle = &part;
      if(part.GetParticleType() == "nucleus" || 
	 part.GetPDGMass() > GeV) {isIon = true;} 
    } 
    // limitations for ions
    if(isIon) {
      SetStepLimitType(fMinimal);
      SetLateralDisplasmentFlag(false);
      SetBuildLambdaTable(false);
    }
    currentParticle = &part;
  }

  (G4LossTableManager::Instance())->PreparePhysicsTable(&part, this);

  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << " local particle " << firstParticle->GetParticleName()
           << G4endl;
  }

  if(firstParticle == &part) {

    InitialiseProcess(firstParticle);

    // initialisation of models
    G4int nmod = modelManager->NumberOfModels();
    for(G4int i=0; i<nmod; ++i) {
      G4VMscModel* msc = static_cast<G4VMscModel*>(modelManager->GetModel(i));
      if(isIon) {
	msc->SetStepLimitType(fMinimal);
	msc->SetLateralDisplasmentFlag(false);
	msc->SetRangeFactor(0.2);
      } else {
	msc->SetStepLimitType(StepLimitType());
	msc->SetLateralDisplasmentFlag(LateralDisplasmentFlag());
	msc->SetSkin(Skin());
	msc->SetRangeFactor(RangeFactor());
	msc->SetGeomFactor(GeomFactor());
      }
      msc->SetPolarAngleLimit(polarAngleLimit);
      if(msc->HighEnergyLimit() > maxKinEnergy) {
	msc->SetHighEnergyLimit(maxKinEnergy);
      }
    }

    modelManager->Initialise(firstParticle, G4Electron::Electron(), 
			     10.0, verboseLevel);

    // prepare tables
    if(buildLambdaTable) {
      theLambdaTable = G4PhysicsTableHelper::PreparePhysicsTable(theLambdaTable);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::BuildPhysicsTable(const G4ParticleDefinition& part)
{
  G4String num = part.GetParticleName();
  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() for "
           << GetProcessName()
           << " and particle " << num
           << G4endl;
  }

  (G4LossTableManager::Instance())->BuildPhysicsTable(firstParticle);

  if (buildLambdaTable && firstParticle == &part) {

    const G4ProductionCutsTable* theCoupleTable=
          G4ProductionCutsTable::GetProductionCutsTable();
    size_t numOfCouples = theCoupleTable->GetTableSize();

    G4bool splineFlag = (G4LossTableManager::Instance())->SplineFlag();

    G4PhysicsLogVector* aVector = 0;
    G4PhysicsLogVector* bVector = 0;

    for (size_t i=0; i<numOfCouples; ++i) {

      if (theLambdaTable->GetFlag(i)) {
        // create physics vector and fill it
        const G4MaterialCutsCouple* couple = 
	  theCoupleTable->GetMaterialCutsCouple(i);
	if(!bVector) {
	  aVector = static_cast<G4PhysicsLogVector*>(PhysicsVector(couple));
	  bVector = aVector;
	} else {
	  aVector = new G4PhysicsLogVector(*bVector);
	}       
	//G4PhysicsVector* aVector = PhysicsVector(couple);
	aVector->SetSpline(splineFlag);
        modelManager->FillLambdaVector(aVector, couple, false);
	if(splineFlag) aVector->FillSecondDerivatives();
        G4PhysicsTableHelper::SetPhysicsVector(theLambdaTable, i, aVector);
      }
    }

    if(1 < verboseLevel) {
      G4cout << "Lambda table is built for "
             << num
             << G4endl;
    }
  }
  if(verboseLevel>0 && ( num == "e-" || num == "mu+" ||  
                         num == "proton" || num == "pi-" || 
			 num == "kaon+" || num == "GenericIon")) {
    PrintInfoDefinition();
    if(2 < verboseLevel && theLambdaTable) G4cout << *theLambdaTable << G4endl;
  }

  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::BuildPhysicsTable() done for "
           << GetProcessName()
           << " and particle " << num
           << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::PrintInfoDefinition()
{
  if (0 < verboseLevel) {
    G4cout << G4endl << GetProcessName() 
	   << ":   for " << firstParticle->GetParticleName()
	   << "    SubType= " << GetProcessSubType() 
	   << G4endl;
    if (theLambdaTable) {
      G4cout << "      Lambda tables from "
	     << G4BestUnit(MinKinEnergy(),"Energy")
	     << " to "
	     << G4BestUnit(MaxKinEnergy(),"Energy")
	     << " in " << nBins << " bins, spline: " 
	     << (G4LossTableManager::Instance())->SplineFlag()
	     << G4endl;
    }
    PrintInfo();
    modelManager->DumpModelList(verboseLevel);
    if (2 < verboseLevel) {
      G4cout << "LambdaTable address= " << theLambdaTable << G4endl;
      if(theLambdaTable) G4cout << (*theLambdaTable) << G4endl;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4VMultipleScattering::AlongStepGetPhysicalInteractionLength(
                             const G4Track& track,
                             G4double,
                             G4double currentMinimalStep,
                             G4double&,
                             G4GPILSelection* selection)
{
  // get Step limit proposed by the process
  *selection = NotCandidateForSelection;
  G4double x = currentMinimalStep;
  DefineMaterial(track.GetMaterialCutsCouple());
  G4double ekin = track.GetKineticEnergy();
  if(isIon) { ekin *= proton_mass_c2/track.GetParticleDefinition()->GetPDGMass(); }
  currentModel = static_cast<G4VMscModel*>(SelectModel(ekin));
  if(x > 0.0 && ekin > 0.0 && currentModel->IsActive(ekin)) {
    G4double tPathLength = 
      currentModel->ComputeTruePathLengthLimit(track, theLambdaTable, x);
    if (tPathLength < x) { *selection = CandidateForSelection; }
    x = currentModel->ComputeGeomPathLength(tPathLength);
    //  G4cout << "tPathLength= " << tPathLength
    //         << " stepLimit= " << x
    //        << " currentMinimalStep= " << currentMinimalStep<< G4endl;
  }
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VMultipleScattering::PostStepGetPhysicalInteractionLength(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4VMultipleScattering::AlongStepDoIt(const G4Track& track, const G4Step& step)
{
  if(currentModel->IsActive(track.GetKineticEnergy())) {
    fParticleChange.ProposeTrueStepLength(currentModel->ComputeTrueStepLength(step.GetStepLength()));
  }
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4VMultipleScattering::PostStepDoIt(const G4Track& track, const G4Step& step)
{
  fParticleChange.Initialize(track);
  if(currentModel->IsActive(track.GetKineticEnergy())) {
    currentModel->SampleScattering(track.GetDynamicParticle(),
				   step.GetPostStepPoint()->GetSafety());
  }
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::GetContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  G4GPILSelection selection = NotCandidateForSelection;
  G4double x = AlongStepGetPhysicalInteractionLength(track,previousStepSize,
						     currentMinimalStep,
						     currentSafety, &selection);
  return x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::ContinuousStepLimit(
                                       const G4Track& track,
                                       G4double previousStepSize,
                                       G4double currentMinimalStep,
                                       G4double& currentSafety)
{
  return GetContinuousStepLimit(track,previousStepSize,currentMinimalStep,
				currentSafety);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4VMultipleScattering::GetMeanFreePath(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4PhysicsVector* G4VMultipleScattering::PhysicsVector(const G4MaterialCutsCouple* couple)
{
  G4int nbins = 3;
  if( couple->IsUsed() ) nbins = nBins;
  G4PhysicsVector* v = new G4PhysicsLogVector(minKinEnergy, maxKinEnergy, nbins);
  v->SetSpline((G4LossTableManager::Instance())->SplineFlag());
  return v;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4VMultipleScattering::StorePhysicsTable(const G4ParticleDefinition* part,
			 	                const G4String& directory,
				                      G4bool ascii)
{
  G4bool yes = true;
  if ( theLambdaTable && part == firstParticle) {
    const G4String name = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
    G4bool yes = theLambdaTable->StorePhysicsTable(name,ascii);

    if ( yes ) {
      if ( verboseLevel>0 ) {
        G4cout << "Physics table are stored for " << part->GetParticleName()
               << " and process " << GetProcessName()
	       << " in the directory <" << directory
	       << "> " << G4endl;
      }
    } else {
      G4cout << "Fail to store Physics Table for " << part->GetParticleName()
             << " and process " << GetProcessName()
	     << " in the directory <" << directory
	     << "> " << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4VMultipleScattering::RetrievePhysicsTable(const G4ParticleDefinition* part,
					    const G4String& directory,
					    G4bool ascii)
{
  if(0 < verboseLevel) {
    G4cout << "G4VMultipleScattering::RetrievePhysicsTable() for "
           << part->GetParticleName() << " and process "
	   << GetProcessName() << G4endl;
  }
  G4bool yes = true;

  if(!buildLambdaTable || firstParticle != part) return yes;

  const G4String particleName = part->GetParticleName();

  G4String filename = GetPhysicsTableFileName(part,directory,"Lambda",ascii);
  yes = 
    G4PhysicsTableHelper::RetrievePhysicsTable(theLambdaTable,filename,ascii);
  if ( yes ) {
    if (0 < verboseLevel) {
        G4cout << "Lambda table for " << part->GetParticleName() 
	       << " is retrieved from <"
               << filename << ">"
               << G4endl;
    }
    if((G4LossTableManager::Instance())->SplineFlag()) {
      size_t n = theLambdaTable->length();
      for(size_t i=0; i<n; ++i) {
        if((* theLambdaTable)[i]) {
	  (* theLambdaTable)[i]->SetSpline(true);
	}
      }
    }
  } else {
    if (1 < verboseLevel) {
        G4cout << "Lambda table for " << part->GetParticleName() 
	       << " in file <"
               << filename << "> is not exist"
               << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

