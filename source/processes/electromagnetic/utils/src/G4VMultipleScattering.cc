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
// $Id: G4VMultipleScattering.cc 66594 2012-12-23 10:17:42Z vnivanch $
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
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"
#include "G4MaterialCutsCouple.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VEmFluctuationModel.hh"
#include "G4UnitsTable.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"
#include "G4TransportationManager.hh"
#include "G4SafetyHelper.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::G4VMultipleScattering(const G4String& name, G4ProcessType):
  G4VContinuousDiscreteProcess("msc", fElectromagnetic),
  numberOfModels(0),
  firstParticle(0),
  currParticle(0),
  stepLimit(fUseSafety),
  skin(1.0),
  facrange(0.04),
  facgeom(2.5),
  latDisplasment(true),
  isIon(false)
{
  SetVerboseLevel(1);
  SetProcessSubType(fMultipleScattering);
  if("ionmsc" == name) { firstParticle = G4GenericIon::GenericIon(); }

  geomMin = 1.e-6*CLHEP::mm;
  lowestKinEnergy = 10*eV;

  // default limit on polar angle
  polarAngleLimit = 0.0;

  physStepLimit = gPathLength = tPathLength = 0.0;
  fIonisation = 0;

  pParticleChange = &fParticleChange;
  safetyHelper = 0;
  fPositionChanged = false;
  isActive = false;

  modelManager = new G4EmModelManager();
  emManager = G4LossTableManager::Instance();
  emManager->Register(this);

  warn = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMultipleScattering::~G4VMultipleScattering()
{
  if(1 < verboseLevel) {
    G4cout << "G4VMultipleScattering destruct " << GetProcessName() 
	   << G4endl;
  }
  delete modelManager;
  emManager->DeRegister(this);
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
  ++warn;
  if(warn < 10) { 
    G4cout << "### G4VMultipleScattering::SetModel is obsolete method "
	   << "and will be removed for the next release." << G4endl;
    G4cout << "    Please, use SetEmModel" << G4endl;
  } 
  G4int n = mscModels.size();
  if(index >= n) { for(G4int i=n; i<=index; ++i) {mscModels.push_back(0);} }
  mscModels[index] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMscModel* G4VMultipleScattering::Model(G4int index)
{
  ++warn;
  if(warn < 10) { 
    G4cout << "### G4VMultipleScattering::Model is obsolete method "
	   << "and will be removed for the next release." << G4endl;
    G4cout << "    Please, use EmModel" << G4endl;
  } 
  G4VMscModel* p = 0;
  if(index >= 0 && index <  G4int(mscModels.size())) { p = mscModels[index]; }
  return p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::SetEmModel(G4VMscModel* p, G4int index)
{
  G4int n = mscModels.size();
  if(index >= n) { for(G4int i=n; i<=index; ++i) { mscModels.push_back(0); } }
  mscModels[index] = p;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VMscModel* G4VMultipleScattering::EmModel(G4int index)
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

void 
G4VMultipleScattering::PreparePhysicsTable(const G4ParticleDefinition& part)
{
  if(!firstParticle) { firstParticle = &part; }
  if(part.GetParticleType() == "nucleus") {
    SetStepLimitType(fMinimal);
    SetLateralDisplasmentFlag(false);
    SetRangeFactor(0.2);
    if(&part == G4GenericIon::GenericIon()) { firstParticle = &part; }
    isIon = true; 
  }

  emManager->PreparePhysicsTable(&part, this);
  currParticle = 0;

  if(1 < verboseLevel) {
    G4cout << "### G4VMultipleScattering::PrepearPhysicsTable() for "
           << GetProcessName()
           << " and particle " << part.GetParticleName()
           << " local particle " << firstParticle->GetParticleName()
	   << " isIon= " << isIon 
           << G4endl;
  }

  if(firstParticle == &part) {

    InitialiseProcess(firstParticle);

    // initialisation of models
    numberOfModels = modelManager->NumberOfModels();
    for(G4int i=0; i<numberOfModels; ++i) {
      G4VMscModel* msc = static_cast<G4VMscModel*>(modelManager->GetModel(i));
      msc->SetIonisation(0, firstParticle);
      if(0 == i) { currentModel = msc; }
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
      G4double emax = 
	std::min(msc->HighEnergyLimit(),emManager->MaxKinEnergy());
      msc->SetHighEnergyLimit(emax);
    }

    modelManager->Initialise(firstParticle, G4Electron::Electron(), 
			     10.0, verboseLevel);

    if(!safetyHelper) {
      safetyHelper = G4TransportationManager::GetTransportationManager()
	->GetSafetyHelper();
      safetyHelper->InitialiseHelper();
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

  emManager->BuildPhysicsTable(firstParticle);

  // explicitly defined printout by particle name
  if(1 < verboseLevel || 
     (0 < verboseLevel && (num == "e-" || 
			   num == "e+"    || num == "mu+" || 
			   num == "mu-"   || num == "proton"|| 
			   num == "pi+"   || num == "pi-" || 
			   num == "kaon+" || num == "kaon-" || 
			   num == "alpha" || num == "anti_proton" || 
			   num == "GenericIon")))
    { 
      G4cout << G4endl << GetProcessName() 
	     << ":   for " << num
	     << "    SubType= " << GetProcessSubType() 
	     << G4endl;
      PrintInfo();
      modelManager->DumpModelList(verboseLevel);
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
    PrintInfo();
    modelManager->DumpModelList(verboseLevel);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::StartTracking(G4Track* track)
{
  G4VEnergyLossProcess* eloss = 0;
  if(track->GetParticleDefinition() != currParticle) {
    currParticle = track->GetParticleDefinition();
    fIonisation = emManager->GetEnergyLossProcess(currParticle);
    eloss = fIonisation;
  }
  // one model
  if(1 == numberOfModels) {
    currentModel->StartTracking(track);
    if(eloss) { currentModel->SetIonisation(fIonisation, currParticle); }

    // many models
  } else { 
    for(G4int i=0; i<numberOfModels; ++i) {
      G4VMscModel* msc = static_cast<G4VMscModel*>(modelManager->GetModel(i));
      msc->StartTracking(track);
      if(eloss) { msc->SetIonisation(fIonisation, currParticle); }
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
  physStepLimit = gPathLength = tPathLength = currentMinimalStep;

  G4double ekin = track.GetKineticEnergy();
  // isIon flag is used only to select a model
  if(isIon) { 
    ekin *= proton_mass_c2/track.GetParticleDefinition()->GetPDGMass(); 
  }
  
  // select new model
  if(1 < numberOfModels) {
    currentModel = static_cast<G4VMscModel*>(
      SelectModel(ekin,track.GetMaterialCutsCouple()->GetIndex()));
  }

  // step limit
  if(currentModel->IsActive(ekin) && gPathLength >= geomMin 
     && ekin >= lowestKinEnergy) {
    isActive = true;
    tPathLength = currentModel->ComputeTruePathLengthLimit(track, gPathLength);
    if (tPathLength < physStepLimit) { 
      *selection = CandidateForSelection; 
    }
  } else { isActive = false; }
  /*  
  if(currParticle->GetPDGMass() > GeV)    
  G4cout << "MSC::AlongStepGPIL: Ekin= " << ekin
	 << " gPathLength= " << gPathLength
	 << " tPathLength= " << tPathLength
	 << " currentMinimalStep= " << currentMinimalStep
	 << " isActive " << isActive << G4endl;
  */
  return gPathLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double 
G4VMultipleScattering::PostStepGetPhysicalInteractionLength(
              const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = Forced;
  //*condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4VMultipleScattering::AlongStepDoIt(const G4Track& track, const G4Step& step)
{
  fParticleChange.ProposeMomentumDirection(step.GetPostStepPoint()->GetMomentumDirection());
  fNewPosition = step.GetPostStepPoint()->GetPosition();
  fParticleChange.ProposePosition(fNewPosition);
  fPositionChanged = false;

  G4double geomLength = step.GetStepLength();

  // very small step - no msc
  if(!isActive) {
    tPathLength = geomLength;

    // sample msc
  } else {
    G4double range = 
      currentModel->GetRange(currParticle,track.GetKineticEnergy(),
			     track.GetMaterialCutsCouple());

    G4double trueLength = currentModel->ComputeTrueStepLength(geomLength);

    
    // protection against wrong t->g->t conversion
    //    if(trueLength > tPathLength) 
    /*
    if(currParticle->GetPDGMass() > GeV)    
    G4cout << "G4VMsc::AlongStepDoIt: GeomLength= " 
	   << geomLength 
	   << " trueLenght= " << trueLength 
	   << " tPathLength= " << tPathLength
	   << " dr= " << range - trueLength 
	   << " ekin= " << track.GetKineticEnergy() << G4endl;
    */
    if (trueLength <= physStepLimit) {
      tPathLength = trueLength; 
    } else {
      tPathLength = physStepLimit - 0.5*geomMin; 
    }

    // do not sample scattering at the last or at a small step
    if(tPathLength + geomMin < range && tPathLength > geomMin) {

      G4double preSafety = step.GetPreStepPoint()->GetSafety();
      G4double postSafety= preSafety - geomLength; 
      G4bool safetyRecomputed = false;
      if( postSafety < geomMin ) {
	safetyRecomputed = true;
	postSafety = currentModel->ComputeSafety(fNewPosition,0.0); 
      } 
      G4ThreeVector displacement = 
	currentModel->SampleScattering(step.GetPostStepPoint()->GetMomentumDirection(),postSafety);

      G4double r2 = displacement.mag2();

      //G4cout << "R= " << sqrt(r2) << " postSafety= " << postSafety << G4endl;

      // make correction for displacement
      if(r2 > 0.0) {

	fPositionChanged = true;
        G4double fac = 1.0;

	// displaced point is definitely within the volume
	if(r2 > postSafety*postSafety) {

          G4double dispR = std::sqrt(r2);
	  if(!safetyRecomputed) {
	    postSafety = currentModel->ComputeSafety(fNewPosition, dispR);
	  } 
	  // add a factor which ensure numerical stability
	  if(dispR > postSafety) { fac = 0.99*postSafety/dispR; }
	}
	// compute new endpoint of the Step
	fNewPosition += fac*displacement;
	//safetyHelper->ReLocateWithinVolume(fNewPosition);
      }
    }
  }
  fParticleChange.ProposeTrueStepLength(tPathLength);
  //fParticleChange.ProposePosition(fNewPosition);
  return &fParticleChange;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* 
G4VMultipleScattering::PostStepDoIt(const G4Track& track, const G4Step&)
{
  fParticleChange.Initialize(track);  
  
  if(fPositionChanged) { 
    safetyHelper->ReLocateWithinVolume(fNewPosition);
    fParticleChange.ProposePosition(fNewPosition); 
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4VMultipleScattering::StorePhysicsTable(const G4ParticleDefinition* part,
					 const G4String& directory,
					 G4bool ascii)
{
  G4bool yes = true;
  if(part != firstParticle) { return yes; }
  G4int nmod = modelManager->NumberOfModels();
  const G4String ss[4] = {"1","2","3","4"};
  for(G4int i=0; i<nmod; ++i) {
    G4VEmModel* msc = modelManager->GetModel(i);
    yes = true;
    G4PhysicsTable* table = msc->GetCrossSectionTable();
    if (table) {
      G4int j = std::min(i,3); 
      G4String name = 
	GetPhysicsTableFileName(part,directory,"LambdaMod"+ss[j],ascii);
      yes = table->StorePhysicsTable(name,ascii);

      if ( yes ) {
	if ( verboseLevel>0 ) {
	  G4cout << "Physics table are stored for " << part->GetParticleName()
		 << " and process " << GetProcessName()
		 << " with a name <" << name << "> " << G4endl;
	}
      } else {
	G4cout << "Fail to store Physics Table for " << part->GetParticleName()
	       << " and process " << GetProcessName()
	       << " in the directory <" << directory
	       << "> " << G4endl;
      }
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool 
G4VMultipleScattering::RetrievePhysicsTable(const G4ParticleDefinition*,
					    const G4String&,
					    G4bool)
{
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4VMultipleScattering::SetIonisation(G4VEnergyLossProcess* p)
{
  for(G4int i=0; i<numberOfModels; ++i) {
    G4VMscModel* msc = static_cast<G4VMscModel*>(modelManager->GetModel(i));
    msc->SetIonisation(p, firstParticle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

