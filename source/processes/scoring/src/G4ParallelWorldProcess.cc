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
// $Id: G4ParallelWorldProcess.cc 103731 2017-04-25 08:01:05Z gcosmo $
// GEANT4 tag $Name: geant4-09-04-ref-00 $
//
//

#include "G4ios.hh"
#include "G4ParallelWorldProcess.hh"
#include "G4ParallelWorldProcessStore.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4Navigator.hh"
#include "G4VTouchable.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleChange.hh"
#include "G4PathFinder.hh"
#include "G4TransportationManager.hh"
#include "G4ParticleChange.hh"
#include "G4StepPoint.hh"
#include "G4FieldTrackUpdator.hh"
#include "G4Material.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"

#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

G4ThreadLocal G4Step* G4ParallelWorldProcess::fpHyperStep = 0;
G4ThreadLocal G4int   G4ParallelWorldProcess::nParallelWorlds = 0;
G4ThreadLocal G4int   G4ParallelWorldProcess::fNavIDHyp = 0;
const G4Step* G4ParallelWorldProcess::GetHyperStep()
{ return fpHyperStep; }
G4int G4ParallelWorldProcess::GetHypNavigatorID()
{ return fNavIDHyp; }

G4ParallelWorldProcess::
G4ParallelWorldProcess(const G4String& processName,G4ProcessType theType)
:G4VProcess(processName,theType),fGhostWorld(nullptr),fGhostNavigator(nullptr),
 fNavigatorID(-1),fFieldTrack('0'),fGhostSafety(0.),fOnBoundary(false),
 layeredMaterialFlag(false)
{
  SetProcessSubType(491);
  if(!fpHyperStep) fpHyperStep = new G4Step();
  iParallelWorld = ++nParallelWorlds;

  pParticleChange = &aDummyParticleChange;

  fGhostStep = new G4Step();
  fGhostPreStepPoint = fGhostStep->GetPreStepPoint();
  fGhostPostStepPoint = fGhostStep->GetPostStepPoint();

  fTransportationManager = G4TransportationManager::GetTransportationManager();
  fTransportationManager->GetNavigatorForTracking()->SetPushVerbosity(false);
  fPathFinder = G4PathFinder::GetInstance();

  fGhostWorldName = "** NotDefined **";
  G4ParallelWorldProcessStore::GetInstance()->SetParallelWorld(this,processName);

  if (verboseLevel>0)
  {
    G4cout << GetProcessName() << " is created " << G4endl;
  }
}

G4ParallelWorldProcess::~G4ParallelWorldProcess()
{
  delete fGhostStep;
  nParallelWorlds--;
  if(nParallelWorlds==0)
  {
    delete fpHyperStep;
    fpHyperStep = 0;
  }
}

void G4ParallelWorldProcess::
SetParallelWorld(G4String parallelWorldName)
{
  fGhostWorldName = parallelWorldName;
  fGhostWorld = fTransportationManager->GetParallelWorld(fGhostWorldName);
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
  fGhostNavigator->SetPushVerbosity(false);
}

void G4ParallelWorldProcess::
SetParallelWorld(G4VPhysicalVolume* parallelWorld)
{
  fGhostWorldName = parallelWorld->GetName();
  fGhostWorld = parallelWorld;
  fGhostNavigator = fTransportationManager->GetNavigator(fGhostWorld);
  fGhostNavigator->SetPushVerbosity(false);
}

void G4ParallelWorldProcess::StartTracking(G4Track* trk)
{
  if(fGhostNavigator)
  { fNavigatorID = fTransportationManager->ActivateNavigator(fGhostNavigator); }
  else
  {
    G4Exception("G4ParallelWorldProcess::StartTracking",
       "ProcParaWorld000",FatalException,
       "G4ParallelWorldProcess is used for tracking without having a parallel world assigned");
  }
  fPathFinder->PrepareNewTrack(trk->GetPosition(),trk->GetMomentumDirection());

  fOldGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
  fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
  fNewGhostTouchable = fOldGhostTouchable;
  fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);

  fGhostSafety = -1.;
  fOnBoundary = false;
  fGhostPreStepPoint->SetStepStatus(fUndefined);
  fGhostPostStepPoint->SetStepStatus(fUndefined);

//  G4VPhysicalVolume* thePhys = fNewGhostTouchable->GetVolume();
//  if(thePhys)
//  {
//    G4Material* ghostMaterial = thePhys->GetLogicalVolume()->GetMaterial();
//    if(ghostMaterial)
//    { G4cout << " --- Material : " << ghostMaterial->GetName() << G4endl; }
//  }

  *(fpHyperStep->GetPostStepPoint()) = *(trk->GetStep()->GetPostStepPoint());
  if(layeredMaterialFlag)
  {
    G4StepPoint* realWorldPostStepPoint = trk->GetStep()->GetPostStepPoint();
    SwitchMaterial(realWorldPostStepPoint);
  }
  *(fpHyperStep->GetPreStepPoint()) = *(fpHyperStep->GetPostStepPoint());
}

G4double 
G4ParallelWorldProcess::AtRestGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4ForceCondition* condition)
{
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// At Rest must be registered ONLY for the particle which has other At Rest
// process(es).
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  *condition = Forced;
  return DBL_MAX;
}

G4VParticleChange* G4ParallelWorldProcess::AtRestDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// At Rest must be registered ONLY for the particle which has other At Rest
// process(es).
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
  G4VSensitiveDetector* aSD = 0;
  if(fOldGhostTouchable->GetVolume())
  { aSD = fOldGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector(); }
  fOnBoundary = false;
  if(aSD)
  {
    CopyStep(step);
    fGhostPreStepPoint->SetSensitiveDetector(aSD);

    fNewGhostTouchable = fOldGhostTouchable;
  
    fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
    fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);
    if(fNewGhostTouchable->GetVolume())
    {
      fGhostPostStepPoint->SetSensitiveDetector(
        fNewGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector());
    }
    else
    { fGhostPostStepPoint->SetSensitiveDetector(0); }

    aSD->Hit(fGhostStep);
  }

  pParticleChange->Initialize(track);
  return pParticleChange;
}

G4double 
G4ParallelWorldProcess::PostStepGetPhysicalInteractionLength(
         const G4Track& /*track*/, 
         G4double   /*previousStepSize*/, 
         G4ForceCondition* condition)
{
  *condition = StronglyForced;
  return DBL_MAX;
}

G4VParticleChange* G4ParallelWorldProcess::PostStepDoIt(
     const G4Track& track,
     const G4Step& step)
{ 
  fOldGhostTouchable = fGhostPostStepPoint->GetTouchableHandle();
  G4VSensitiveDetector* aSD = 0;
  if(fOldGhostTouchable->GetVolume())
  { aSD = fOldGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector(); }
  CopyStep(step);
  fGhostPreStepPoint->SetSensitiveDetector(aSD);

  if(fOnBoundary)
  {
    fNewGhostTouchable = fPathFinder->CreateTouchableHandle(fNavigatorID);
  }
  else
  {
    fNewGhostTouchable = fOldGhostTouchable;
  }
    
  fGhostPreStepPoint->SetTouchableHandle(fOldGhostTouchable);
  fGhostPostStepPoint->SetTouchableHandle(fNewGhostTouchable);

  if(fNewGhostTouchable->GetVolume())
  {
    fGhostPostStepPoint->SetSensitiveDetector(
      fNewGhostTouchable->GetVolume()->GetLogicalVolume()->GetSensitiveDetector());
  }
  else
  { fGhostPostStepPoint->SetSensitiveDetector(0); }

  G4VSensitiveDetector* sd = fGhostPreStepPoint->GetSensitiveDetector();
  if(sd)
  {
    sd->Hit(fGhostStep);
  }

  pParticleChange->Initialize(track); 
  if(layeredMaterialFlag)
  {
    G4StepPoint* realWorldPostStepPoint =
     ((G4Step*)(track.GetStep()))->GetPostStepPoint();
    SwitchMaterial(realWorldPostStepPoint);
  }
  return pParticleChange;
}

G4double G4ParallelWorldProcess::AlongStepGetPhysicalInteractionLength(
            const G4Track& track, G4double  previousStepSize, G4double  currentMinimumStep,
            G4double& proposedSafety, G4GPILSelection* selection)
{
  static G4ThreadLocal G4FieldTrack *endTrack_G4MT_TLS_ = 0 ; if (!endTrack_G4MT_TLS_) endTrack_G4MT_TLS_ = new  G4FieldTrack ('0') ;  G4FieldTrack &endTrack = *endTrack_G4MT_TLS_;
  //static ELimited eLimited;
  ELimited eLimited;
  ELimited eLim = kUndefLimited;
  
  *selection = NotCandidateForSelection;
  G4double returnedStep = DBL_MAX;

  if (previousStepSize > 0.)
  { fGhostSafety -= previousStepSize; }
  if (fGhostSafety < 0.) fGhostSafety = 0.0;
      
  if (currentMinimumStep <= fGhostSafety && currentMinimumStep > 0.)
  {
    // I have no chance to limit
    returnedStep = currentMinimumStep;
    fOnBoundary = false;
    proposedSafety = fGhostSafety - currentMinimumStep;
    eLim = kDoNot;
  }
  else 
  {
    G4FieldTrackUpdator::Update(&fFieldTrack,&track);

#ifdef G4DEBUG_PARALLEL_WORLD_PROCESS    
    if( verboseLevel > 0 ){
       int localVerb = verboseLevel-1; 

       if( localVerb == 1 ) { 
          G4cout << " Pll Wrl  proc::AlongStepGPIL " << this->GetProcessName() << G4endl;
       }else if( localVerb > 1 ) { 
          G4cout << "----------------------------------------------" << G4endl;
          G4cout << " ParallelWorldProcess: field Track set to : " << G4endl;
          G4cout << "----------------------------------------------" << G4endl;
          G4cout << fFieldTrack << G4endl;
          G4cout << "----------------------------------------------" << G4endl;
       }
    }
#endif
    
    returnedStep
      = fPathFinder->ComputeStep(fFieldTrack,currentMinimumStep,fNavigatorID,
                     track.GetCurrentStepNumber(),fGhostSafety,eLimited,
                     endTrack,track.GetVolume());
    if(eLimited == kDoNot)
    {
      fOnBoundary = false;
      fGhostSafety = fGhostNavigator->ComputeSafety(endTrack.GetPosition());      
    }
    else
    {
      fOnBoundary = true;
      // fGhostSafetyEnd = 0.0;    // At end-point of expected step only
    }
    proposedSafety = fGhostSafety;
    if(eLimited == kUnique || eLimited == kSharedOther) {
       *selection = CandidateForSelection;
    }
    else if (eLimited == kSharedTransport) { 
       returnedStep *= (1.0 + 1.0e-9);  
    }
    eLim = eLimited;
  }

  if(iParallelWorld==nParallelWorlds) fNavIDHyp = 0;
  if(eLim == kUnique || eLim == kSharedOther) fNavIDHyp = fNavigatorID;
  return returnedStep;
}

G4VParticleChange* G4ParallelWorldProcess::AlongStepDoIt(
    const G4Track& track, const G4Step& )
{
  pParticleChange->Initialize(track);
  return pParticleChange; 
}

void G4ParallelWorldProcess::CopyStep(const G4Step & step)
{
  G4StepStatus prevStat = fGhostPostStepPoint->GetStepStatus();

  fGhostStep->SetTrack(step.GetTrack());
  fGhostStep->SetStepLength(step.GetStepLength());
  fGhostStep->SetTotalEnergyDeposit(step.GetTotalEnergyDeposit());
  fGhostStep->SetNonIonizingEnergyDeposit(step.GetNonIonizingEnergyDeposit());
  fGhostStep->SetControlFlag(step.GetControlFlag());
  fGhostStep->SetSecondary((const_cast<G4Step&>(step)).GetfSecondary());

  *fGhostPreStepPoint = *(step.GetPreStepPoint());
  *fGhostPostStepPoint = *(step.GetPostStepPoint());

  fGhostPreStepPoint->SetStepStatus(prevStat);
  if(fOnBoundary)
  { fGhostPostStepPoint->SetStepStatus(fGeomBoundary); }
  else if(fGhostPostStepPoint->GetStepStatus()==fGeomBoundary)
  { fGhostPostStepPoint->SetStepStatus(fPostStepDoItProc); }

  if(iParallelWorld==1)
  {
    G4StepStatus prevStatHyp = fpHyperStep->GetPostStepPoint()->GetStepStatus();

    fpHyperStep->SetTrack(step.GetTrack());
    fpHyperStep->SetStepLength(step.GetStepLength());
    fpHyperStep->SetTotalEnergyDeposit(step.GetTotalEnergyDeposit());
    fpHyperStep->SetNonIonizingEnergyDeposit(step.GetNonIonizingEnergyDeposit());
    fpHyperStep->SetControlFlag(step.GetControlFlag());

    *(fpHyperStep->GetPreStepPoint()) = *(fpHyperStep->GetPostStepPoint());
    *(fpHyperStep->GetPostStepPoint()) = *(step.GetPostStepPoint());
  
    fpHyperStep->GetPreStepPoint()->SetStepStatus(prevStatHyp);
  }

  if(fOnBoundary)
  { fpHyperStep->GetPostStepPoint()->SetStepStatus(fGeomBoundary); }
}

void G4ParallelWorldProcess::SwitchMaterial(G4StepPoint* realWorldStepPoint)
{
  if(realWorldStepPoint->GetStepStatus()==fWorldBoundary) return;
  G4VPhysicalVolume* thePhys = fNewGhostTouchable->GetVolume();
  if(thePhys)
  {
    G4Material* ghostMaterial = thePhys->GetLogicalVolume()->GetMaterial();
    if(ghostMaterial)
    {
      G4Region* ghostRegion = thePhys->GetLogicalVolume()->GetRegion();
      G4ProductionCuts* prodCuts =
          realWorldStepPoint->GetMaterialCutsCouple()->GetProductionCuts();
      if(ghostRegion)
      {
        G4ProductionCuts* ghostProdCuts = ghostRegion->GetProductionCuts();
        if(ghostProdCuts) prodCuts = ghostProdCuts;
      }
      const G4MaterialCutsCouple* ghostMCCouple =
          G4ProductionCutsTable::GetProductionCutsTable()
          ->GetMaterialCutsCouple(ghostMaterial,prodCuts);
      if(ghostMCCouple)
      {
        realWorldStepPoint->SetMaterial(ghostMaterial);
        realWorldStepPoint->SetMaterialCutsCouple(ghostMCCouple);
        *(fpHyperStep->GetPostStepPoint()) = *(fGhostPostStepPoint);
        fpHyperStep->GetPostStepPoint()->SetMaterial(ghostMaterial);
        fpHyperStep->GetPostStepPoint()->SetMaterialCutsCouple(ghostMCCouple);
      }
      else
      {
        G4cout << "!!! MaterialCutsCouple is not found for "
               << ghostMaterial->GetName() << "." << G4endl
               << "    Material in real world ("
               << realWorldStepPoint->GetMaterial()->GetName()
               << ") is used." << G4endl;
      }
    }
  }
}

G4bool G4ParallelWorldProcess::IsAtRestRequired(G4ParticleDefinition* partDef)
{
  G4int pdgCode = partDef->GetPDGEncoding();
  if(pdgCode==0)
  {
    G4String partName = partDef->GetParticleName();
    if(partName=="opticalphoton") return false;
    if(partName=="geantino") return false;
    if(partName=="chargedgeantino") return false;
  }
  else
  {
    if(pdgCode==22) return false; // gamma
    if(pdgCode==11) return false; // electron
    if(pdgCode==2212) return false; // proton
    if(pdgCode==-12) return false; // anti_nu_e
    if(pdgCode==12) return false; // nu_e
    if(pdgCode==-14) return false; // anti_nu_mu
    if(pdgCode==14) return false; // nu_mu
    if(pdgCode==-16) return false; // anti_nu_tau
    if(pdgCode==16) return false; // nu_tau
  }
  return true;
}

