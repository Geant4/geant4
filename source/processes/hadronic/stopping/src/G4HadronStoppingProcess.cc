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
//---------------------------------------------------------------------
//
// GEANT4 Class 
//
// File name:     G4HadronStoppingProcess
//
// Author V.Ivanchenko 21 April 2012 
//
//
// Class Description:
//
// Base process class for nuclear capture of negatively charged particles
//
// Modifications:
//
//  20120522  M. Kelsey -- Set enableAtRestDoIt flag for G4ProcessManager
//  20120914  M. Kelsey -- Pass subType in base ctor, remove enable flags
//  20121004  K. Genser -- use G4HadronicProcessType in the constructor
//  20121016  K. Genser -- Reverting to use one argument c'tor
//  20140818  K. Genser -- Labeled tracks using G4PhysicsModelCatalog
//
//------------------------------------------------------------------------

#include "G4HadronStoppingProcess.hh"
#include "G4HadronicProcessStore.hh"
#include "G4HadronicProcessType.hh"
#include "G4EmCaptureCascade.hh"
#include "G4Nucleus.hh"
#include "G4HadFinalState.hh"
#include "G4HadProjectile.hh"
#include "G4HadSecondary.hh"
#include "G4Material.hh"
#include "G4Element.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadronStoppingProcess::G4HadronStoppingProcess(const G4String& name)
    : G4HadronicProcess(name, fHadronAtRest),
      fElementSelector(new G4ElementSelector()),
      fEmCascade(new G4EmCaptureCascade()),  // Owned by InteractionRegistry
      fBoundDecay(0),
      emcID(-1),
      ncID(-1),
      dioID(-1)
{
  // Modify G4VProcess flags to emulate G4VRest instead of G4VDiscrete
  enableAtRestDoIt = true;
  enablePostStepDoIt = false;

  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4HadronStoppingProcess::~G4HadronStoppingProcess()
{
  //G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
  delete fElementSelector;
  // NOTE: fEmCascade and fEmBoundDecay owned by registry, not locally
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4HadronStoppingProcess::IsApplicable(const G4ParticleDefinition& p)
{
  return (p.GetPDGCharge() < 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4HadronStoppingProcess::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,&p);
  emcID = G4PhysicsModelCatalog::GetModelID(G4String("model_" + (GetProcessName() + "_EMCascade")));
  ncID  = G4PhysicsModelCatalog::GetModelID(G4String("model_" + (GetProcessName() + "_NuclearCapture")));
  dioID = G4PhysicsModelCatalog::GetModelID(G4String("model_" + (GetProcessName() + "_DIO")));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4HadronStoppingProcess::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4HadronStoppingProcess::AtRestGetPhysicalInteractionLength(
    const G4Track&, G4ForceCondition* condition)
{
  *condition = NotForced;
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4HadronStoppingProcess::PostStepGetPhysicalInteractionLength(
    const G4Track&, G4double, G4ForceCondition* condition)
{
  *condition = NotForced;
  return DBL_MAX;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4HadronStoppingProcess::AtRestDoIt(const G4Track& track, 
						       const G4Step&)
{
  // if primary is not Alive then do nothing
  theTotalResult->Initialize(track);

  G4Nucleus* nucleus = GetTargetNucleusPointer();
  const G4Element* elm = fElementSelector->SelectZandA(track, nucleus);

  G4HadFinalState* result = 0;
  thePro.Initialise(track);

  // save track time an dstart capture from zero time
  thePro.SetGlobalTime(0.0);
  G4double time0 = track.GetGlobalTime();

  G4bool nuclearCapture = true;

  // Do the electromagnetic cascade in the nuclear field.
  // EM cascade should keep G4HadFinalState object,
  // because it will not be deleted at the end of this method
  //
  result = fEmCascade->ApplyYourself(thePro, *nucleus);
  G4double ebound = result->GetLocalEnergyDeposit(); 
  G4double edep = 0.0; 
  G4int nSecondaries = (G4int)result->GetNumberOfSecondaries();
  G4int nEmCascadeSec = nSecondaries;

  // Try decay from bound level 
  // For mu- the time of projectile should be changed.
  // Decay should keep G4HadFinalState object,
  // because it will not be deleted at the end of this method.
  //
  thePro.SetBoundEnergy(ebound);
  if(fBoundDecay) {
    G4HadFinalState* resultDecay = 
      fBoundDecay->ApplyYourself(thePro, *nucleus);
    G4int n = (G4int)resultDecay->GetNumberOfSecondaries();
    if(0 < n) {
      nSecondaries += n;
      result->AddSecondaries(resultDecay);
    } 
    if(resultDecay->GetStatusChange() == stopAndKill) {
      nuclearCapture = false; 
    }
    resultDecay->Clear();
  }

  if(nuclearCapture) {
 
    // delay of capture
    G4double capTime = thePro.GetGlobalTime();
    thePro.SetGlobalTime(0.0);

    // select model
    G4HadronicInteraction* model = 0;
    try {
      model = ChooseHadronicInteraction(thePro, *nucleus, 
					track.GetMaterial(), elm);
    }
    catch(G4HadronicException & aE) {
      G4ExceptionDescription ed;
      ed << "Target element "<<elm->GetName()<<"  Z= " 
	 << nucleus->GetZ_asInt() << "  A= " 
	 << nucleus->GetA_asInt() << G4endl;
      DumpState(track,"ChooseHadronicInteraction",ed);
      ed << " No HadronicInteraction found out" << G4endl;
      G4Exception("G4HadronStoppingProcess::AtRestDoIt", "had005", 
		  FatalException, ed);
    }

    G4HadFinalState* resultNuc = 0;
    G4int reentryCount = 0;
    do {
      // sample final state
      // nuclear interaction should keep G4HadFinalState object
      // model should define time of each secondary particle
      try {
	resultNuc = model->ApplyYourself(thePro, *nucleus);
	++reentryCount;
      }
      catch(G4HadronicException & aR) {
	G4ExceptionDescription ed;
	ed << "Call for " << model->GetModelName() << G4endl;
	ed << "Target element "<<elm->GetName()<<"  Z= " 
	   << nucleus->GetZ_asInt() 
	   << "  A= " << nucleus->GetA_asInt() << G4endl;
	DumpState(track,"ApplyYourself",ed);
	ed << " ApplyYourself failed" << G4endl;
	G4Exception("G4HadronStoppingProcess::AtRestDoIt", "had006", 
		    FatalException, ed);
      }

      // Check the result for catastrophic energy non-conservation
      resultNuc = CheckResult(thePro, *nucleus, resultNuc);

      if(reentryCount>100) {
	G4ExceptionDescription ed;
	ed << "Call for " << model->GetModelName() << G4endl;
	ed << "Target element "<<elm->GetName()<<"  Z= " 
	   << nucleus->GetZ_asInt() 
	   << "  A= " << nucleus->GetA_asInt() << G4endl;
	DumpState(track,"ApplyYourself",ed);
	ed << " ApplyYourself does not completed after 100 attempts" << G4endl;
	G4Exception("G4HadronStoppingProcess::AtRestDoIt", "had006", 
		    FatalException, ed);  
      }
      // Loop checking, 06-Aug-2015, Vladimir Ivanchenko
    } while(!resultNuc);

    edep = resultNuc->GetLocalEnergyDeposit();
    std::size_t nnuc = resultNuc->GetNumberOfSecondaries();

    // add delay time of capture
    for(std::size_t i=0; i<nnuc; ++i) { 
      G4HadSecondary* sec = resultNuc->GetSecondary(i);
      sec->SetTime(capTime + sec->GetTime());
    }

    nSecondaries += nnuc;
    result->AddSecondaries(resultNuc); 
    resultNuc->Clear();
  }

  // Fill results
  //
  theTotalResult->ProposeTrackStatus(fStopAndKill);
  theTotalResult->ProposeLocalEnergyDeposit(edep);  
  theTotalResult->SetNumberOfSecondaries(nSecondaries);
  G4double w  = track.GetWeight();
  theTotalResult->ProposeWeight(w);
  for(G4int i=0; i<nSecondaries; ++i) {
    G4HadSecondary* sec = result->GetSecondary(i);

    // add track global time to the reaction time
    G4double time = sec->GetTime();
    if(time < 0.0) { time = 0.0; }
    time += time0;

    // create secondary track
    G4Track* t = new G4Track(sec->GetParticle(),
			     time, 
			     track.GetPosition());
    t->SetWeight(w*sec->GetWeight());

    // use SetCreatorModelID to "label" the track
    if (i<nEmCascadeSec) {
      t->SetCreatorModelID(emcID);
    } else if (nuclearCapture) {
      t->SetCreatorModelID(ncID);
    } else {
      t->SetCreatorModelID(dioID);
    }

    t->SetTouchableHandle(track.GetTouchableHandle());
    theTotalResult->AddSecondary(t);
  }
  result->Clear();

  if (epReportLevel != 0) { 
    CheckEnergyMomentumConservation(track, *nucleus);
  }
  return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4HadronStoppingProcess::ProcessDescription(std::ostream& outFile) const
{
  outFile << "Base process for negatively charged particle capture at rest.\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
