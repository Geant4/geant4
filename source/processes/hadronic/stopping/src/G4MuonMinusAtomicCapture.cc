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
// GEANT4 Class header file
//
// File name:     G4MuonMinusAtomicCapture
//
// 20160912 K.L. Genser - New process using G4MuonicAtom somewhat
//                        based on G4HadronStoppingProcess
//
// Class Description: 
//
// Stopping of mu-
//
// G4VParticleChange will contain gammas from G4EmCaptureCascade and
// resulting G4MuonicAtom
//
//
//------------------------------------------------------------------------

#include "G4MuonMinusAtomicCapture.hh"
#include "G4ParticleDefinition.hh"
#include "G4HadronicProcessType.hh"
#include "G4MuonMinusBoundDecay.hh"
#include "G4HadronicInteraction.hh"
#include "G4HadProjectile.hh"
#include "G4HadronicProcessStore.hh"
#include "G4EmCaptureCascade.hh"
#include "G4MuonMinus.hh"
#include "G4IonTable.hh"
#include "G4RandomDirection.hh"
#include "G4HadSecondary.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusAtomicCapture::G4MuonMinusAtomicCapture(const G4String& name)
  : G4VRestProcess(name, fHadronic),
    fElementSelector(new G4ElementSelector()),
    fEmCascade(new G4EmCaptureCascade()),  // Owned by InteractionRegistry
    theTotalResult(new G4ParticleChange()),
    result(nullptr)
{
  SetProcessSubType(fMuAtomicCapture);
  G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4MuonMinusAtomicCapture::~G4MuonMinusAtomicCapture()
{
  G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
  delete theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool G4MuonMinusAtomicCapture::IsApplicable(const G4ParticleDefinition& p)
{
  return (&p == G4MuonMinus::MuonMinus());
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void 
G4MuonMinusAtomicCapture::PreparePhysicsTable(const G4ParticleDefinition& p)
{
  G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this,&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusAtomicCapture::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
  G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double G4MuonMinusAtomicCapture::AtRestGetPhysicalInteractionLength(
                                                              const G4Track&, G4ForceCondition* condition)
{
  *condition = NotForced;
  return 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* G4MuonMinusAtomicCapture::AtRestDoIt(const G4Track& track, 
                                                const G4Step&)
{
  // if primary is not Alive then do nothing (how?)
  theTotalResult->Initialize(track);

  G4Nucleus* nucleus = &targetNucleus;
  // the call below actually sets the nucleus params; 
  // G4Nucleus targetNucleus; is a member of G4HadronicProcess
  // G4Element* elm = 
  fElementSelector->SelectZandA(track, nucleus);

  thePro.Initialise(track); // thePro was G4HadProjectile from G4HadronicProcess

  // save track time an dstart capture from zero time
  thePro.SetGlobalTime(0.0);
  G4double time0 = track.GetGlobalTime();

  // Do the electromagnetic cascade in the nuclear field.
  // EM cascade should keep G4HadFinalState object,
  // because it will not be deleted at the end of this method
  //
  result = fEmCascade->ApplyYourself(thePro, *nucleus);
  G4double ebound = result->GetLocalEnergyDeposit(); // may need to carry this over; review
  G4double edep = 0.0;
  G4int nSecondaries = (G4int)result->GetNumberOfSecondaries();
  thePro.SetBoundEnergy(ebound);

  // creating the muonic atom
  ++nSecondaries;

  G4IonTable* itp = G4IonTable::GetIonTable();
  G4ParticleDefinition* muonicAtom = itp->GetMuonicAtom(nucleus->GetZ_asInt(),
                                                        nucleus->GetA_asInt());

  G4DynamicParticle* dp = new G4DynamicParticle(muonicAtom,G4RandomDirection(),0.);
  G4HadSecondary hadSec(dp);
  hadSec.SetTime(time0);
  result->AddSecondary(hadSec);

  // Fill results
  //
  theTotalResult->ProposeTrackStatus(fStopAndKill);
  theTotalResult->ProposeLocalEnergyDeposit(edep);  
  theTotalResult->SetNumberOfSecondaries(nSecondaries);
  G4double w  = track.GetWeight();
  theTotalResult->ProposeWeight(w);

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << __func__
             << " nSecondaries "
             << nSecondaries
             << G4endl;
    }
#endif

  for(G4int i=0; i<nSecondaries; ++i) {
    G4HadSecondary* sec = result->GetSecondary(i);

    // add track global time to the reaction time
    G4double time = sec->GetTime();
    if(time < 0.0) { time = 0.0; }
    time += time0;

#ifdef G4VERBOSE
    if (GetVerboseLevel() > 1) {
      G4cout << __func__
             << " "
             << i
             << " Resulting secondary "
             << sec->GetParticle()->GetPDGcode()
             << " "
             << sec->GetParticle()->GetDefinition()->GetParticleName()
             << G4endl;
    }
#endif

    // create secondary track
    G4Track* t = new G4Track(sec->GetParticle(),
			     time, 
			     track.GetPosition());
    t->SetWeight(w*sec->GetWeight());

    t->SetTouchableHandle(track.GetTouchableHandle());
    theTotalResult->AddSecondary(t);
  }
  result->Clear();

  // fixme: needs to be done at the MuonicAtom level
  // if (epReportLevel != 0) { // G4HadronicProcess::
  //   CheckEnergyMomentumConservation(track, *nucleus);
  // }
  return theTotalResult;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void G4MuonMinusAtomicCapture::ProcessDescription(std::ostream& outFile) const
{
  outFile << "Stopping of mu- using default element selector, EM cascade"
          << "G4MuonicAtom is created\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
