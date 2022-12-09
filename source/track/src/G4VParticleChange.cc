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
// G4VParticleChange class implementation
//
// Author: Hisaya Kurashige, 23 March 1998
// --------------------------------------------------------------------

#include "G4VParticleChange.hh"
#include "G4SystemOfUnits.hh"
#include "G4ExceptionSeverity.hh"

const G4double G4VParticleChange::accuracyForWarning   = 1.0e-9;
const G4double G4VParticleChange::accuracyForException = 0.001;
const G4int G4VParticleChange::maxError = 10;

// --------------------------------------------------------------------
G4VParticleChange::G4VParticleChange()
{
#ifdef G4VERBOSE
  // activate CheckIt if in VERBOSE mode
  debugFlag = true;
#endif
}

// --------------------------------------------------------------------
void G4VParticleChange::AddSecondary(G4Track* aTrack)
{
  if(debugFlag)
    CheckSecondary(*aTrack);

  if(!fSetSecondaryWeightByProcess)
    aTrack->SetWeight(theParentWeight);

  // add a secondary after size check
  if(theSizeOftheListOfSecondaries > theNumberOfSecondaries)
  {
    theListOfSecondaries[theNumberOfSecondaries] = aTrack;
  }
  else
  {
    theListOfSecondaries.push_back(aTrack);
    ++theSizeOftheListOfSecondaries;
  }
  ++theNumberOfSecondaries;
}

// --------------------------------------------------------------------
G4Step* G4VParticleChange::UpdateStepInfo(G4Step* pStep)
{
  // Update the G4Step specific attributes
  pStep->SetStepLength(theTrueStepLength);
  pStep->AddTotalEnergyDeposit(theLocalEnergyDeposit);
  pStep->AddNonIonizingEnergyDeposit(theNonIonizingEnergyDeposit);
  pStep->SetControlFlag(theSteppingControlFlag);

  if(theFirstStepInVolume)
  {
    pStep->SetFirstStepFlag();
  }
  else
  {
    pStep->ClearFirstStepFlag();
  }
  if(theLastStepInVolume)
  {
    pStep->SetLastStepFlag();
  }
  else
  {
    pStep->ClearLastStepFlag();
  }

  return pStep;
}

// --------------------------------------------------------------------
G4Step* G4VParticleChange::UpdateStepForAtRest(G4Step* Step)
{
  if(isParentWeightProposed)
  {
    Step->GetPostStepPoint()->SetWeight(theParentWeight);
  }
  return UpdateStepInfo(Step);
}

// --------------------------------------------------------------------
G4Step* G4VParticleChange::UpdateStepForAlongStep(G4Step* Step)
{
  if(isParentWeightProposed)
  {
    G4double initialWeight = Step->GetPreStepPoint()->GetWeight();
    G4double currentWeight = Step->GetPostStepPoint()->GetWeight();
    G4double finalWeight   = (theParentWeight / initialWeight) * currentWeight;
    Step->GetPostStepPoint()->SetWeight(finalWeight);
  }
  return UpdateStepInfo(Step);
}

// --------------------------------------------------------------------
G4Step* G4VParticleChange::UpdateStepForPostStep(G4Step* Step)
{
  if(isParentWeightProposed)
  {
    Step->GetPostStepPoint()->SetWeight(theParentWeight);
  }
  return UpdateStepInfo(Step);
}

// --------------------------------------------------------------------
void G4VParticleChange::DumpInfo() const
{
  auto vol = theCurrentTrack->GetVolume();
  G4String vname = (nullptr == vol) ? "" : vol->GetName();
  G4long olprc = G4cout.precision(8);
  G4cout << "      -----------------------------------------------" << G4endl;
  G4cout << "        G4VParticleChange Information " << G4endl;
  G4cout << "        TrackID             : " << theCurrentTrack->GetTrackID()
	 << G4endl;
  G4cout << "        ParentID            : " << theCurrentTrack->GetParentID()
	 << G4endl;
  G4cout << "        Particle            : " 
	 << theCurrentTrack->GetParticleDefinition()->GetParticleName()
         << G4endl;
  G4cout << "        Kinetic energy (MeV): " 
         << theCurrentTrack->GetKineticEnergy() << G4endl;
  G4cout << "        Position (mm)       : " 
         << theCurrentTrack->GetPosition() << G4endl;
  G4cout << "        Direction           : "
         << theCurrentTrack->GetMomentumDirection() << G4endl;
  G4cout << "        PhysicsVolume       : " << vname << G4endl;
  G4cout << "        Material            : " 
	 << theCurrentTrack->GetMaterial()->GetName() << G4endl;
  G4cout << "      -----------------------------------------------" << G4endl;

  G4cout << "        # of secondaries    : " << std::setw(20)
         << theNumberOfSecondaries << G4endl;

  G4cout << "      -----------------------------------------------" << G4endl;

  G4cout << "        Energy Deposit (MeV): " << std::setw(20)
         << theLocalEnergyDeposit / MeV << G4endl;

  G4cout << "   NIEL Energy Deposit (MeV): " << std::setw(20)
         << theNonIonizingEnergyDeposit / MeV << G4endl;

  G4cout << "        Track Status        : " << std::setw(20);
  if(theStatusChange == fAlive)
  {
    G4cout << " Alive";
  }
  else if(theStatusChange == fStopButAlive)
  {
    G4cout << " StopButAlive";
  }
  else if(theStatusChange == fStopAndKill)
  {
    G4cout << " StopAndKill";
  }
  else if(theStatusChange == fKillTrackAndSecondaries)
  {
    G4cout << " KillTrackAndSecondaries";
  }
  else if(theStatusChange == fSuspend)
  {
    G4cout << " Suspend";
  }
  else if(theStatusChange == fPostponeToNextEvent)
  {
    G4cout << " PostponeToNextEvent";
  }
  G4cout << G4endl;
  G4cout << "        TruePathLength (mm) : " << std::setw(20)
         << theTrueStepLength / mm << G4endl;
  G4cout << "        Stepping Control    : " << std::setw(20)
         << theSteppingControlFlag << G4endl;
  if(theFirstStepInVolume)
  {
    G4cout << "       First step in volume" << G4endl;
  }
  if(theLastStepInVolume)
  {
    G4cout << "       Last step in volume" << G4endl;
  }

#ifdef G4VERBOSE
  if(nError == maxError)
  {
    G4cout << "      -----------------------------------------------" << G4endl;
    G4cout << "        G4VParticleChange warnings closed " << G4endl;
    G4cout << "      -----------------------------------------------" << G4endl;
  }
#endif
  
  G4cout.precision(olprc);
}

// --------------------------------------------------------------------
G4bool G4VParticleChange::CheckIt([[maybe_unused]] const G4Track& aTrack)
{
  G4bool isOK = true;

  // Energy deposit should not be negative
  if(theLocalEnergyDeposit < 0.0)
  {
    isOK = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4cout << "  G4VParticleChange::CheckIt : ";
      G4cout << "the energy deposit " << theLocalEnergyDeposit/MeV
	     << " MeV is negative !!" << G4endl;
    }
#endif
    theLocalEnergyDeposit = 0.0;
  }

  // true path length should not be negative
  if(theTrueStepLength < 0.0)
  {
    isOK = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4cout << "  G4VParticleChange::CheckIt : ";
      G4cout << "true path length " << theTrueStepLength/mm
	     << " mm is negative !!" << G4endl;
    }
#endif
    theTrueStepLength = (1.e-12) * mm;
  }

  if(!isOK)
  {
    if(nError < maxError)
    {
#ifdef G4VERBOSE
      // dump out information of this particle change
      DumpInfo();
#endif
      G4Exception("G4VParticleChange::CheckIt()", "TRACK001", JustWarning,
		  "Step length and/or energy deposit are illegal");
    }
  }
  return isOK;
}

// --------------------------------------------------------------------
G4bool G4VParticleChange::CheckSecondary(G4Track& aTrack)
{
  G4bool isOK = true;

  // MomentumDirection should be unit vector
  G4double ekin = aTrack.GetKineticEnergy();
  auto dir = aTrack.GetMomentumDirection();
  G4double accuracy = std::abs(dir.mag2() - 1.0);
  if(accuracy > accuracyForWarning)
  {
    isOK = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4String mname = aTrack.GetCreatorModelName();
      G4cout << " G4VParticleChange::CheckSecondary : " << G4endl;
      G4cout << " the momentum direction " << dir
	     << " is not unit vector !!" << G4endl;
      G4cout << " Difference=" << accuracy 
	     << " Ekin(MeV)=" << ekin/MeV 
	     << "  " << aTrack.GetParticleDefinition()->GetParticleName()
	     << " created by " << mname
             << G4endl;
    }
#endif
    aTrack.SetMomentumDirection(dir.unit());
  }

  // Kinetic Energy should not be negative
  if(ekin < 0.0)
  {
    isOK = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4String mname = aTrack.GetCreatorModelName();
      G4cout << " G4VParticleChange::CheckSecondary : " << G4endl;
      G4cout << " Ekin(MeV)=" << ekin << " is negative !!  "
	     << aTrack.GetParticleDefinition()->GetParticleName()
	     << " created by " << mname
             << G4endl;
    }
#endif
    aTrack.SetKineticEnergy(0.0);
  }

  // Check timing of secondaries
  G4double time = aTrack.GetGlobalTime();
  if(time < theParentGlobalTime)
  {
    isOK = false;
    ++nError;
#ifdef G4VERBOSE
    if(nError < maxError)
    {
      G4String mname = aTrack.GetCreatorModelName();
      G4cout << " G4VParticleChange::CheckSecondary : " << G4endl;
      G4cout << " The global time of secondary goes back compared to the parent !!" << G4endl;
      G4cout << " ParentTime(ns)=" << theParentGlobalTime/ns
	     << " SecondaryTime(ns)= " << time/ns
	     << " Difference(ns)=" << (theParentGlobalTime - time)/ns
	     << G4endl;
      G4cout << " Ekin(MeV)=" << ekin
	     << aTrack.GetParticleDefinition()->GetParticleName()
	     << " created by " << mname << G4endl;
    }
#endif
    aTrack.SetGlobalTime(theParentGlobalTime);
  }

  // Exit with error
  if(!isOK)
  {
    if(nError < maxError)
    {
#ifdef G4VERBOSE
      DumpInfo();
#endif
      G4Exception("G4VParticleChange::CheckSecondary()", "TRACK001",
		  JustWarning, "Secondary with illegal time and/or energy and/or momentum");
    }
  }
  return isOK;
}

// --------------------------------------------------------------------
G4double G4VParticleChange::GetAccuracyForWarning() const
{
  return accuracyForWarning;
}

// --------------------------------------------------------------------
G4double G4VParticleChange::GetAccuracyForException() const
{
  return accuracyForException;
}

// --------------------------------------------------------------------
// Obsolete methods for parent weight
//
void G4VParticleChange::SetParentWeightByProcess(G4bool) {}
G4bool G4VParticleChange::IsParentWeightSetByProcess() const { return true; }
