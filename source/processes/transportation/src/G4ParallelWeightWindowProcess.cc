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
//
// $Id: G4ParallelWeightWindowProcess.cc,v 1.2 2002-05-31 10:16:02 dressel Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ParallelWeightWindowProcess.cc
//
// ----------------------------------------------------------------------

#include "G4ParallelWeightWindowProcess.hh"
#include "G4VIStore.hh"
#include "G4VParallelStepper.hh"
#include "G4VWeightWindowAlgorithm.hh"

#include "G4Pstring.hh"

G4ParallelWeightWindowProcess::G4ParallelWeightWindowProcess(G4VIStore &aIstore,
					       G4VParallelStepper &aStepper,
					       G4VWeightWindowAlgorithm &aWWAlgorithm,
					       const G4String &aName)
  : G4VProcess(aName), 
  fIStore(aIstore),
  fPStepper(aStepper),
  fWWAlgorithm(aWWAlgorithm),
  fNsplit_Weight(G4Nsplit_Weight(0,0)),
  fInitStep(false)
{
  fParticleChange = new G4ParticleChange;
  G4VProcess::pParticleChange = fParticleChange;
}

G4ParallelWeightWindowProcess::~G4ParallelWeightWindowProcess()
{
  delete fParticleChange;
}

void G4ParallelWeightWindowProcess::StartTracking(){
  fInitStep = true;
}
void G4ParallelWeightWindowProcess::EndTracking(){}

G4double 
G4ParallelWeightWindowProcess::
PostStepGetPhysicalInteractionLength(const G4Track& aTrack,
				     G4double   previousStepSize,
				     G4ForceCondition* condition)
{
  *condition = NotForced;
  G4double stepLength = kInfinity;
  
  if (fInitStep && aTrack.GetTrackStatus() != fStopAndKill) {
    fInitStep = false;
    G4double importance = fIStore.GetImportance(fPStepper.GetPStep().fPostTouchableKey);
    G4double weight = aTrack.GetWeight();
    
    fNsplit_Weight = fWWAlgorithm.Calculate(weight, importance);
    
    if (aTrack.GetWeight() != fNsplit_Weight.fW || 
	fNsplit_Weight.fN != 1) {
      *condition = ExclusivelyForced;
    }
  }
  
  return stepLength;
}

G4VParticleChange * 
G4ParallelWeightWindowProcess::PostStepDoIt(const G4Track& aTrack,
				  const G4Step& aStep)
{
  fParticleChange->Initialize(aTrack);

  fImportancePostStepDoIt.DoIt(aTrack, fParticleChange, fNsplit_Weight);
  
  return fParticleChange;
}
    
void G4ParallelWeightWindowProcess::Error(const G4String &m)
{
  G4cout << "ERROR - G4ParallelWeightWindowProcess::" << m << G4endl;
  G4Exception("Program aborted.");
}

void G4ParallelWeightWindowProcess::Warning(const G4String &m)
{
  G4cout << "WARNING - G4ParallelWeightWindowProcess: " << G4endl;
  G4cout << m << G4endl;
}
