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

#include "Test29HadronProduction.hh"
#include "Test29VSecondaryGenerator.hh"

#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4DynamicParticle.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"
#include "G4Tubs.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4ParticleChange.hh"

// Initializes the Hadronic Process by name (Only G4VDiscreteProcess, not G4VRestProcess)
Test29HadronProduction::Test29HadronProduction(const G4String& aName)
 :G4VDiscreteProcess(aName), theGenerator(0) { InitializeMe(); }

// Model Navigation tool for initialization of the particular model (?)
void Test29HadronProduction::InitializeMe() {}; // M.K.? What for is this member function?

// the G4VParticleChange must be cleaned up by caller! - Bad G4 style (?)
Test29HadronProduction::~Test29HadronProduction() {if(theGenerator) delete theGenerator;}

// Trivial definition of the PostStep InteractionLength definition
G4double Test29HadronProduction::PostStepGetPhysicalInteractionLength
                                 (const G4Track&, G4double, G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  return DBL_MAX;
}

// Set the secondary generator for the Hadronic Process (not used for ElectroWeak Process)
void Test29HadronProduction::SetSecondaryGenerator(Test29VSecondaryGenerator* gen)
{
  if(theGenerator) delete theGenerator;
  theGenerator = gen;
}

// Converts Result of Hadronic Processes to G4VParticleChange form (not for ElectroWeak)
G4VParticleChange* Test29HadronProduction::PostStepDoIt(const G4Track& track,const G4Step&)
{
  G4HadFinalState* result = theGenerator->Secondaries(track);
  ClearNumberOfInteractionLengthLeft();

  theChange.Initialize(track);

  G4int ns = result->GetNumberOfSecondaries();
  G4int nb = ns;
  if(result->GetStatusChange() == isAlive) nb++;
  
  theChange.SetStatusChange(fStopAndKill);
  theChange.SetNumberOfSecondaries(nb);

  for(G4int i=0; i<ns; i++)
  {
    G4Track* tr = new G4Track(result->GetSecondary(i)->GetParticle(),track.GetGlobalTime(),
	                          track.GetPosition());
    theChange.AddSecondary(tr);
  }

  if(result->GetStatusChange() == isAlive)
  {
    G4DynamicParticle* dp = new G4DynamicParticle(*(track.GetDynamicParticle()));
    G4Track* tr = new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
    tr->SetKineticEnergy(result->GetEnergyChange());
    tr->SetMomentumDirection(result->GetMomentumChange());
    theChange.AddSecondary(tr);
  }
  result->Clear();
  return &theChange;
}
