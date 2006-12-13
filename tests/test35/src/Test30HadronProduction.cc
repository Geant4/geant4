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
// -------------------------------------------------------------
//      GEANT 4 class 
//
//      ---------- Test30HadronProduction -------
//                by Vladimir Ivanchenko, 12 March 2002 
// 
//    Modified:
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#include "Test30HadronProduction.hh"
#include "Test30VSecondaryGenerator.hh"

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30HadronProduction::Test30HadronProduction(const G4String& aName)
 :G4VDiscreteProcess(aName),
  theGenerator(0)
{
  InitializeMe();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test30HadronProduction::InitializeMe()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Test30HadronProduction::~Test30HadronProduction()
{
  if(theGenerator) delete theGenerator;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4double Test30HadronProduction::PostStepGetPhysicalInteractionLength(
                       const G4Track&,
                             G4double,
                             G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double z = DBL_MAX;

  return z;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Test30HadronProduction::SetSecondaryGenerator(Test30VSecondaryGenerator* gen)
{
  if(theGenerator) delete theGenerator;
  theGenerator = gen;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* Test30HadronProduction::PostStepDoIt(
                    const G4Track& track,
                    const G4Step&)
{
  G4HadFinalState* result = theGenerator->Secondaries(track);
  ClearNumberOfInteractionLengthLeft();

  theChange.Initialize(track);

  G4int ns = result->GetNumberOfSecondaries();
  G4int nb = ns;
  if(result->GetStatusChange() == isAlive) nb++;
  
  theChange.ProposeTrackStatus(fStopAndKill);
  theChange.SetNumberOfSecondaries(nb);

  for(G4int i=0; i<ns; i++) {
    G4Track* tr = new G4Track(result->GetSecondary(i)->GetParticle(),
                              track.GetGlobalTime(),
	                      track.GetPosition());
    theChange.AddSecondary(tr);
  }

  if(result->GetStatusChange() == isAlive) {
    G4DynamicParticle* dp = new G4DynamicParticle(*(track.GetDynamicParticle()));
    G4Track* tr = new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
    tr->SetKineticEnergy(result->GetEnergyChange());
    tr->SetMomentumDirection(result->GetMomentumChange());
    theChange.AddSecondary(tr);
  }
  result->Clear();

  return &theChange;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
