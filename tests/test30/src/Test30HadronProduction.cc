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
                       const G4Track& aTrack,
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
  if(gen) {
    G4cout << "Test30HadronProduction: the Secondary Generator <"
           << theGenerator->GeneratorName()
           << "> is initialized"
           << G4endl;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4VParticleChange* Test30HadronProduction::PostStepDoIt(
                    const G4Track& track, 
                    const G4Step& step)
{
  pParticleChange = theGenerator->Secondaries(track);
  ClearNumberOfInteractionLengthLeft();

  if(pParticleChange->GetStatusChange() != fStopAndKill) {

    theChange.Initialize(track);
    G4int ns = pParticleChange->GetNumberOfSecondaries();
    ns++;
    theChange.SetNumberOfSecondaries(ns);
    for(G4int i=0; i<ns-1; i++) {
      theChange.AddSecondary(pParticleChange->GetSecondary(i));
    }
    G4ParticleChange* pc = dynamic_cast<G4ParticleChange*>(pParticleChange);
    G4Track* tr = new G4Track(track);
    tr->SetKineticEnergy(pc->GetEnergyChange());
    tr->SetMomentumDirection(*(pc->GetMomentumDirectionChange()));
    theChange.AddSecondary(tr);
    theChange.SetStatusChange(fStopAndKill);
    pParticleChange->Clear();
    pParticleChange = &theChange;     
  }

  return pParticleChange;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
