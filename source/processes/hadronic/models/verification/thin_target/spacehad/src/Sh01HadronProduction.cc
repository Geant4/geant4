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
//   Sh01HadronProduction 
//      
// 
//    Modified:
//
//    06.03.03 V. Grichine (based on test30 of V. Ivanchenko)


#include "Sh01HadronProduction.hh"
#include "Sh01SecondaryGenerator.hh"

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

////////////////////////////////////////////////////////////////////////////

Sh01HadronProduction::Sh01HadronProduction(const G4String& aName)  
 :G4VDiscreteProcess(aName),
  theGenerator(0)
{
  InitializeMe();
}

//////////////////////////////////////////////////////////////////////////

void Sh01HadronProduction::InitializeMe()
{}

//////////////////////////////////////////////////////////////////////////

Sh01HadronProduction::~Sh01HadronProduction()
{
  if(theGenerator) delete theGenerator;
}

/////////////////////////////////////////////////////////////////////////

G4double Sh01HadronProduction::PostStepGetPhysicalInteractionLength(
                       const G4Track& aTrack,
                             G4double,
                             G4ForceCondition* condition)
{
  // condition is set to "Not Forced"
  *condition = NotForced;
  G4double z = DBL_MAX;

  return z;
}
    
/////////////////////////////////////////////////////////////////////////

void Sh01HadronProduction::SetSecondaryGenerator(Sh01SecondaryGenerator* gen)
{
  if(theGenerator) delete theGenerator;
  theGenerator = gen;
}

/////////////////////////////////////////////////////////////////////////

G4VParticleChange* Sh01HadronProduction::PostStepDoIt(
                    const G4Track& track, 
                    const G4Step& step)
{
  // pParticleChange = theGenerator->Secondaries(track);


  G4HadFinalState* result = theGenerator->Secondaries(track);
  ClearNumberOfInteractionLengthLeft();


  theChange.Initialize(track);

  G4int ns = result->GetNumberOfSecondaries();
  G4int nb = ns;
  if(result->GetStatusChange() == isAlive) nb++;

  theChange.SetStatusChange(fStopAndKill);
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
 
  /*
  if(pParticleChange->GetStatusChange() != fStopAndKill) 
  {

    theChange.Initialize(track);
    G4int ns = pParticleChange->GetNumberOfSecondaries();
    ns++;
    theChange.SetNumberOfSecondaries(ns);

    for(G4int i=0; i<ns-1; i++) 
    {
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
  */
}


//////////////////////////////////////////////////////////////////////////





