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
//---------------------------------------------------------------------

#include "G4PiMinusAbsorptionBertini.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleTypes.hh"
#include "Randomize.hh" 
#include "G4HadronicProcessStore.hh"
#include <string.h>
#include <cmath>
#include <stdio.h>


// constructor

G4PiMinusAbsorptionBertini::G4PiMinusAbsorptionBertini(const G4String& processName,
                                                       G4ProcessType   aType ) :
G4VRestProcess (processName, aType),
pdefPionMinus(G4PionMinus::PionMinus())

{
    if (verboseLevel>0) {
        G4cout << GetProcessName() << " is created "<< G4endl;
    }
    SetProcessSubType(fHadronAtRest);
    
    cascade = new G4CascadeInterface;
    cascade->usePreCompoundDeexcitation();
    
    G4HadronicProcessStore::Instance()->RegisterExtraProcess(this);
}

// destructor

G4PiMinusAbsorptionBertini::~G4PiMinusAbsorptionBertini()
{
    G4HadronicProcessStore::Instance()->DeRegisterExtraProcess(this);
}

void G4PiMinusAbsorptionBertini::PreparePhysicsTable(const G4ParticleDefinition& p) 
{
    G4HadronicProcessStore::Instance()->RegisterParticleForExtraProcess(this, &p);
}

void G4PiMinusAbsorptionBertini::BuildPhysicsTable(const G4ParticleDefinition& p) 
{
    G4HadronicProcessStore::Instance()->PrintInfo(&p);
}

// methods.............................................................................

G4bool G4PiMinusAbsorptionBertini::IsApplicable(const G4ParticleDefinition& particle)
{
    return ( &particle == pdefPionMinus );
}


G4VParticleChange* G4PiMinusAbsorptionBertini::AtRestDoIt(const G4Track& track,
                                                          const G4Step&)
//
// Handles PionMinuss at rest
//
{
    // We construct a fake track and we set the kinetic energy to 1keV in order to be able to use Bertini+PrecompoundDeexciation
    
    G4Track faketrack(track);
    faketrack.SetStep(track.GetStep());
    faketrack.SetKineticEnergy(1*keV);
    
    
    //   Initialize ParticleChange
    //     all members of G4VParticleChange are set to equal to 
    //     corresponding member in G4Track
    
    aParticleChange.Initialize(track);
    
    //   Store some global quantities that depend on current material and particle
    
    G4Material * aMaterial = track.GetMaterial();
    
    const G4int numberOfElements = aMaterial->GetNumberOfElements();
    const G4ElementVector* theElementVector = aMaterial->GetElementVector();
    
    const G4double* theAtomicNumberDensity = aMaterial->GetAtomicNumDensityVector();
    G4double normalization = 0;
    for ( G4int i1=0; i1 < numberOfElements; i1++ )
    {
        normalization += theAtomicNumberDensity[i1] ; // change when nucleon specific
        // probabilities are included.
    }
    G4double runningSum= 0.;
    G4double random = G4UniformRand()*normalization;
    for ( G4int i2=0; i2 < numberOfElements; i2++ )
    {
        runningSum += theAtomicNumberDensity[i2]; // change when nucleon specific
        // probabilities are included.
        if (random<=runningSum)
        {
            targetZ = G4double((*theElementVector)[i2]->GetZ());
            currentN = (*theElementVector)[i2]->GetN();
        }
    }
    if (random>runningSum)
    {
        targetZ = G4double((*theElementVector)[numberOfElements-1]->GetZ());
        currentN = (*theElementVector)[numberOfElements-1]->GetN();
        
    }
    
    targetNucleus.SetParameters(currentN, targetZ);
    
    if (verboseLevel>1) {
        G4cout << "G4PiMinusAbsorptionBertini::AtRestDoIt is invoked " <<G4endl;
    }
    
    G4HadFinalState* result = cascade->ApplyYourself(faketrack, targetNucleus);
    
    ClearNumberOfInteractionLengthLeft();
    
    
    G4int ns = result->GetNumberOfSecondaries();
    G4int nb = ns;
    if(result->GetStatusChange() == isAlive) nb++;
    
    aParticleChange.ProposeTrackStatus(fStopAndKill);
    aParticleChange.SetNumberOfSecondaries(nb);
    
    for(G4int i=0; i<ns; i++) {
        G4Track* tr = new G4Track(result->GetSecondary(i)->GetParticle(),
                                  track.GetGlobalTime(),
                                  track.GetPosition());
        aParticleChange.AddSecondary(tr);
    }
    
    if(result->GetStatusChange() == isAlive) {
        G4DynamicParticle* dp = new G4DynamicParticle(*(track.GetDynamicParticle()));
        G4Track* tr = new G4Track(dp,track.GetGlobalTime(),track.GetPosition());
        tr->SetKineticEnergy(result->GetEnergyChange());
        tr->SetMomentumDirection(result->GetMomentumChange());
        aParticleChange.AddSecondary(tr);
    }
    result->Clear();
    
    
    return &aParticleChange;
    
}


