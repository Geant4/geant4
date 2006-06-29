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
#include "Tst33WeightChangeProcess.hh"
#include "Randomize.hh"

Tst33WeightChangeProcess::Tst33WeightChangeProcess()
 : 
  G4VProcess("Tst33WeightChangeProcess"),
  aParticleChange(new G4ParticleChange)
{
  G4VProcess::pParticleChange = aParticleChange;
}

Tst33WeightChangeProcess::~Tst33WeightChangeProcess()
{
  delete aParticleChange;
  aParticleChange = 0;
  G4VProcess::pParticleChange = 0;
}

G4double Tst33WeightChangeProcess::
PostStepGetPhysicalInteractionLength(const G4Track&,
				     G4double  ,
				     G4ForceCondition* condition)
{
  *condition = Forced;
  return kInfinity;
}
  
G4VParticleChange * 
Tst33WeightChangeProcess::PostStepDoIt(const G4Track& aTrack, const G4Step &)
{
  aParticleChange->Initialize(aTrack);
  G4double w = aTrack.GetWeight();

  G4double actionProb = 0.1;
  G4double lowestRelativeWeight = 0.1;
  G4double relativeWeightRange = 1 - lowestRelativeWeight;

  if (G4UniformRand() < actionProb) {
    G4double newWeight = w * (lowestRelativeWeight + 
			      relativeWeightRange * G4UniformRand());
    
    aParticleChange->ProposeWeight(newWeight);
  }
  return aParticleChange;
}

const G4String &Tst33WeightChangeProcess::GetName() const {
  return theProcessName;
}


G4double Tst33WeightChangeProcess::
AlongStepGetPhysicalInteractionLength(const G4Track&,
				      G4double  ,
				      G4double  ,
				      G4double& ,
				      G4GPILSelection*) {
  return -1.0;
}

G4double Tst33WeightChangeProcess::
AtRestGetPhysicalInteractionLength(const G4Track&,
				   G4ForceCondition*) {
  return -1.0;
}

G4VParticleChange* Tst33WeightChangeProcess::AtRestDoIt(const G4Track&,
					       const G4Step&) {
  return 0;
}

G4VParticleChange* Tst33WeightChangeProcess::AlongStepDoIt(const G4Track&,
						  const G4Step&) {
  return 0;
}

