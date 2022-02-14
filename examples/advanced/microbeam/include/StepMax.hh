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
// This example is provided by the Geant4-DNA collaboration
// Any report or published results obtained using the Geant4-DNA software 
// shall cite the following Geant4-DNA collaboration publication:
// Med. Phys. 37 (2010) 4692-4708
// The Geant4-DNA web site is available at http://geant4-dna.org
// 
// If you use this example, please cite the following publication:
// Rad. Prot. Dos. 133 (2009) 2-11
//
#ifndef StepMax_h
#define StepMax_h 1

#include "globals.hh"
#include "G4VEmProcess.hh"
#include "G4ParticleDefinition.hh"
#include "G4Step.hh"

class PhysicsListMessenger;

class StepMax : public G4VEmProcess
{
public:

  explicit StepMax(PhysicsListMessenger* mess);
  ~StepMax() override;

  G4bool IsApplicable(const G4ParticleDefinition&) override;

  void PreparePhysicsTable(const G4ParticleDefinition&)override;

  void BuildPhysicsTable(const G4ParticleDefinition&)override;

  void InitialiseProcess(const G4ParticleDefinition*)override;

  G4double PostStepGetPhysicalInteractionLength(const G4Track& track,
                                                        G4double previousStep,
                                                        G4ForceCondition* cond)override;

  G4VParticleChange* PostStepDoIt(const G4Track&, const G4Step&) override;

private:

  PhysicsListMessenger* fMessenger;

  G4double fMaxChargedStep;
  G4bool isInitialised;
};
#endif

