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
// $Id: G4ProcessTest.cc,v 1.9 2001-11-12 09:14:50 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 07 Oct 2001   MGP        Created
//
// -------------------------------------------------------------------
// Class description:
// Test DoIt method of electromagnetic physics processes
// Further documentation available from http://www.ge.infn.it/geant4/lowE/index.html

// -------------------------------------------------------------------

#include "globals.hh"
#include "G4ProcessTest.hh"
#include "G4VProcess.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eIonisation.hh"
#include "G4eplusAnnihilation.hh"
#include "SystemOfUnits.h"

#include "G4ProcessTestAnalysis.hh"

G4ProcessTest::G4ProcessTest() : 
  process(0), ioni(0), brem(0), eProcessManager(0), gProcessManager(0), def(0)
{ }

G4ProcessTest:: ~G4ProcessTest()
{
  delete process;
  process = 0;
  delete ioni;
  ioni = 0;
  delete brem;
  brem = 0;
  delete eProcessManager;
  eProcessManager = 0;
  delete gProcessManager;
  gProcessManager = 0;
 }

void G4ProcessTest::buildTables(const G4String& type, G4bool isPolarised) 
{ 
  G4ParticleDefinition* gamma = G4Gamma::GammaDefinition();
  gamma->SetCuts(1e-3*mm);
  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  electron->SetCuts(1e-3*mm);

  def = createIncidentParticle();

  process = createProcess();
  G4cout << process->GetProcessName() << " created" << G4endl;
 
  brem = createBremsstrahlung();
  if (brem != 0) G4cout << brem->GetProcessName() << " created" << G4endl;

  ioni = createElectronIonisation();
  if (ioni != 0) G4cout << ioni->GetProcessName() << " created" << G4endl;

  eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);

 G4cout << "Now building physics tables, it will take a while..." << G4endl;

  if (def == G4Gamma::GammaDefinition())
    {
      gProcessManager = new G4ProcessManager(def);
      def->SetProcessManager(gProcessManager);
      gProcessManager->AddProcess(process);
      process->BuildPhysicsTable(*def);
    }

  if (def == G4Electron::ElectronDefinition())
    {
      eProcessManager->AddProcess(process);
      process->BuildPhysicsTable(*def);
    }

  // Electron processes are always created; they are needed in photon
  // processes for range tests 
  if (ioni != 0) eProcessManager->AddProcess(ioni);
  if (brem != 0) eProcessManager->AddProcess(brem);
  if (ioni != 0) ioni->BuildPhysicsTable(*electron);       
  if (brem != 0) brem->BuildPhysicsTable(*electron);

  // Build (standard) physics tables for positrons, needed by GammaConversion

  G4ParticleDefinition* positron = G4Positron::PositronDefinition();
  positron->SetCuts(1e-3*mm);
  G4VProcess* theeplusIonisation          = new G4eIonisation();
  G4VProcess* theeplusBremsstrahlung      = new G4eBremsstrahlung();
  G4VProcess* theeplusAnnihilation        = new G4eplusAnnihilation();
  G4ProcessManager* posProcessManager = new G4ProcessManager(positron);
  positron->SetProcessManager(posProcessManager);
  posProcessManager->AddProcess(theeplusIonisation);
  posProcessManager->AddProcess(theeplusBremsstrahlung);
  posProcessManager->AddProcess(theeplusAnnihilation);
  theeplusIonisation->BuildPhysicsTable(*positron);
  theeplusBremsstrahlung->BuildPhysicsTable(*positron);     
  theeplusAnnihilation->BuildPhysicsTable(*positron) ;
}


void G4ProcessTest::postStepTest(const G4Track& track,
				 const G4Step& step) const
{
  G4ParticleChange* particleChange = (G4ParticleChange*) process->PostStepDoIt(track, step);

  G4ProcessTestAnalysis* analysis = G4ProcessTestAnalysis::getInstance();
  analysis->analyseGeneral(track,particleChange); 
  analysis->analyseSecondaries(particleChange); 

  for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
    {
      delete particleChange->GetSecondary(i);
    }
  
  particleChange->Clear();
}

void G4ProcessTest::alongStepTest(const G4Track& track,
				  const G4Step& step) const 
{
  G4ParticleChange* particleChange = (G4ParticleChange*) process->PostStepDoIt(track, step);

  G4ProcessTestAnalysis* analysis = G4ProcessTestAnalysis::getInstance();
  analysis->analyseGeneral(track,particleChange); 
  analysis->analyseSecondaries(particleChange); 

  for (G4int i = 0; i < (particleChange->GetNumberOfSecondaries()); i++) 
    {
      delete particleChange->GetSecondary(i);
    }
  
  particleChange->Clear();
}

G4ParticleDefinition* G4ProcessTest::createIncidentParticle()
{
  return G4Gamma::GammaDefinition();
}
