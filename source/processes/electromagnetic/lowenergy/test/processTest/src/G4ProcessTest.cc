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
// $Id: G4ProcessTest.cc,v 1.8 2001-11-06 01:38:24 pia Exp $
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
// Test DoIt method of physics processes
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
  process = createProcess();
  G4cout << process->GetProcessName() << " created" << G4endl;

  def = createIncidentParticle();

  brem = createBremsstrahlung();
  if (brem != 0)  
    {
      G4cout << brem->GetProcessName() << " created" << G4endl;
    }

  ioni = createElectronIonisation();
  if (ioni != 0)  
    {
      G4cout << ioni->GetProcessName() << " created" << G4endl;
    }

  if (def == G4Gamma::GammaDefinition())
    {
      def->SetCuts(1e-3*mm);
      gProcessManager = new G4ProcessManager(def);
      def->SetProcessManager(gProcessManager);
      gProcessManager->AddProcess(process);
      process->BuildPhysicsTable(*def);
    }

  G4ParticleDefinition* electron = G4Electron::ElectronDefinition();
  electron->SetCuts(1e-3*mm);
  eProcessManager = new G4ProcessManager(electron);
  electron->SetProcessManager(eProcessManager);

  G4cout << "Now building physics tables, it will take a while..." << G4endl;

  if (def == G4Electron::ElectronDefinition())
    {
      eProcessManager->AddProcess(process);
      process->BuildPhysicsTable(*def); 
    }
  if (ioni != 0)  
    {
      eProcessManager->AddProcess(ioni);
      ioni->BuildPhysicsTable(*electron);
    }
  if (brem != 0)  
    {
      eProcessManager->AddProcess(ioni);
      brem->BuildPhysicsTable(*electron);
    }
}
 
void G4ProcessTest::postStepTest(const G4Track& track,
				 const G4Step& step) const
{
  G4ParticleChange* particleChange = (G4ParticleChange*) process->PostStepDoIt(track, step);

  G4ProcessTestAnalysis* analyser = G4ProcessTestAnalysis::getInstance();
  analyser->analyseGeneral(track,particleChange); 
  analyser->analyseSecondaries(particleChange); 

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

  G4ProcessTestAnalysis* analyser = G4ProcessTestAnalysis::getInstance();
  analyser->analyseGeneral(track,particleChange); 
  analyser->analyseSecondaries(particleChange); 

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
