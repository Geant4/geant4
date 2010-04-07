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
// $Id: RE01LeptonPhysics.cc,v 1.3 2010-04-07 01:27:53 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "RE01LeptonPhysics.hh"

#include "G4ParticleTable.hh"
#include "G4ProcessManager.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuMultipleScattering.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4TauMinus.hh"
#include "G4TauPlus.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"
#include "G4NeutrinoTau.hh"
#include "G4AntiNeutrinoTau.hh"


RE01LeptonPhysics::RE01LeptonPhysics(const G4String& name)
               :  G4VPhysicsConstructor(name)
{;}


RE01LeptonPhysics::~RE01LeptonPhysics()
{;}


void RE01LeptonPhysics::ConstructParticle()
{ 
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4TauPlus::TauPlusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();
}


void RE01LeptonPhysics::ConstructProcess()
{
  G4ProcessManager * pManager = 0;
  
  // Electron physics

  pManager = G4Electron::Electron()->GetProcessManager();
  pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);  

  //Positron physics

  pManager = G4Positron::Positron()->GetProcessManager(); 
  pManager->AddProcess(new G4eMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4eIonisation(),        -1, 2, 2);
  pManager->AddProcess(new G4eBremsstrahlung(),    -1, 3, 3);  
  pManager->AddProcess(new G4eplusAnnihilation(),   0,-1, 4);

  // Muon-

  pManager = G4MuonMinus::MuonMinus()->GetProcessManager(); 
  pManager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4MuIonisation(),       -1, 2, 2);
  pManager->AddProcess(new G4MuBremsstrahlung(),   -1, 3, 3);  
  pManager->AddProcess(new G4MuPairProduction(),   -1, 4, 4);
  
  // Muon+

  pManager = G4MuonPlus::MuonPlus()->GetProcessManager(); 
  pManager->AddProcess(new G4MuMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4MuIonisation(),       -1, 2, 2);
  pManager->AddProcess(new G4MuBremsstrahlung(),   -1, 3, 3);  
  pManager->AddProcess(new G4MuPairProduction(),   -1, 4, 4);

  // Tau-

  pManager = G4TauMinus::TauMinus()->GetProcessManager();
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);
 
  // Tau+
  
  pManager = G4TauPlus::TauPlus()->GetProcessManager();
  pManager->AddProcess(new G4hMultipleScattering(), -1, 1, 1);
  pManager->AddProcess(new G4hIonisation(),        -1, 2, 2);

}
