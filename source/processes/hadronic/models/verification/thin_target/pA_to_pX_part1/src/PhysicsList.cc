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
// $Id: PhysicsList.cc,v 1.2 2003-05-29 16:56:33 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#include "PhysicsList.hh"
#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"

#include "G4LEProtonInelastic.hh"
#include "G4BinaryCascade.hh"
#include "G4CascadeInterface.hh"


PhysicsList::PhysicsList()
{;}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::ConstructParticle()
{
  G4Proton::ProtonDefinition();

  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();
  G4TauPlus::TauPlusDefinition();
  G4TauMinus::TauMinusDefinition();
  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
  G4NeutrinoTau::NeutrinoTauDefinition();
  G4AntiNeutrinoTau::AntiNeutrinoTauDefinition();

  // mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();

}

void PhysicsList::ConstructProcess()
{
  // Define transportation process

  AddTransportation();
  ConstructHad();
}

void PhysicsList::ConstructHad()
{
  G4HadronicInteraction * theModel(0);
  if(getenv("TestBinary"))
  {
    G4BinaryCascade* theBinaryCascade = new G4BinaryCascade();
    theModel = theBinaryCascade;
  }
  //  theBinaryCascade->SetMinEnergy(100.0*MeV);
  else if(getenv("TestBertini"))
  {
    G4CascadeInterface* theBertiniCascade = new G4CascadeInterface();
    theModel = theBertiniCascade;
  }
  //  theBertiniCascade->SetMinEnergy(100.0*MeV);
  else  if(getenv("TestLEP"))
  {
    G4LEProtonInelastic* theLEProtonModel = new G4LEProtonInelastic();
    theModel = theLEProtonModel;
  }
  else
  {
    G4Exception("PhysicsList::ConstructHad(): ERROR - No model selected");
  }
  //  theLEProtonModel->SetMaxEnergy(100.0*MeV);


  // proton
  //
  G4ProcessManager* pManager = G4Proton::Proton()->GetProcessManager();
  //  theProtonInelastic.RegisterMe(theBinaryCascade);
  //  theProtonInelastic.RegisterMe(theBertiniCascade);
  theProtonInelastic.RegisterMe(theModel);
  pManager->AddDiscreteProcess(&theProtonInelastic);

}

void PhysicsList::SetCuts()
{
  // suppress error messages even in case e/gamma/proton do not exist    
  G4int temp = GetVerboseLevel();
  SetVerboseLevel(0);                                    
  SetCutsWithDefault();   

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}




