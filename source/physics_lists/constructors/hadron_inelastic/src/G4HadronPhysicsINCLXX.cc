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
// $Id: G4HadronPhysicsINCLXX.cc 66892 2013-01-17 10:57:59Z gunter $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsINCLXX
//
// Author: 2011 P. Kaitaniemi
//
// Modified:
// 08.03.2013 D. Mancusi: Fix a problem with overlapping model ranges
// 01.03.2013 D. Mancusi: Rename to G4HadronPhysicsINCLXX and introduce
//                        parameters for FTFP and NeutronHP
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 23.03.2012 D. Mancusi: Extended INCL++ to incident heavy ions up to 16O
// 27.11.2011 P.Kaitaniemi: Created physics list for INCL++ using QGSP_INCL_ABLA as a template
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsINCLXX.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsINCLXX);

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(G4int)
    :  G4VPhysicsConstructor("hInelastic INCLXX")
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theHighEnergyNeutron(0)
    , theBertiniNeutron(0)
    , theINCLXXNeutron(0)
    , theNeutronHP(0)
    , thePiK(0)
    , theLEPPiK(0)
    , theHighEnergyPiK(0)
    , theBertiniPiK(0)
    , theINCLXXPiK(0)
    , thePro(0)
    , theLEPPro(0)
    , theHighEnergyPro(0)
    , theBertiniPro(0)
    , theINCLXXPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , theMisc(0)
    , QuasiElastic(true)
    , withNeutronHP(false)
    , withFTFP(false)
{
}

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(const G4String& name, const G4bool quasiElastic, const G4bool neutronHP, const G4bool ftfp)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theHighEnergyNeutron(0)
    , theBertiniNeutron(0)
    , theINCLXXNeutron(0)
    , theNeutronHP(0)
    , thePiK(0)
    , theLEPPiK(0)
    , theHighEnergyPiK(0)
    , theBertiniPiK(0)
    , theINCLXXPiK(0)
    , thePro(0)
    , theLEPPro(0)
    , theHighEnergyPro(0)
    , theBertiniPro(0)
    , theINCLXXPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , theMisc(0)
    , QuasiElastic(quasiElastic)
    , withNeutronHP(neutronHP)
    , withFTFP(ftfp)
{
}

void G4HadronPhysicsINCLXX::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  thePro=new G4ProtonBuilder;
  thePiK=new G4PiKBuilder;
  if(withFTFP) {
    theNeutrons->RegisterMe(theHighEnergyNeutron=new G4FTFPNeutronBuilder(QuasiElastic));
    theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
    theLEPNeutron->SetMinInelasticEnergy(0.*eV);
    theLEPNeutron->SetMaxInelasticEnergy(0.*eV);  

    thePro->RegisterMe(theHighEnergyPro=new G4FTFPProtonBuilder(QuasiElastic));

    thePiK->RegisterMe(theHighEnergyPiK=new G4FTFPPiKBuilder(QuasiElastic));

    theHyperon=new G4HyperonFTFPBuilder;

    theAntiBaryon=new G4AntiBarionBuilder;
    theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
  } else {
    theNeutrons->RegisterMe(theHighEnergyNeutron=new G4QGSPNeutronBuilder(QuasiElastic));
    theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
    theLEPNeutron->SetMinInelasticEnergy(9.5*GeV);
    theLEPNeutron->SetMaxInelasticEnergy(25*GeV);  
    if(withNeutronHP)
      theLEPNeutron->SetMinEnergy(19.9*MeV);

    thePro->RegisterMe(theHighEnergyPro=new G4QGSPProtonBuilder(QuasiElastic));
    thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
    theLEPPro->SetMinEnergy(9.5*GeV);
    theLEPPro->SetMaxEnergy(25*GeV);

    thePiK->RegisterMe(theHighEnergyPiK=new G4QGSPPiKBuilder(QuasiElastic));
    thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
    theLEPPiK->SetMaxEnergy(25*GeV);
    theLEPPiK->SetMinEnergy(9.5*GeV);
  }

  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(2.9*GeV);
  theBertiniNeutron->SetMaxEnergy(9.9*GeV);

  theNeutrons->RegisterMe(theINCLXXNeutron=new G4INCLXXNeutronBuilder);
  theINCLXXNeutron->SetMaxEnergy(3.0*GeV);

  if(withNeutronHP) {
    theINCLXXNeutron->UsePreCompound(false);
    theINCLXXNeutron->SetMinEnergy(19.9*MeV);
    theNeutrons->RegisterMe(theNeutronHP=new G4NeutronHPBuilder);
  } else {
    theINCLXXNeutron->UsePreCompound(true);
    theINCLXXNeutron->SetMinPreCompoundEnergy(0.0*MeV);
    theINCLXXNeutron->SetMaxPreCompoundEnergy(2.0*MeV);
    theINCLXXNeutron->SetMinEnergy(1.0*MeV);
  }

  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMinEnergy(2.9*GeV);
  theBertiniPro->SetMaxEnergy(9.9*GeV);

  thePro->RegisterMe(theINCLXXPro=new G4INCLXXProtonBuilder);
  theINCLXXPro->SetMinEnergy(1.0*MeV);
  theINCLXXPro->SetMaxEnergy(3.0*GeV);
  
  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMinEnergy(2.9*GeV);
  theBertiniPiK->SetMaxEnergy(9.9*GeV);

  thePiK->RegisterMe(theINCLXXPiK=new G4INCLXXPiKBuilder);
  theINCLXXPiK->SetMinEnergy(0.0*GeV);
  theINCLXXPiK->SetMaxEnergy(3.0*GeV);
  
  theMisc = new G4MiscBuilder;
}

G4HadronPhysicsINCLXX::~G4HadronPhysicsINCLXX()
{
   delete theMisc;
   delete theHighEnergyNeutron;
   delete theLEPNeutron;
   delete theBertiniNeutron;
   delete theINCLXXNeutron;
   delete theNeutronHP;
   delete theHighEnergyPro;
   delete theLEPPro;
   delete thePro;
   delete theBertiniPro;
   delete theINCLXXPro;
   delete theHighEnergyPiK;
   delete theLEPPiK;
   delete theINCLXXPiK;
   delete theBertiniPiK;
   delete thePiK;
   delete theHyperon;
   delete theAntiBaryon;
   delete theFTFPAntiBaryon;
}

void G4HadronPhysicsINCLXX::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void G4HadronPhysicsINCLXX::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  if(theHyperon)
    theHyperon->Build();
  if(theAntiBaryon)
    theAntiBaryon->Build();
  theMisc->Build();
}

