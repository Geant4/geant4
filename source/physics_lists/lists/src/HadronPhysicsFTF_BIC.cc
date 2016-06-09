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
// $Id$
//
//---------------------------------------------------------------------------
//
// ClassName:   HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsFTF_BIC.hh"

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
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsFTF_BIC);

HadronPhysicsFTF_BIC::HadronPhysicsFTF_BIC(G4int)
    :  G4VPhysicsConstructor("hInelastic FTF_BIC")
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theFTFBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theKaon(0)
    , theBICPion(0)
    , theBertiniKaon(0)
    , theFTFBinaryPion(0)
    , theFTFBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(false)
{}

HadronPhysicsFTF_BIC::HadronPhysicsFTF_BIC(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theFTFBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theKaon(0)
    , theBICPion(0)
    , theBertiniKaon(0)
    , theFTFBinaryPion(0)
    , theFTFBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(quasiElastic)
{}

void HadronPhysicsFTF_BIC::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;

  theNeutrons->RegisterMe(theFTFBinaryNeutron=new G4FTFBinaryNeutronBuilder(QuasiElastic));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.*eV);
  theLEPNeutron->SetMaxInelasticEnergy(0.*eV);  

  theNeutrons->RegisterMe(theBinaryNeutron=new G4BinaryNeutronBuilder);
  theBinaryNeutron->SetMaxEnergy(5.0*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theFTFBinaryPro=new G4FTFBinaryProtonBuilder(QuasiElastic));

  thePro->RegisterMe(theBinaryPro=new G4BinaryProtonBuilder);
  theBinaryPro->SetMaxEnergy(5.0*GeV);
  
  thePion=new G4PionBuilder;
  thePion->RegisterMe(theFTFBinaryPion=new G4FTFBinaryPionBuilder(QuasiElastic));
  thePion->RegisterMe(theBICPion = new G4BinaryPionBuilder);
  theBICPion->SetMaxEnergy(5*GeV);        //  use Binary up to 5GeV for pion

  theKaon=new G4KaonBuilder;
  theKaon->RegisterMe(theFTFBinaryKaon=new G4FTFBinaryKaonBuilder(QuasiElastic));
  theKaon->RegisterMe(theBertiniKaon=new G4BertiniKaonBuilder);
  theBertiniKaon->SetMaxEnergy(5*GeV);

  
  theHyperon=new G4HyperonFTFPBuilder;
    
  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
}

HadronPhysicsFTF_BIC::~HadronPhysicsFTF_BIC() 
{
   delete theFTFBinaryNeutron;
   delete theLEPNeutron;
   delete theBinaryNeutron;
   delete theNeutrons;

   delete theFTFBinaryPro;
   delete theBinaryPro;
   delete thePro;

   delete theFTFBinaryPion;
   delete theBICPion;
   delete thePion;

   delete theFTFBinaryKaon;
   delete theBertiniKaon;
   delete theKaon;

   delete theHyperon;
   delete theAntiBaryon;
   delete theFTFPAntiBaryon;
  
}

void HadronPhysicsFTF_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

//#include "G4ProcessManager.hh"
#include "G4PhysListUtil.hh"
void HadronPhysicsFTF_BIC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePion->Build();
  theKaon->Build();
  
  theHyperon->Build();
  theAntiBaryon->Build();
}

