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
// ClassName:   G4HadronPhysicsQGS_BIC
//
// Author: 2007 Gunter Folger
//     created from G4HadronPhysicsQGSP_BIC  by H.P.Wellisch
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGS_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGS_BIC);

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(G4int)
    :  G4VPhysicsConstructor("hInelastic QGS_BIC")
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theFTFBinaryNeutron(0)
    , theQGSBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theBinaryPion(0)
    , theBertiniPion(0)
    , theFTFBinaryPion(0)
    , theQGSBinaryPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFBinaryKaon(0)
    , theQGSBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theQGSBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(true)
{}

G4HadronPhysicsQGS_BIC::G4HadronPhysicsQGS_BIC(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theFTFBinaryNeutron(0)
    , theQGSBinaryNeutron(0)
    , theBinaryNeutron(0)
    , thePion(0)
    , theBinaryPion(0)
    , theBertiniPion(0)
    , theFTFBinaryPion(0)
    , theQGSBinaryPion(0)
    , theKaon(0)
    , theBertiniKaon(0)
    , theFTFBinaryKaon(0)
    , theQGSBinaryKaon(0)
    , thePro(0)
    , theFTFBinaryPro(0)
    , theQGSBinaryPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , QuasiElastic(quasiElastic)
{}

void G4HadronPhysicsQGS_BIC::CreateModels()
{
  G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool quasiElasticQGS= true;    // For QGS, it must use it.

  const G4double maxFTFP =     25.0*GeV;
  const G4double minFTFP =      9.5*GeV;
  const G4double maxBIC  =      9.9*GeV;
  const G4double maxPionBIC =   1.3*GeV;
  const G4double maxPionBERT =  5.0*GeV;
  const G4double minPionBERT =  1.2*GeV;
  const G4double maxKaonBERT =  5.0*GeV;

  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSBinaryNeutron=new G4QGSBinaryNeutronBuilder(quasiElasticQGS));
  theNeutrons->RegisterMe(theFTFBinaryNeutron=new G4FTFBinaryNeutronBuilder(quasiElasticFTF));
  theFTFBinaryNeutron->SetMinEnergy(minFTFP);
  theFTFBinaryNeutron->SetMaxEnergy(maxFTFP);
  // Exclude LEP only from Inelastic 
  //  -- Register it for other processes: Capture, Elastic
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(0.0*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(0.0*GeV);

  theNeutrons->RegisterMe(theBinaryNeutron=new G4BinaryNeutronBuilder);
  theBinaryNeutron->SetMaxEnergy(maxBIC);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSBinaryPro=new G4QGSBinaryProtonBuilder(quasiElasticQGS));
  thePro->RegisterMe(theFTFBinaryPro=new G4FTFBinaryProtonBuilder(quasiElasticFTF));
  theFTFBinaryPro->SetMinEnergy(minFTFP);
  theFTFBinaryPro->SetMaxEnergy(maxFTFP);

  thePro->RegisterMe(theBinaryPro=new G4BinaryProtonBuilder);
  theBinaryPro->SetMaxEnergy(maxBIC);

  thePion=new G4PionBuilder;
  thePion->RegisterMe(theQGSBinaryPion=new G4QGSBinaryPionBuilder(quasiElasticQGS));
  thePion->RegisterMe(theFTFBinaryPion=new G4FTFBinaryPionBuilder(quasiElasticFTF));
  theFTFBinaryPion->SetMaxEnergy(maxFTFP);
  thePion->RegisterMe(theBertiniPion=new G4BertiniPionBuilder);
  theBertiniPion->SetMinEnergy(minPionBERT);
  theBertiniPion->SetMaxEnergy(maxPionBERT);
  thePion->RegisterMe(theBinaryPion = new G4BinaryPionBuilder);
  theBinaryPion->SetMaxEnergy(maxPionBIC);

  theKaon=new G4KaonBuilder;
  theKaon->RegisterMe(theQGSBinaryKaon=new G4QGSBinaryKaonBuilder(quasiElasticQGS));
  theKaon->RegisterMe(theFTFBinaryKaon=new G4FTFBinaryKaonBuilder(quasiElasticFTF));
  theFTFBinaryKaon->SetMaxEnergy(maxFTFP);
  theKaon->RegisterMe(theBertiniKaon=new G4BertiniKaonBuilder);
  theBertiniKaon->SetMaxEnergy(maxKaonBERT);

  theHyperon=new G4HyperonFTFPBuilder;

  theAntiBaryon=new G4AntiBarionBuilder;
  theAntiBaryon->RegisterMe(theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

G4HadronPhysicsQGS_BIC::~G4HadronPhysicsQGS_BIC() 
{
   delete theBinaryNeutron;
   delete theQGSBinaryNeutron;
   delete theFTFBinaryNeutron;
   delete theLEPNeutron;
   delete theNeutrons;
   delete theQGSBinaryPion;
   delete theFTFBinaryPion;
   delete theBertiniPion;
   delete theBinaryPion;
   delete thePion;
   delete theQGSBinaryKaon;
   delete theFTFBinaryKaon;
   delete theBertiniKaon;
   delete theKaon;
   delete theBinaryPro;
   delete theQGSBinaryPro;
   delete theFTFBinaryPro;
   delete thePro;
   delete theFTFPAntiBaryon;
   delete theAntiBaryon;
   delete theHyperon;
}

void G4HadronPhysicsQGS_BIC::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();

  G4IonConstructor pIonConstructor;
  pIonConstructor.ConstructParticle();
}

#include "G4ProcessManager.hh"
void G4HadronPhysicsQGS_BIC::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePion->Build();
  theKaon->Build();
  theHyperon->Build();
  theAntiBaryon->Build();
}

