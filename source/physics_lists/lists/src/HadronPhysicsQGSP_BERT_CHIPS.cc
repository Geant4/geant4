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
// ClassName:   HadronPhysicsQGSP_BERT_CHIPS
//
// Author: 2010 G.Folger
//   derived from HadronPhysicsQGSP_BERT
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsQGSP_BERT_CHIPS.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4HadronicProcessType.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4QHadronInelasticDataSet.hh"

#include "G4HadronInelasticDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4PhysListUtil.hh"


HadronPhysicsQGSP_BERT_CHIPS::HadronPhysicsQGSP_BERT_CHIPS(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSP_BERT_CHIPS")
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePion(0)
    , theLEPPion(0)
    , theQGSPPion(0)
    , theBertiniPion(0)
    , theKaon(0)
    , thePro(0)
    , theLEPPro(0)
    , theQGSPPro(0)
    , theBertiniPro(0)
    , theMiscCHIPS(0)
    , QuasiElastic(true)
    , ProjectileDiffraction(false)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
{
}

HadronPhysicsQGSP_BERT_CHIPS::HadronPhysicsQGSP_BERT_CHIPS(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name) 
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePion(0)
    , theLEPPion(0)
    , theQGSPPion(0)
    , theBertiniPion(0)
    , theKaon(0)
    , thePro(0)
    , theLEPPro(0)
    , theQGSPPro(0)
    , theBertiniPro(0)
    , theMiscCHIPS(0)
    , QuasiElastic(quasiElastic)
    , ProjectileDiffraction(false)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
{
}

void HadronPhysicsQGSP_BERT_CHIPS::CreateModels()
{
  theNeutrons=new G4NeutronBuilder;
  theNeutrons->RegisterMe(theQGSPNeutron=new G4QGSPNeutronBuilder(QuasiElastic, ProjectileDiffraction));
  theNeutrons->RegisterMe(theLEPNeutron=new G4LEPNeutronBuilder);
  theLEPNeutron->SetMinInelasticEnergy(9.5*GeV);
  theLEPNeutron->SetMaxInelasticEnergy(25*GeV);  

  theNeutrons->RegisterMe(theBertiniNeutron=new G4BertiniNeutronBuilder);
  theBertiniNeutron->SetMinEnergy(0.0*GeV);
  theBertiniNeutron->SetMaxEnergy(9.9*GeV);

  thePro=new G4ProtonBuilder;
  thePro->RegisterMe(theQGSPPro=new G4QGSPProtonBuilder(QuasiElastic, ProjectileDiffraction));
  thePro->RegisterMe(theLEPPro=new G4LEPProtonBuilder);
  theLEPPro->SetMinEnergy(9.5*GeV);
  theLEPPro->SetMaxEnergy(25*GeV);

  thePro->RegisterMe(theBertiniPro=new G4BertiniProtonBuilder);
  theBertiniPro->SetMaxEnergy(9.9*GeV);
  
  thePion=new G4PionBuilder;
  thePion->RegisterMe(theQGSPPion=new G4QGSPPionBuilder(QuasiElastic));
  thePion->RegisterMe(theLEPPion=new G4LEPPionBuilder);
  theLEPPion->SetMaxEnergy(25*GeV);
  theLEPPion->SetMinEnergy(9.5*GeV);

  thePion->RegisterMe(theBertiniPion=new G4BertiniPionBuilder);
  theBertiniPion->SetMaxEnergy(9.9*GeV);

  G4int verbosity(0);
  theKaon=new G4ChipsKaonBuilder(verbosity);  // is self contained, use G4QInelastic   
  
  theMiscCHIPS=new G4MiscCHIPSBuilder;
}

HadronPhysicsQGSP_BERT_CHIPS::~HadronPhysicsQGSP_BERT_CHIPS()
{
   delete theMiscCHIPS;
   delete theQGSPNeutron;
   delete theLEPNeutron;
   delete theBertiniNeutron;
   delete theQGSPPro;
   delete theLEPPro;
   delete thePro;
   delete theBertiniPro;
   delete theQGSPPion;
   delete theLEPPion;
   delete theBertiniPion;
   delete thePion;
   delete theKaon;

   delete xsAxenWellischGGProton;
   delete xsLaidlawWellischGGNeutron;
}

void HadronPhysicsQGSP_BERT_CHIPS::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_BERT_CHIPS::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePion->Build();
  theKaon->Build();   // has CHIPS cross sections for Kaons
  theMiscCHIPS->Build();

  // Inelastic cross sections

  // --- Protons ---
  // Use Axen-Wellisch inelastic proton cross section up to 91 GeV,
  // and Glauber-Gribov above
  xsAxenWellischGGProton = new G4CrossSectionPairGG(new G4ProtonInelasticCrossSection(), 91*GeV); 
  G4PhysListUtil::FindInelasticProcess(G4Proton::Proton())->AddDataSet(xsAxenWellischGGProton);

  // --- Neutrons ---
  // Use Laidlaw-Wellisch inelastic neutron cross section up to 91 GeV,
  // and Glauber-Gribov above
  xsLaidlawWellischGGNeutron = new G4CrossSectionPairGG(new G4NeutronInelasticCrossSection(), 91*GeV);
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(xsLaidlawWellischGGNeutron);
}
