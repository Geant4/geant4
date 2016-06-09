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
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "HadronPhysicsQGSP_BERT_95.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4PiNuclearCrossSection.hh"
#include "G4HadronInelasticDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(HadronPhysicsQGSP_BERT_95);

HadronPhysicsQGSP_BERT_95::HadronPhysicsQGSP_BERT_95(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSP_BERT_95")
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePiK(0)
    , theLEPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theLEPPro(0)
    , theQGSPPro(0) 
    , theBertiniPro(0)
    , theMiscLHEP(0)
    , QuasiElastic(true)
    , ProjectileDiffraction(false)
    , xsBarashenkovGGPion(0)
    , xsGeisha(0)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
{
}

HadronPhysicsQGSP_BERT_95::HadronPhysicsQGSP_BERT_95(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
    , theNeutrons(0)
    , theLEPNeutron(0)
    , theQGSPNeutron(0)
    , theBertiniNeutron(0)
    , thePiK(0)
    , theLEPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theLEPPro(0)
    , theQGSPPro(0) 
    , theBertiniPro(0)
    , theMiscLHEP(0)
    , QuasiElastic(quasiElastic)
    , ProjectileDiffraction(false)
    , xsBarashenkovGGPion(0)
    , xsGeisha(0)
    , xsAxenWellischGGProton(0)
    , xsLaidlawWellischGGNeutron(0)
{
}

void HadronPhysicsQGSP_BERT_95::CreateModels()
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
  
  thePiK=new G4PiKBuilder;
  thePiK->RegisterMe(theQGSPPiK=new G4QGSPPiKBuilder(QuasiElastic));
  thePiK->RegisterMe(theLEPPiK=new G4LEPPiKBuilder);
  theLEPPiK->SetMaxEnergy(25*GeV);
  theLEPPiK->SetMinEnergy(9.5*GeV);

  thePiK->RegisterMe(theBertiniPiK=new G4BertiniPiKBuilder);
  theBertiniPiK->SetMaxEnergy(9.9*GeV);
  
  theMiscLHEP=new G4MiscLHEPBuilder;
}

HadronPhysicsQGSP_BERT_95::~HadronPhysicsQGSP_BERT_95()
{
   delete theMiscLHEP;
   delete theQGSPNeutron;
   delete theLEPNeutron;
   delete theNeutrons;
   delete theBertiniNeutron;
   delete theQGSPPro;
   delete theLEPPro;
   delete thePro;
   delete theBertiniPro;
   delete theQGSPPiK;
   delete theLEPPiK;
   delete theBertiniPiK;
   delete thePiK;

   delete xsBarashenkovGGPion;
   delete xsGeisha;
   delete xsAxenWellischGGProton;
   delete xsLaidlawWellischGGNeutron;
}

void HadronPhysicsQGSP_BERT_95::ConstructParticle()
{
  G4MesonConstructor pMesonConstructor;
  pMesonConstructor.ConstructParticle();

  G4BaryonConstructor pBaryonConstructor;
  pBaryonConstructor.ConstructParticle();

  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();  
}

#include "G4ProcessManager.hh"
void HadronPhysicsQGSP_BERT_95::ConstructProcess()
{
  CreateModels();
  theNeutrons->Build();
  thePro->Build();
  thePiK->Build();
  theMiscLHEP->Build();

  // Inelastic cross sections

  // --- Pions ---
  // Use Barashenkov inelastic pion cross section up to 91 GeV, 
  // and Glauber-Gribov above
  xsBarashenkovGGPion = new G4CrossSectionPairGG(new G4PiNuclearCrossSection(), 91*GeV);
  G4PhysListUtil::FindInelasticProcess(G4PionPlus::PionPlus())->AddDataSet(xsBarashenkovGGPion);
  G4PhysListUtil::FindInelasticProcess(G4PionMinus::PionMinus())->AddDataSet(xsBarashenkovGGPion);

  // --- Kaons ---
  // Use Geisha inelastic cross sections
  xsGeisha = new G4HadronInelasticDataSet();
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(xsGeisha);

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

  // --- Hyperons ---
  // Use Geisha inelastic cross sections
  G4PhysListUtil::FindInelasticProcess(G4Lambda::Lambda())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiLambda::AntiLambda())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4SigmaMinus::SigmaMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiSigmaMinus::AntiSigmaMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4SigmaPlus::SigmaPlus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiSigmaPlus::AntiSigmaPlus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4XiMinus::XiMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiXiMinus::AntiXiMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4XiZero::XiZero())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiXiZero::AntiXiZero())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4OmegaMinus::OmegaMinus())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiOmegaMinus::AntiOmegaMinus())->AddDataSet(xsGeisha);

  // --- AntiBaryons ---
  // Use Geisha inelastic cross sections
  G4PhysListUtil::FindInelasticProcess(G4AntiProton::AntiProton())->AddDataSet(xsGeisha);
  G4PhysListUtil::FindInelasticProcess(G4AntiNeutron::AntiNeutron())->AddDataSet(xsGeisha);

}
