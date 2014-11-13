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
// $Id: G4HadronPhysicsQGSP_BIC_AllHP.cc 73040 2013-08-15 09:36:57Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsQGSP_BIC_AllHP
//
// Author: 2006 G.Folger
//
// Based on G4HadronPhysicsQGSP_BIC
//
// Modified:
// 25.04.2007 G.Folger: Add code for quasielastic
// 31.10.2012 A.Ribon: Use G4MiscBuilder
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsQGSP_BIC_AllHP.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4LFission.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsQGSP_BIC_AllHP);

G4ThreadLocal G4HadronPhysicsQGSP_BIC_AllHP::ThreadPrivate*
G4HadronPhysicsQGSP_BIC_AllHP::tpdata = 0;

G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP(G4int)
    :  G4VPhysicsConstructor("hInelastic QGSP_BIC_HP")
/*    , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBinaryNeutron(0)
    , theHPNeutron(0)
    , thePiKB(0)
    , theFTFPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(true)
{}

G4HadronPhysicsQGSP_BIC_AllHP::G4HadronPhysicsQGSP_BIC_AllHP(const G4String& name, G4bool /* quasiElastic */)
    :  G4VPhysicsConstructor(name)  
/*    , theNeutrons(0)
    , theFTFPNeutron(0)
    , theQGSPNeutron(0)
    , theBinaryNeutron(0)
    , theHPNeutron(0)
    , thePiKB(0)
    , theFTFPPiK(0)
    , theQGSPPiK(0)
    , theBertiniPiK(0)
    , thePro(0)
    , theFTFPPro(0)
    , theQGSPPro(0)
    , theBinaryPro(0)
    , theHyperon(0)
    , theAntiBaryon(0)
    , theFTFPAntiBaryon(0)
    , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(quasiElastic)
{}

void G4HadronPhysicsQGSP_BIC_AllHP::CreateModels()
{
  G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool quasiElasticQGS= true;    // For QGS, it must use it.

  const G4double maxFTFP = 25.0*GeV;
  const G4double minFTFP =  9.5*GeV;
  const G4double maxBIC  =  9.9*GeV;
  const G4double maxBERT =  5.0*GeV;
  const G4double maxHP   = 19.9*MeV;

  tpdata->theNeutronB=new G4NeutronBuilder( true ); // Fission on
  tpdata->theNeutronB->RegisterMe(tpdata->theQGSPNeutron=new G4QGSPNeutronBuilder(quasiElasticQGS));
  tpdata->theNeutronB->RegisterMe(tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(quasiElasticFTF));
  tpdata->theFTFPNeutron->SetMinEnergy(minFTFP);
  tpdata->theFTFPNeutron->SetMaxEnergy(maxFTFP);

  tpdata->theNeutronB->RegisterMe(tpdata->theBinaryNeutron=new G4BinaryNeutronBuilder);
  tpdata->theBinaryNeutron->SetMinEnergy(maxHP);
  tpdata->theBinaryNeutron->SetMaxEnergy(maxBIC);

  //ParticleHP
  tpdata->theNeutronB->RegisterMe(tpdata->thePHPNeutron=new G4NeutronPHPBuilder);

  tpdata->theProtonB=new G4ProtonBuilder;
  tpdata->theProtonB->RegisterMe(tpdata->theQGSPProton=new G4QGSPProtonBuilder(quasiElasticQGS));
  tpdata->theProtonB->RegisterMe(tpdata->theFTFPProton=new G4FTFPProtonBuilder(quasiElasticFTF));
  tpdata->theFTFPProton->SetMinEnergy(minFTFP);
  tpdata->theFTFPProton->SetMaxEnergy(maxFTFP);
  tpdata->theProtonB->RegisterMe(tpdata->theBinaryProton=new G4BinaryProtonBuilder);
  tpdata->theBinaryProton->SetMaxEnergy(maxBIC);

  //ParticleHP
  tpdata->theBinaryProton->SetMinEnergy(200*MeV);
  tpdata->thePHPProton=new G4ProtonPHPBuilder;
  tpdata->theProtonB->RegisterMe(tpdata->thePHPProton);
  tpdata->thePHPProton->SetMinEnergy(0.*MeV); 
  tpdata->thePHPProton->SetMaxEnergy(200*MeV); 

  //ParticleHP
  tpdata->theDeuteronB=new G4DeuteronBuilder;
  tpdata->thePHPDeuteron=new G4DeuteronPHPBuilder;
  tpdata->theDeuteronB->RegisterMe(tpdata->thePHPDeuteron);
  tpdata->thePHPDeuteron->SetMaxEnergy(200*MeV);
  tpdata->theDeuteronB->RegisterMe(tpdata->theBinaryDeuteron=new G4BinaryDeuteronBuilder);
  tpdata->theBinaryDeuteron->SetMaxEnergy(maxBIC);


  //ParticleHP
  tpdata->theTritonB=new G4TritonBuilder;
  tpdata->thePHPTriton=new G4TritonPHPBuilder;
  tpdata->theTritonB->RegisterMe(tpdata->thePHPTriton);
  tpdata->thePHPTriton->SetMaxEnergy(200*MeV);
  tpdata->theTritonB->RegisterMe(tpdata->theBinaryTriton=new G4BinaryTritonBuilder);
  tpdata->theBinaryTriton->SetMaxEnergy(maxBIC);

  //ParticleHP
  tpdata->theHe3B=new G4He3Builder;
  tpdata->thePHPHe3=new G4He3PHPBuilder;
  tpdata->theHe3B->RegisterMe(tpdata->thePHPHe3);
  tpdata->thePHPHe3->SetMaxEnergy(200*MeV);
  tpdata->theHe3B->RegisterMe(tpdata->theBinaryHe3=new G4BinaryHe3Builder);
  tpdata->theBinaryHe3->SetMaxEnergy(maxBIC);

  //ParticleHP
  tpdata->theAlphaB=new G4AlphaBuilder;
  tpdata->thePHPAlpha=new G4AlphaPHPBuilder;
  tpdata->theAlphaB->RegisterMe(tpdata->thePHPAlpha);
  tpdata->thePHPAlpha->SetMaxEnergy(200*MeV);
  tpdata->theAlphaB->RegisterMe(tpdata->theBinaryAlpha=new G4BinaryAlphaBuilder);
  tpdata->theBinaryAlpha->SetMaxEnergy(maxBIC);

  
  tpdata->thePiKB=new G4PiKBuilder;
  tpdata->thePiKB->RegisterMe(tpdata->theQGSPPiK=new G4QGSPPiKBuilder(quasiElasticQGS));
  tpdata->thePiKB->RegisterMe(tpdata->theFTFPPiK=new G4FTFPPiKBuilder(quasiElasticFTF));
  tpdata->theFTFPPiK->SetMaxEnergy(maxFTFP);
  tpdata->thePiKB->RegisterMe(tpdata->theBertiniPiK=new G4BertiniPiKBuilder);
  tpdata->theBertiniPiK->SetMaxEnergy(maxBERT);

  tpdata->theHyperon=new G4HyperonFTFPBuilder;

  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

G4HadronPhysicsQGSP_BIC_AllHP::~G4HadronPhysicsQGSP_BIC_AllHP() 
{
  //ParticleHP
   delete tpdata->thePHPNeutron;
   delete tpdata->thePHPProton;
   delete tpdata->thePHPDeuteron;
   delete tpdata->thePHPTriton;
   delete tpdata->thePHPHe3;
   delete tpdata->thePHPAlpha;
   delete tpdata->theBinaryNeutron;
   delete tpdata->theQGSPNeutron;
   delete tpdata->theFTFPNeutron;
   delete tpdata->theBertiniPiK;
   delete tpdata->theQGSPPiK;
   delete tpdata->theFTFPPiK;
   delete tpdata->thePiKB;
   delete tpdata->theBinaryProton;
   delete tpdata->theQGSPProton;
   delete tpdata->theFTFPProton;
   delete tpdata->theProtonB;
   delete tpdata->theFTFPAntiBaryon;
   delete tpdata->theAntiBaryon;
   delete tpdata->theHyperon;
   delete tpdata->xsNeutronCaptureXS;

   delete tpdata; tpdata = 0;
}

void G4HadronPhysicsQGSP_BIC_AllHP::ConstructParticle()
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
void G4HadronPhysicsQGSP_BIC_AllHP::ConstructProcess()
{
  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  tpdata->theNeutronB->Build();
  tpdata->theProtonB->Build();
  tpdata->thePiKB->Build();
  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  G4HadronicProcess* capture = 0;
  G4HadronicProcess* fission = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    } else if ( fFission == ((*pv)[i])->GetProcessSubType() ) {
      fission = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  tpdata->xsNeutronCaptureXS = new G4NeutronCaptureXS();
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->AddDataSet( new G4NeutronHPCaptureData );
  G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
  theNeutronRadCapture->SetMinEnergy( 19.9*MeV ); 
  capture->RegisterMe( theNeutronRadCapture );
  if ( ! fission ) {
    fission = new G4HadronFissionProcess("nFission");
    pmanager->AddDiscreteProcess(fission);
  }
  G4LFission* theNeutronLEPFission = new G4LFission();
  theNeutronLEPFission->SetMinEnergy( 19.9*MeV );
  fission->RegisterMe( theNeutronLEPFission );
}

