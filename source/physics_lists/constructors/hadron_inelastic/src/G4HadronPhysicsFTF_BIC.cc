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
// $Id: G4HadronPhysicsFTF_BIC.cc 83699 2014-09-10 07:18:25Z gcosmo $
//
//---------------------------------------------------------------------------
//
// ClassName:   G4HadronPhysicsFTF_BIC
//
// Author: 2007 Gunter Folger
//
// Modified:
//
//----------------------------------------------------------------------------
//
#include <iomanip>   

#include "G4HadronPhysicsFTF_BIC.hh"

#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4CrossSectionDataSetRegistry.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsFTF_BIC);

G4ThreadLocal G4HadronPhysicsFTF_BIC::ThreadPrivate* 
G4HadronPhysicsFTF_BIC::tpdata = 0;

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(G4int)
    :  G4VPhysicsConstructor("hInelastic FTF_BIC")
/*    , theNeutrons(0)
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
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
    , QuasiElastic(false)
{}

G4HadronPhysicsFTF_BIC::G4HadronPhysicsFTF_BIC(const G4String& name, G4bool quasiElastic)
    :  G4VPhysicsConstructor(name)
/*    , theNeutrons(0)
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
    , xsNeutronInelasticXS(0)
    , xsNeutronCaptureXS(0)*/
    , QuasiElastic(quasiElastic)
{}

void G4HadronPhysicsFTF_BIC::CreateModels()
{
  tpdata->theNeutrons=new G4NeutronBuilder;

  tpdata->theNeutrons->RegisterMe(tpdata->theFTFBinaryNeutron=new G4FTFBinaryNeutronBuilder(QuasiElastic));

  tpdata->theNeutrons->RegisterMe(tpdata->theBinaryNeutron=new G4BinaryNeutronBuilder);
  tpdata->theBinaryNeutron->SetMaxEnergy(5.0*GeV);

  tpdata->thePro=new G4ProtonBuilder;
  tpdata->thePro->RegisterMe(tpdata->theFTFBinaryPro=new G4FTFBinaryProtonBuilder(QuasiElastic));

  tpdata->thePro->RegisterMe(tpdata->theBinaryPro=new G4BinaryProtonBuilder);
  tpdata->theBinaryPro->SetMaxEnergy(5.0*GeV);
  
  tpdata->thePion=new G4PionBuilder;
  tpdata->thePion->RegisterMe(tpdata->theFTFBinaryPion=new G4FTFBinaryPionBuilder(QuasiElastic));
  tpdata->thePion->RegisterMe(tpdata->theBICPion = new G4BinaryPionBuilder);
  tpdata->theBICPion->SetMaxEnergy(5*GeV);        //  use Binary up to 5GeV for pion

  tpdata->theKaon=new G4KaonBuilder;
  tpdata->theKaon->RegisterMe(tpdata->theFTFBinaryKaon=new G4FTFBinaryKaonBuilder(QuasiElastic));
  tpdata->theKaon->RegisterMe(tpdata->theBertiniKaon=new G4BertiniKaonBuilder);
  tpdata->theBertiniKaon->SetMaxEnergy(5*GeV);
 
  tpdata->theHyperon=new G4HyperonFTFPBuilder;
    
  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new  G4FTFPAntiBarionBuilder(QuasiElastic));
}

G4HadronPhysicsFTF_BIC::~G4HadronPhysicsFTF_BIC() 
{
  if (!tpdata) return;

   delete tpdata->theFTFBinaryNeutron;
   delete tpdata->theBinaryNeutron;
   delete tpdata->theNeutrons;

   delete tpdata->theFTFBinaryPro;
   delete tpdata->theBinaryPro;
   delete tpdata->thePro;

   delete tpdata->theFTFBinaryPion;
   delete tpdata->theBICPion;
   delete tpdata->thePion;

   delete tpdata->theFTFBinaryKaon;
   delete tpdata->theBertiniKaon;
   delete tpdata->theKaon;

   delete tpdata->theHyperon;
   delete tpdata->theAntiBaryon;
   delete tpdata->theFTFPAntiBaryon;

   delete tpdata; tpdata = 0;
}

void G4HadronPhysicsFTF_BIC::ConstructParticle()
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
void G4HadronPhysicsFTF_BIC::ConstructProcess()
{
  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePion->Build();
  tpdata->theKaon->Build();
  
  tpdata->theHyperon->Build();
  tpdata->theAntiBaryon->Build();

  // --- Neutrons ---
  tpdata->xsNeutronInelasticXS = (G4NeutronInelasticXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronInelasticXS::Default_Name());
  G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron())->AddDataSet(tpdata->xsNeutronInelasticXS);

  G4HadronicProcess* capture = 0;
  G4ProcessManager* pmanager = G4Neutron::Neutron()->GetProcessManager();
  G4ProcessVector*  pv = pmanager->GetProcessList();
  for ( size_t i=0; i < static_cast<size_t>(pv->size()); ++i ) {
    if ( fCapture == ((*pv)[i])->GetProcessSubType() ) {
      capture = static_cast<G4HadronicProcess*>((*pv)[i]);
    }
  }
  if ( ! capture ) {
    capture = new G4HadronCaptureProcess("nCapture");
    pmanager->AddDiscreteProcess(capture);
  }
  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  capture->RegisterMe(new G4NeutronRadCapture());
}

