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
// 22.05.2014 D. Mancusi: Extend INCL++ to 20 GeV
// 19.03.2013 A.Ribon: Replace LEP with FTFP and BERT
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
#include "G4IonConstructor.hh"

#include "G4ComponentGGHadronNucleusXsc.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4LFission.hh"

#include "G4CrossSectionDataSetRegistry.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsINCLXX);

G4ThreadLocal G4HadronPhysicsINCLXX::ThreadPrivate* 
G4HadronPhysicsINCLXX::tpdata = 0;

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(G4int)
    :  G4VPhysicsConstructor("hInelastic INCLXX")
    , QuasiElastic(true)
    , withNeutronHP(false)
    , withFTFP(false)
{
}

G4HadronPhysicsINCLXX::G4HadronPhysicsINCLXX(const G4String& name, const G4bool quasiElastic, const G4bool neutronHP, const G4bool ftfp)
    :  G4VPhysicsConstructor(name) 
    , QuasiElastic(quasiElastic)
    , withNeutronHP(neutronHP)
    , withFTFP(ftfp)
{
}

void G4HadronPhysicsINCLXX::CreateModels()
{
  G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
  G4bool quasiElasticQGS= true;    // For QGS, it must use it.

// initialise fields in tpdata where assignment is optional below.
  tpdata->theNeutronHP=0;
  tpdata->theQGSPNeutron=0;
  tpdata->theQGSPPro=0;
  tpdata->theQGSPPion=0;
  tpdata->theQGSPKaon=0;
  tpdata->theFTFPNeutron=0;
  tpdata->theFTFPPro=0;
  tpdata->theFTFPPion=0;
  tpdata->theFTFPKaon=0;

  tpdata->theNeutrons=new G4NeutronBuilder( withNeutronHP );
  if(withFTFP) {
    tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron=new G4FTFPNeutronBuilder(quasiElasticFTF));
    tpdata->theFTFPNeutron->SetMinEnergy(15.*GeV);
  } else {
    tpdata->theNeutrons->RegisterMe(tpdata->theQGSPNeutron=new G4QGSPNeutronBuilder(quasiElasticQGS));
    tpdata->theQGSPNeutron->SetMinEnergy(15.*GeV);
  }
  tpdata->theNeutrons->RegisterMe(tpdata->theINCLXXNeutron=new G4INCLXXNeutronBuilder);
  tpdata->theINCLXXNeutron->SetMaxEnergy(20.0*GeV);
  if(withNeutronHP) {
    tpdata->theINCLXXNeutron->UsePreCompound(false);
    tpdata->theINCLXXNeutron->SetMinEnergy(19.9*MeV);
    tpdata->theNeutrons->RegisterMe(tpdata->theNeutronHP=new G4NeutronPHPBuilder);
  } else {
    tpdata->theINCLXXNeutron->UsePreCompound(true);
    tpdata->theINCLXXNeutron->SetMinPreCompoundEnergy(0.0*MeV);
    tpdata->theINCLXXNeutron->SetMaxPreCompoundEnergy(2.0*MeV);
    tpdata->theINCLXXNeutron->SetMinEnergy(1.0*MeV);
  }

  tpdata->thePro=new G4ProtonBuilder;
  if(withFTFP) {
    tpdata->thePro->RegisterMe(tpdata->theFTFPPro=new G4FTFPProtonBuilder(quasiElasticFTF));
    tpdata->theFTFPPro->SetMinEnergy(15.*GeV);
  } else {
    tpdata->thePro->RegisterMe(tpdata->theQGSPPro=new G4QGSPProtonBuilder(quasiElasticQGS));
    tpdata->theQGSPPro->SetMinEnergy(15.*GeV);
  }
  tpdata->thePro->RegisterMe(tpdata->theINCLXXPro=new G4INCLXXProtonBuilder);
  tpdata->theINCLXXPro->SetMinEnergy(1.0*MeV);
  tpdata->theINCLXXPro->SetMaxEnergy(20.0*GeV);

  tpdata->thePion=new G4PionBuilder;
  if(withFTFP) {
    tpdata->thePion->RegisterMe(tpdata->theFTFPPion=new G4FTFPPionBuilder(quasiElasticFTF));
    tpdata->theFTFPPion->SetMinEnergy(15.*GeV);
  } else {
    tpdata->thePion->RegisterMe(tpdata->theQGSPPion=new G4QGSPPionBuilder(quasiElasticQGS));
    tpdata->theQGSPPion->SetMinEnergy(15.*GeV);
  }
  tpdata->thePion->RegisterMe(tpdata->theINCLXXPion=new G4INCLXXPionBuilder);
  tpdata->theINCLXXPion->SetMinEnergy(0.0*GeV);
  tpdata->theINCLXXPion->SetMaxEnergy(20.0*GeV);

  tpdata->theKaon=new G4KaonBuilder;
  if(withFTFP) {
    tpdata->theKaon->RegisterMe(tpdata->theFTFPKaon=new G4FTFPKaonBuilder(quasiElasticFTF));
    tpdata->theFTFPKaon->SetMinEnergy(14.*GeV);
  } else {
    tpdata->theKaon->RegisterMe(tpdata->theQGSPKaon=new G4QGSPKaonBuilder(quasiElasticQGS));
    tpdata->theQGSPKaon->SetMinEnergy(14.*GeV);
  }
  tpdata->theKaon->RegisterMe(tpdata->theBertiniKaon=new G4BertiniKaonBuilder);
  tpdata->theBertiniKaon->SetMinEnergy(0.0*GeV);
  tpdata->theBertiniKaon->SetMaxEnergy(15.0*GeV);

  tpdata->theHyperon=new G4HyperonFTFPBuilder;

  tpdata->theAntiBaryon=new G4AntiBarionBuilder;
  tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon=new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

G4HadronPhysicsINCLXX::~G4HadronPhysicsINCLXX()
{
  if(tpdata) {
    delete tpdata->theFTFPNeutron;
    delete tpdata->theQGSPNeutron;
    delete tpdata->theINCLXXNeutron;
    delete tpdata->theNeutronHP;
    delete tpdata->theFTFPPro;
    delete tpdata->theQGSPPro;
    delete tpdata->thePro;
    delete tpdata->theINCLXXPro;
    delete tpdata->theFTFPPion;
    delete tpdata->theQGSPPion;
    delete tpdata->theFTFPKaon;
    delete tpdata->theQGSPKaon;
    delete tpdata->theINCLXXPion;
    delete tpdata->thePion;
    delete tpdata->theKaon;
    delete tpdata->theHyperon;
    delete tpdata->theAntiBaryon;
    delete tpdata->theFTFPAntiBaryon;

    delete tpdata; tpdata = 0;
  }
}

void G4HadronPhysicsINCLXX::ConstructParticle()
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
void G4HadronPhysicsINCLXX::ConstructProcess()
{
  if ( tpdata == 0 ) tpdata = new ThreadPrivate;
  CreateModels();
  tpdata->theNeutrons->Build();
  tpdata->thePro->Build();
  tpdata->thePion->Build();
  tpdata->theKaon->Build();

  // --- Kaons ---
  tpdata->xsKaon = new G4ComponentGGHadronNucleusXsc();
  G4VCrossSectionDataSet * kaonxs = new G4CrossSectionInelastic(tpdata->xsKaon);
  G4PhysListUtil::FindInelasticProcess(G4KaonMinus::KaonMinus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonPlus::KaonPlus())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroShort::KaonZeroShort())->AddDataSet(kaonxs);
  G4PhysListUtil::FindInelasticProcess(G4KaonZeroLong::KaonZeroLong())->AddDataSet(kaonxs);

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
  tpdata->xsNeutronCaptureXS = (G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name());
  capture->AddDataSet(tpdata->xsNeutronCaptureXS);
  G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
  capture->RegisterMe( theNeutronRadCapture );
  if ( withNeutronHP ) {
    capture->AddDataSet( new G4ParticleHPCaptureData );
    theNeutronRadCapture->SetMinEnergy( 19.9*MeV ); 
    if ( ! fission ) {
      fission = new G4HadronFissionProcess("nFission");
      pmanager->AddDiscreteProcess(fission);
    }
    G4LFission* theNeutronLEPFission = new G4LFission();
    theNeutronLEPFission->SetMinEnergy( 19.9*MeV );
    fission->RegisterMe( theNeutronLEPFission );
  }
}
