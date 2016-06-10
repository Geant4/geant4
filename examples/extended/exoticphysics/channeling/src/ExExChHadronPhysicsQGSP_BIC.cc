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

#include <iomanip>

#include "ExExChHadronPhysicsQGSP_BIC.hh"

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
#include "G4NeutronInelasticXS.hh"
#include "G4NeutronCaptureXS.hh"

#include "G4PhysListUtil.hh"

// factory
#include "G4PhysicsConstructorFactory.hh"
//

G4_DECLARE_PHYSCONSTR_FACTORY(ExExChHadronPhysicsQGSP_BIC);

G4ThreadLocal ExExChHadronPhysicsQGSP_BIC::ThreadPrivate*
ExExChHadronPhysicsQGSP_BIC::tpdata = 0;

ExExChHadronPhysicsQGSP_BIC::ExExChHadronPhysicsQGSP_BIC(G4int)
:  G4VPhysicsConstructor("hInelastic Xtal QGSP_BIC")
/*    , theNeutrons(0)
 , theFTFPNeutron(0)
 , theQGSPNeutron(0)
 , theBinaryNeutron(0)
 , thePiK(0)
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
 , xsNeutronInelasticXS(0)
 , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(true)
{}

ExExChHadronPhysicsQGSP_BIC::ExExChHadronPhysicsQGSP_BIC(
                            const G4String& name, G4bool /* quasiElastic */)
:  G4VPhysicsConstructor(name)
/*    , theNeutrons(0)
 , theFTFPNeutron(0)
 , theQGSPNeutron(0)
 , theBinaryNeutron(0)
 , thePiK(0)
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
 , xsNeutronInelasticXS(0)
 , xsNeutronCaptureXS(0)*/
//    , QuasiElastic(quasiElastic)
{}

void ExExChHadronPhysicsQGSP_BIC::CreateModels()
{
    G4bool quasiElasticFTF= false;   // Use built-in quasi-elastic (not add-on)
    G4bool quasiElasticQGS= true;    // For QGS, it must use it.
    
    const G4double maxFTFP = 25.0*GeV;
    const G4double minFTFP =  9.5*GeV;
    const G4double maxBIC  =  9.9*GeV;
    const G4double maxBERT =  5.0*GeV;
    
    tpdata->theNeutrons=new G4NeutronBuilder;
    tpdata->theNeutrons->RegisterMe(tpdata->theQGSPNeutron =
                                    new G4QGSPNeutronBuilder(quasiElasticQGS));
    tpdata->theNeutrons->RegisterMe(tpdata->theFTFPNeutron =
                                    new G4FTFPNeutronBuilder(quasiElasticFTF));
    tpdata->theFTFPNeutron->SetMinEnergy(minFTFP);
    tpdata->theFTFPNeutron->SetMaxEnergy(maxFTFP);
    
    tpdata->theNeutrons->RegisterMe(tpdata->theBinaryNeutron =
                                    new G4BinaryNeutronBuilder);
    tpdata->theBinaryNeutron->SetMaxEnergy(maxBIC);
    
    tpdata->thePro=new ExExChProtonBuilder;
    tpdata->thePro->RegisterMe(tpdata->theQGSPPro =
                               new G4QGSPProtonBuilder(quasiElasticQGS));
    tpdata->thePro->RegisterMe(tpdata->theFTFPPro =
                               new G4FTFPProtonBuilder(quasiElasticFTF));
    tpdata->theFTFPPro->SetMinEnergy(minFTFP);
    tpdata->theFTFPPro->SetMaxEnergy(maxFTFP);
    
    tpdata->thePro->RegisterMe(tpdata->theBinaryPro=new G4BinaryProtonBuilder);
    tpdata->theBinaryPro->SetMaxEnergy(maxBIC);
    
    tpdata->thePiK=new ExExChPiKBuilder;
    tpdata->thePiK->RegisterMe(tpdata->theQGSPPiK =
                               new G4QGSPPiKBuilder(quasiElasticQGS));
    tpdata->thePiK->RegisterMe(tpdata->theFTFPPiK =
                               new G4FTFPPiKBuilder(quasiElasticFTF));
    tpdata->theFTFPPiK->SetMaxEnergy(maxFTFP);
    tpdata->thePiK->RegisterMe(tpdata->theBertiniPiK =
                               new G4BertiniPiKBuilder);
    tpdata->theBertiniPiK->SetMaxEnergy(maxBERT);
    
    tpdata->theHyperon=new ExExChHyperonFTFPBuilder;
    
    tpdata->theAntiBaryon=new ExExChAntiBarionBuilder;
    tpdata->theAntiBaryon->RegisterMe(tpdata->theFTFPAntiBaryon =
                                new G4FTFPAntiBarionBuilder(quasiElasticFTF));
}

ExExChHadronPhysicsQGSP_BIC::~ExExChHadronPhysicsQGSP_BIC()
{
    if(tpdata) delete tpdata;
}

void ExExChHadronPhysicsQGSP_BIC::ConstructParticle()
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
void ExExChHadronPhysicsQGSP_BIC::ConstructProcess()
{
    if ( tpdata == 0 ) tpdata = new ThreadPrivate;
    CreateModels();
    tpdata->theNeutrons->Build();
    tpdata->thePro->Build();
    tpdata->thePiK->Build();
    tpdata->theHyperon->Build();
    tpdata->theAntiBaryon->Build();
    
    // --- Neutrons ---
    tpdata->xsNeutronInelasticXS = new G4NeutronInelasticXS();
    G4PhysListUtil::FindInelasticProcess(
            G4Neutron::Neutron())->AddDataSet(tpdata->xsNeutronInelasticXS);
    
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
    tpdata->xsNeutronCaptureXS = new G4NeutronCaptureXS();
    capture->AddDataSet(tpdata->xsNeutronCaptureXS);
    capture->RegisterMe(new G4NeutronRadCapture());
}

