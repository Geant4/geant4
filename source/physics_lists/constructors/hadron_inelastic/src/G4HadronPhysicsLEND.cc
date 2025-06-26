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
//
//---------------------------------------------------------------------------
//
// Description: Create all hadronic (neutron and photon) processes using LEND which is valid up to 20 MeV
//
// Author: Douglas M Wright, LLNL 2022-04-25 
//
//----------------------------------------------------------------------------
//

#include "G4HadronPhysicsLEND.hh"
#include "G4PhysListUtil.hh"

#include "G4NeutronBuilder.hh"
#include "G4NeutronLENDBuilder.hh"

#include "G4PhysicsListHelper.hh"
#include "G4LENDorBERTModel.hh"
//#include "G4LENDCombinedModel.hh"
#include "G4LENDCombinedCrossSection.hh"

#include "G4CrossSectionDataSetRegistry.hh"
#include "G4GammaNuclearXS.hh"

#include "G4LossTableManager.hh"
#include "G4GammaGeneralProcess.hh"

#include "G4PhysicsConstructorFactory.hh"
G4_DECLARE_PHYSCONSTR_FACTORY(G4HadronPhysicsLEND);

G4HadronPhysicsLEND::G4HadronPhysicsLEND(G4int verb, const G4String& eva)
 : G4VPhysicsConstructor("hadron inelastic LEND")
{
    verbose = verb;
    evaluation = eva;
}

void G4HadronPhysicsLEND::ConstructProcess()
{
    //....neutron-induced reactions
    
    //....create neutron processes if they do not exist
    if( G4PhysListUtil::FindInelasticProcess(G4Neutron::Neutron()) == nullptr ){
        auto neutron_processes = new G4NeutronBuilder( true ); // Fission on
        AddBuilder(neutron_processes);
        neutron_processes->Build();
    }

    //....get pointers to neutron processes so can add LEND model
    auto* neutron_inelastic = (G4HadronInelasticProcess*)G4PhysListUtil::FindInelasticProcess( G4Neutron::Neutron() );
    auto* neutron_capture   = (G4NeutronCaptureProcess*)G4PhysListUtil::FindCaptureProcess( G4Neutron::Neutron() );
    auto* neutron_fission   = (G4NeutronFissionProcess*)G4PhysListUtil::FindFissionProcess( G4Neutron::Neutron() );
        
    auto* neutronLENDBuilder = new G4NeutronLENDBuilder(evaluation);
    neutronLENDBuilder->Build(neutron_inelastic);
    neutronLENDBuilder->Build(neutron_capture);
    neutronLENDBuilder->Build(neutron_fission);

    //....photon-induced reactions
    
    //....use gamma hadronic process if it exists
    //          check Inelastic process using G4PhysListUtil
    //          then check for GammaNuclear in G4GammaGeneralProcess
    auto* gamma_inelastic = (G4HadronInelasticProcess*)G4PhysListUtil::FindInelasticProcess( G4Gamma::Gamma() );
    
    if( gamma_inelastic == nullptr){
        G4LossTableManager* emManager = G4LossTableManager::Instance();
        auto gproc = dynamic_cast<G4GammaGeneralProcess*>(emManager->GetGammaGeneralProcess());
        if (gproc != nullptr) {
            G4HadronicProcess* gnuc = gproc->GetGammaNuclear();
            if (gnuc != nullptr )
                gamma_inelastic = (G4HadronInelasticProcess*)gnuc;
        }
    }
    //....did not find existing process, so create photonuclear process
    if( gamma_inelastic == nullptr){
        //std::cout << "DMW: HadPhysLEND = gamma inelastic pointer not found\n";
        gamma_inelastic = new G4HadronInelasticProcess( "photonNuclear", G4Gamma::Gamma() );
        G4PhysicsListHelper* plHelper = G4PhysicsListHelper::GetPhysicsListHelper();
        plHelper->RegisterProcess(gamma_inelastic, G4Gamma::Gamma());

        //....need these XS or LENDonly crashes
        auto xsreg = G4CrossSectionDataSetRegistry::Instance();
        G4VCrossSectionDataSet* xs = nullptr;
        xs = xsreg->GetCrossSectionDataSet("GammaNuclearXS");
        if(nullptr == xs) xs = new G4GammaNuclearXS();
        gamma_inelastic->AddDataSet(xs);
    }

    //....add LEND photonuclear models
    auto* theGammaReactionLowE = new G4LENDorBERTModel( G4Gamma::Gamma() ); //  checks if LEND has data for specified reaction
                                                                            //  (note uses G4LENDCombinedModel)
                                                                            //  if not, uses Bertini cascade
    //auto* theGammaReactionLowE = new G4LENDCombinedModel( G4Gamma::Gamma() );  // uses LEND only
    theGammaReactionLowE->SetMaxEnergy(maxLEND_Energy);
    theGammaReactionLowE->DumpLENDTargetInfo(true);
    gamma_inelastic->RegisterMe(theGammaReactionLowE);
    
    G4LENDCombinedCrossSection* theGammaCrossSectionLowE = new G4LENDCombinedCrossSection( G4Gamma::Gamma() );
    gamma_inelastic->AddDataSet(theGammaCrossSectionLowE);
}

void G4HadronPhysicsLEND::ConstructParticle()
{
    G4Neutron::Neutron();
    G4Gamma::Gamma();
}
