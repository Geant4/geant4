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
// $Id$
//

#include <iomanip>                

#include "Tst33PhysicsList.hh"

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4Material.hh"
#include "G4MaterialTable.hh"

Tst33PhysicsList::Tst33PhysicsList() : G4VUserPhysicsList()
{
  paraWorldName.clear();
  SetVerboseLevel(1);
}

Tst33PhysicsList::~Tst33PhysicsList()
{
  paraWorldName.clear();
}

void Tst33PhysicsList::ConstructParticle()
{
  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructAllBosons();
  ConstructAllLeptons();
  ConstructAllMesons();
  ConstructAllBaryons();
  ConstructAllIons();
  ConstructAllShortLiveds();
}

void Tst33PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst33PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst33PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst33PhysicsList::ConstructProcess()
{
  AddTransportation();
  AddScoringProcess();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"

#include "G4hIonisation.hh"

void Tst33PhysicsList::ConstructEM()
{
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
     
    if (particleName == "gamma") {
    // gamma
      // Construct processes for gamma
      pmanager->AddDiscreteProcess(new G4GammaConversion());
      pmanager->AddDiscreteProcess(new G4ComptonScattering());      
      pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());

    } else if (particleName == "e-") {
    //electron
      // Construct processes for electron
      pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4eIonisation(),-1,2,2);
      pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);
  
    } else if (particleName == "e+") {
    //positron
      // Construct processes for positron
     pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
     
     pmanager->AddProcess(new G4eIonisation(),-1,2,2);
     pmanager->AddProcess(new G4eBremsstrahlung(),-1,-1,3);      
     pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);
  
    } else if( particleName == "mu+" || 
               particleName == "mu-"    ) {
    //muon  
     // Construct processes for muon+
     pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
     pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
     pmanager->AddProcess(new G4MuBremsstrahlung(),-1,-1,3);
     pmanager->AddProcess(new G4MuPairProduction(),-1,-1,4);       
     
    } else if( particleName == "GenericIon" ) {
      pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
      pmanager->AddProcess(new G4hIonisation(),-1,2,2); 
    } else { 
      if ((particle->GetPDGCharge() != 0.0) && 
          (particle->GetParticleName() != "chargedgeantino")&&
          (!particle->IsShortLived()) ) {
     // all others charged particles except geantino
       pmanager->AddProcess(new G4hMultipleScattering(),-1,1,1);
       pmanager->AddProcess(new G4hIonisation(),-1,2,2);       
     }
    }
  }
}

// Hadron Processes

#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

#include "G4PionPlusInelasticProcess.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4KaonPlusInelasticProcess.hh"
#include "G4KaonZeroSInelasticProcess.hh"
#include "G4KaonZeroLInelasticProcess.hh"
#include "G4KaonMinusInelasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4AntiNeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

// Elastic models
#include "G4HadronElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"

// Inelastic models
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4TheoFSGenerator.hh"
#include "G4CascadeInterface.hh"
#include "G4NeutronRadCapture.hh"

// Cross sections
#include "G4VCrossSectionDataSet.hh"
#include "G4CrossSectionDataSetRegistry.hh"

#include "G4CrossSectionElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4AntiNuclElastic.hh"

#include "G4CrossSectionInelastic.hh"
#include "G4PiNuclearCrossSection.hh"
#include "G4CrossSectionPairGG.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4GGNuclNuclCrossSection.hh"
#include "G4NeutronCaptureXS.hh"

// Stopping processes
#include "G4HadronCaptureProcess.hh"
#include "G4PiMinusAbsorptionBertini.hh"
#include "G4KaonMinusAbsorptionBertini.hh"
#include "G4AntiProtonAbsorptionFritiof.hh"


void Tst33PhysicsList::ConstructHad()
{
  //Elastic models
  const G4double elastic_elimitPi = 1.0*GeV;

  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4HadronElastic* elastic_lhep1 = new G4HadronElastic();
  elastic_lhep1->SetMaxEnergy( elastic_elimitPi );
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 
  elastic_he->SetMinEnergy( elastic_elimitPi );
  
  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*GeV;
  const G4double theFTFMin1 =    4.0*GeV;
  const G4double theFTFMax =   100.0*TeV;
  const G4double theBERTMin  =   0.0*GeV;
  const G4double theBERTMax =    5.0*GeV;

  G4FTFModel * theStringModel = new G4FTFModel;
  G4ExcitedStringDecay * theStringDecay = new G4ExcitedStringDecay( new G4LundStringFragmentation );
  theStringModel->SetFragmentationModel( theStringDecay );
  G4PreCompoundModel * thePreEquilib = new G4PreCompoundModel( new G4ExcitationHandler );
  G4GeneratorPrecompoundInterface * theCascade = new G4GeneratorPrecompoundInterface( thePreEquilib );

  G4TheoFSGenerator * theFTFModel0 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel0->SetHighEnergyGenerator( theStringModel );
  theFTFModel0->SetTransport( theCascade );
  theFTFModel0->SetMinEnergy( theFTFMin0 );
  theFTFModel0->SetMaxEnergy( theFTFMax );

  G4TheoFSGenerator * theFTFModel1 = new G4TheoFSGenerator( "FTFP" );
  theFTFModel1->SetHighEnergyGenerator( theStringModel );
  theFTFModel1->SetTransport( theCascade );
  theFTFModel1->SetMinEnergy( theFTFMin1 );
  theFTFModel1->SetMaxEnergy( theFTFMax );

  G4CascadeInterface * theBERTModel = new G4CascadeInterface;
  theBERTModel->SetMinEnergy( theBERTMin );
  theBERTModel->SetMaxEnergy( theBERTMax );

  G4VCrossSectionDataSet * thePiData = new G4CrossSectionPairGG( new G4PiNuclearCrossSection, 91*GeV );
  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4VCrossSectionDataSet * theGGNuclNuclData = G4CrossSectionDataSetRegistry::Instance()->
    GetCrossSectionDataSet(G4GGNuclNuclCrossSection::Default_Name());

  theParticleIterator->reset();
  while ((*theParticleIterator)()) 
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4PionPlusInelasticProcess* theInelasticProcess = 
            new G4PionPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( thePiData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        } 

      else if (particleName == "pi-") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_lhep1 );
          theElasticProcess->RegisterMe( elastic_he );
          pmanager->AddDiscreteProcess( theElasticProcess );
          //Inelastic scattering
          G4PionMinusInelasticProcess* theInelasticProcess = 
            new G4PionMinusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( thePiData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );    
          //Absorption
          pmanager->AddRestProcess(new G4PiMinusAbsorptionBertini, ordDefault);
        }
  
      else if (particleName == "kaon+") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering       
          G4KaonPlusInelasticProcess* theInelasticProcess = 
            new G4KaonPlusInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
      
      else if (particleName == "kaon0S") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering        
          G4KaonZeroSInelasticProcess* theInelasticProcess = 
            new G4KaonZeroSInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );    
        }

      else if (particleName == "kaon0L") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonZeroLInelasticProcess* theInelasticProcess = 
            new G4KaonZeroLInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel ); 
          pmanager->AddDiscreteProcess( theInelasticProcess );    
        }

      else if (particleName == "kaon-") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4KaonMinusInelasticProcess* theInelasticProcess = 
            new G4KaonMinusInelasticProcess("inelastic");       
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->
                                           GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          //Absorption
          pmanager->AddRestProcess(new G4KaonMinusAbsorptionBertini, ordDefault);
        }

      else if (particleName == "proton") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
                                        GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->RegisterMe( elastic_chip );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4ProtonInelasticProcess* theInelasticProcess = 
            new G4ProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
 
      else if (particleName == "anti_proton") 
        {
          // Elastic scattering
          const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
          G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
          elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
          G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
          G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
          elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4AntiProtonInelasticProcess* theInelasticProcess = 
            new G4AntiProtonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
          theInelasticProcess->RegisterMe( theFTFModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );
          // Absorption
          pmanager->AddRestProcess(new G4AntiProtonAbsorptionFritiof, ordDefault);
        }

      else if (particleName == "neutron") {
        // elastic scattering
        G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->
                                       GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        theElasticProcess->RegisterMe( elastic_chip );
        pmanager->AddDiscreteProcess( theElasticProcess );
        // inelastic scattering         
        G4NeutronInelasticProcess* theInelasticProcess =
          new G4NeutronInelasticProcess("inelastic");
        theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
        theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel );
        pmanager->AddDiscreteProcess(theInelasticProcess);
        // capture
        G4HadronCaptureProcess* theCaptureProcess =
          new G4HadronCaptureProcess; 
        theCaptureProcess->AddDataSet( new G4NeutronCaptureXS() );
        G4NeutronRadCapture* theCaptureModel = new G4NeutronRadCapture;
        theCaptureProcess->RegisterMe(theCaptureModel);
        pmanager->AddDiscreteProcess(theCaptureProcess);
      }

      else if (particleName == "anti_neutron") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering (include annihilation on-fly)
          G4AntiNeutronInelasticProcess* theInelasticProcess = 
            new G4AntiNeutronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theAntiNucleonData );
          theInelasticProcess->RegisterMe( theFTFModel0 );
          pmanager->AddDiscreteProcess( theInelasticProcess );    
        }

      else if (particleName == "deuteron") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4DeuteronInelasticProcess* theInelasticProcess = 
            new G4DeuteronInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
      
      else if (particleName == "triton") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4TritonInelasticProcess* theInelasticProcess = 
            new G4TritonInelasticProcess("inelastic");
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }
 
      else if (particleName == "alpha") 
        {
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4AlphaInelasticProcess* theInelasticProcess = 
            new G4AlphaInelasticProcess("inelastic");    
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
          theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel );
          pmanager->AddDiscreteProcess( theInelasticProcess );
        }

    }
}

void Tst33PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void Tst33PhysicsList::ConstructGeneral()
{
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    if (theDecayProcess->IsApplicable(*particle)) { 
      pmanager ->AddProcess(theDecayProcess);
      pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
      pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
    }
  }
}

void Tst33PhysicsList::SetCuts()
{
  if (verboseLevel >0)
  {
    G4cout << "Tst33PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
  //   "G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}

#include "G4ParallelWorldScoringProcess.hh"
void Tst33PhysicsList::AddScoringProcess(){

  G4int npw = paraWorldName.size();
  for ( G4int i = 0; i < npw; i++){
    G4ParallelWorldScoringProcess* theParallelWorldScoringProcess
      = new G4ParallelWorldScoringProcess("ParaWorldScoringProc");
    theParallelWorldScoringProcess->SetParallelWorld(paraWorldName[i]);

    theParticleIterator->reset();
    while( (*theParticleIterator)() ){
      G4ParticleDefinition* particle = theParticleIterator->value();
      if ( !particle->IsShortLived() ){
	G4ProcessManager* pmanager = particle->GetProcessManager();
	pmanager->AddProcess(theParallelWorldScoringProcess);
	pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess,idxAtRest);
	pmanager->SetProcessOrdering(theParallelWorldScoringProcess,idxAlongStep,1);
	pmanager->SetProcessOrderingToLast(theParallelWorldScoringProcess,idxPostStep);
      }
    }
  }

}
