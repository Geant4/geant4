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
// 080901 Add dump neutron Cross Section
//        Add Thermal Scattering by T. Koi
// 091118 Change multiple scattering processes to particle dedicated by T. Koi
// 110906 Migrate to new interface "hadr-man-V09-04-10"
//        From Process::GetMicroscopicCrossSection
//        To Process::GetElementCrossSection
//
#include <iomanip>                

#include "Tst11PhysicsList.hh"
#include "globals.hh"
#include "G4ios.hh"
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


Tst11PhysicsList::Tst11PhysicsList():  G4VUserPhysicsList()
{
  SetVerboseLevel(1);
}

Tst11PhysicsList::~Tst11PhysicsList()
{
}

void Tst11PhysicsList::ConstructParticle()
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

void Tst11PhysicsList::ConstructAllBosons()
{
  // Construct all bosons
  G4BosonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst11PhysicsList::ConstructAllLeptons()
{
  // Construct all leptons
  G4LeptonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst11PhysicsList::ConstructAllMesons()
{
  //  Construct all mesons
  G4MesonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst11PhysicsList::ConstructAllBaryons()
{
  //  Construct all barions
  G4BaryonConstructor pConstructor;
  pConstructor.ConstructParticle();
}

void Tst11PhysicsList::ConstructAllIons()
{
  //  Construct light ions
  G4IonConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst11PhysicsList::ConstructAllShortLiveds()
{
  //  Construct  resonaces and quarks
  G4ShortLivedConstructor pConstructor;
  pConstructor.ConstructParticle();  
}

void Tst11PhysicsList::ConstructProcess()
{
  AddTransportation();
  ConstructEM();
  ConstructLeptHad();
  ConstructHad();
  ConstructGeneral();
}

#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"

//#include "G4MultipleScattering.hh"
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

void Tst11PhysicsList::ConstructEM()
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
#include "G4LambdaInelasticProcess.hh"
#include "G4AntiLambdaInelasticProcess.hh"
#include "G4SigmaPlusInelasticProcess.hh"
#include "G4SigmaMinusInelasticProcess.hh"
#include "G4AntiSigmaPlusInelasticProcess.hh"
#include "G4AntiSigmaMinusInelasticProcess.hh"
#include "G4XiZeroInelasticProcess.hh"
#include "G4XiMinusInelasticProcess.hh"
#include "G4AntiXiZeroInelasticProcess.hh"
#include "G4AntiXiMinusInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4OmegaMinusInelasticProcess.hh"
#include "G4AntiOmegaMinusInelasticProcess.hh"

// Low-energy Models

#include "G4LElastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"

// Low energy neutron models
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4NeutronHPFission.hh"
#include "G4NeutronHPFissionData.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4CrossSectionDataStore.hh"

#include "G4NeutronHPThermalScattering.hh"
#include "G4NeutronHPThermalScatteringData.hh"

// ConstructHad()
//
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering, Inelastic scattering,
// Fission (for neutron only), and Capture (neutron).
//
// F.W.Jones  06-JUL-1998
//

#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryLightIonReaction.hh"

void Tst11PhysicsList::ConstructHad()
{

  //Prepairing models

  //Most hadrons
  //Bertini at low energies, then FTFP
  //theFTFPwtBERT + theBertini
   G4TheoFSGenerator * theFTFPwtBERT;
   G4TheoFSGenerator * theFTFP;
   G4PreCompoundModel * thePreEquilib;
   G4ExcitationHandler * theHandler;
   G4GeneratorPrecompoundInterface * theCascade;
   G4FTFModel * theStringModel;
   G4ExcitedStringDecay * theStringDecay;
   G4LundStringFragmentation * theLund;
   G4CascadeInterface * theBertini;

   theFTFPwtBERT = new G4TheoFSGenerator("FTFP");
   
   theFTFPwtBERT->SetMinEnergy( 2.*GeV );
   theFTFPwtBERT->SetMaxEnergy( 100.*TeV );
 
   theStringModel = new G4FTFModel;
   theStringDecay = new G4ExcitedStringDecay(theLund = new G4LundStringFragmentation);
   theStringModel->SetFragmentationModel(theStringDecay);
 
   theCascade = new G4GeneratorPrecompoundInterface;
   thePreEquilib = new G4PreCompoundModel(theHandler = new G4ExcitationHandler);
   theCascade->SetDeExcitation(thePreEquilib);  
 
   theFTFPwtBERT->SetTransport(theCascade);
   theFTFPwtBERT->SetHighEnergyGenerator(theStringModel);
   
   theBertini = new G4CascadeInterface;
   theBertini->SetMinEnergy( 0.*GeV );
   theBertini->SetMaxEnergy( 6.*GeV );
 
//AntiHyperons:
//Use FTFP for full energy range
//theFTFP
// 
   theFTFP = new G4TheoFSGenerator("FTFP");
   theFTFP->SetMinEnergy( 0.*GeV );
   theFTFP->SetMaxEnergy( 100.*TeV );
   theFTFP->SetTransport(theCascade);
   theFTFP->SetHighEnergyGenerator(theStringModel);

//Ions
//Binary Cascade + FTFP
//theIonBC + theFTFPforIon
//
   G4ExcitationHandler* handler = new G4ExcitationHandler();
   G4PreCompoundModel* thePreCompound = new G4PreCompoundModel(handler);

   G4TheoFSGenerator * theFTFPforIon;

// Binary Cascade
   G4BinaryLightIonReaction* theIonBC = new G4BinaryLightIonReaction(thePreCompound);
   theIonBC->SetMinEnergy(0.0);
   theIonBC->SetMaxEnergy(4*GeV);

// FTFP
   theFTFPforIon = theFTFPwtBERT;

//Prepairing models END

   G4HadronElasticProcess* theElasticProcess = 
                                    new G4HadronElasticProcess;
   G4LElastic* theElasticModel = new G4LElastic;
   theElasticProcess->RegisterMe(theElasticModel);
   G4HadronElasticProcess* theElasticProcess1 = 
                                    new G4HadronElasticProcess;
   theParticleIterator->reset();
   while ((*theParticleIterator)()) {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
     
      if (particleName == "pi+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionPlusInelasticProcess* theInelasticProcess = 
                                new G4PionPlusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "pi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4PionMinusInelasticProcess* theInelasticProcess = 
                                new G4PionMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonPlusInelasticProcess* theInelasticProcess = 
                                  new G4KaonPlusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0S") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroSInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroSInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon0L") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonZeroLInelasticProcess* theInelasticProcess = 
                             new G4KaonZeroLInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "kaon-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4KaonMinusInelasticProcess* theInelasticProcess = 
                                 new G4KaonMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4ProtonInelasticProcess* theInelasticProcess = 
                                    new G4ProtonInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_proton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiProtonInelasticProcess* theInelasticProcess = 
                                new G4AntiProtonInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "neutron") {
         
          // elastic scattering
         G4LElastic* theElasticModel1 = new G4LElastic;
         G4NeutronHPElastic * theElasticNeutron = new G4NeutronHPElastic;
         theElasticProcess1->RegisterMe(theElasticModel1);
         theElasticModel1->SetMinEnergy(19*MeV);
         theElasticNeutron->SetMinEnergy(4*eV);
         theElasticProcess1->RegisterMe(theElasticNeutron);
         G4NeutronHPElasticData * theNeutronData = new G4NeutronHPElasticData;
         theElasticProcess1->AddDataSet(theNeutronData);

//080901 TK add Thermal Scattering 
         G4NeutronHPThermalScattering * theThermal = new G4NeutronHPThermalScattering;
         theElasticProcess1->RegisterMe(theThermal);
         G4NeutronHPThermalScatteringData * theThermalData = new G4NeutronHPThermalScatteringData;
         theElasticProcess1->AddDataSet(theThermalData);
         pmanager->AddDiscreteProcess(theElasticProcess1);
         
          // inelastic scattering
         G4NeutronInelasticProcess* theInelasticProcess = 
                                    new G4NeutronInelasticProcess("inelastic");
         
         G4CascadeInterface* theBertiniForN = new G4CascadeInterface;
         theBertiniForN->SetMinEnergy( 19.*MeV );
         theBertiniForN->SetMaxEnergy( 6.*GeV );
         
         theInelasticProcess->RegisterMe(theBertiniForN);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         G4NeutronHPInelastic * theLENeutronInelasticModel = new G4NeutronHPInelastic;
         theInelasticProcess->RegisterMe(theLENeutronInelasticModel);
         G4NeutronHPInelasticData * theNeutronData1 = new G4NeutronHPInelasticData;
         theInelasticProcess->AddDataSet(theNeutronData1);
         pmanager->AddDiscreteProcess(theInelasticProcess);
         
          // fission
         G4HadronFissionProcess* theFissionProcess =
                                    new G4HadronFissionProcess;
         G4LFission* theFissionModel = new G4LFission;
         theFissionModel->SetMinEnergy(19*MeV);
         theFissionProcess->RegisterMe(theFissionModel);
         G4NeutronHPFission * theLENeutronFissionModel = new G4NeutronHPFission;
         theFissionProcess->RegisterMe(theLENeutronFissionModel);
         G4NeutronHPFissionData * theNeutronData2 = new G4NeutronHPFissionData;
         theFissionProcess->AddDataSet(theNeutronData2);
         pmanager->AddDiscreteProcess(theFissionProcess);
         
         // capture
         G4HadronCaptureProcess* theCaptureProcess =
                                    new G4HadronCaptureProcess;
         G4LCapture* theCaptureModel = new G4LCapture;
         theCaptureModel->SetMinEnergy(19*MeV);
         theCaptureProcess->RegisterMe(theCaptureModel);
         G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
         theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
         G4NeutronHPCaptureData * theNeutronData3 = new G4NeutronHPCaptureData;
         theCaptureProcess->AddDataSet(theNeutronData3);
         pmanager->AddDiscreteProcess(theCaptureProcess);


//080901 TK add 

         G4cout << G4endl;
         G4cout << "Cross Section Dump Through HPData->DumpPhysicsTable() " << G4endl;
         G4cout << "NeutronHP Elastic Cross Sections " << G4endl;
//         theNeutronData->DumpPhysicsTable( *G4Neutron::Neutron() );
         G4cout << "NeutronHP Inelastic Cross Sections " << G4endl;
//         theNeutronData1->DumpPhysicsTable( *G4Neutron::Neutron() );
         G4cout << "NeutronHP Capture Cross Sections " << G4endl;
//         theNeutronData2->DumpPhysicsTable( *G4Neutron::Neutron() );
         G4cout << "NeutronHP Fission Cross Sections " << G4endl;
//         theNeutronData3->DumpPhysicsTable( *G4Neutron::Neutron() );

/*
         G4cout << G4endl;
         G4cout << "Cross Section Dump Through Process::GetElementCrossSection " << G4endl;
         G4int nmat = G4Material::GetNumberOfMaterials();
         std::map< G4int , G4double > tmp_map;
         for ( G4int imat = 0 ; imat < nmat ; imat ++ )
         {
            const G4Material* mat = (*(G4Material::GetMaterialTable()))[ imat ];
            G4int ne = mat->GetNumberOfElements();
            for ( G4int iele = 0 ; iele < ne ; iele ++ )
            {

               // 110906 TK migaration to "hadr-man-V09-04-10"
               //const G4Element* ele = (*(G4Element::GetElementTable()))[ iele ];
               const G4Element* ele = mat->GetElement( iele );
               std::pair < G4int , G4double > apair( ele->GetIndex() , mat->GetTemperature());
               if ( !(tmp_map.insert ( apair )).second ) continue; 

               //G4cout << ele->GetName() << ": " << mat->GetTemperature()/kelvin << " kelvin " << mat->GetName() << G4endl;
               G4cout << ele->GetName() << ": " << mat->GetTemperature()/kelvin << " kelvin " << G4endl;
               G4cout << "Energy[eV]" << '\t' << "Elastic[mb]" << '\t' << "Inelastic[mb]" << '\t' << "Capture[mb]" << '\t' << "Fission[mb]" << G4endl;

               for ( G4int ie = 1 ; ie < 150 ; ie++ )
               {
                  G4double e = 1.0e-5*eV*std::pow( 10.0 , 1.0*ie/10 );  
                  G4DynamicParticle* dp = new G4DynamicParticle( G4Neutron::Neutron(), G4ThreeVector(1,0,0), e);
                  G4cout.precision(7);
                  G4cout << std::scientific << e/eV 
                                    << '\t' << ((G4HadronicProcess*)theElasticProcess1) ->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theInelasticProcess)->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theCaptureProcess)  ->GetElementCrossSection( dp , ele , mat )/millibarn
                                    << '\t' << ((G4HadronicProcess*)theFissionProcess)  ->GetElementCrossSection( dp , ele , mat )/millibarn
                                            << G4endl;
               } 
               // 1 - 20 MeV
//               for ( G4int ie = 1 ; ie < 20 ; ie++ )
//               {
//                  G4DynamicParticle* dp = new G4DynamicParticle( G4Neutron::Neutron(), G4ThreeVector(1,0,0), ie*MeV);
//                  G4cout << ie*MeV/eV   
//                         << " " << ((G4HadronicProcess*)theElasticProcess1)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theInelasticProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theCaptureProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << " " << ((G4HadronicProcess*)theFissionProcess)->GetMicroscopicCrossSection( dp , ele , 300*kelvin )/millibarn
//                         << G4endl;
//               }

            }
         }
         G4cout << G4endl;
*/

      }  
      else if (particleName == "anti_neutron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiNeutronInelasticProcess* theInelasticProcess = 
                               new G4AntiNeutronInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4LambdaInelasticProcess* theInelasticProcess = 
                                    new G4LambdaInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_lambda") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiLambdaInelasticProcess* theInelasticProcess = 
                                new G4AntiLambdaInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaPlusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaPlusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4SigmaMinusInelasticProcess* theInelasticProcess = 
                                 new G4SigmaMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma+") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaPlusInelasticProcess* theInelasticProcess = 
                             new G4AntiSigmaPlusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_sigma-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiSigmaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiSigmaMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiZeroInelasticProcess* theInelasticProcess = 
                            new G4XiZeroInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4XiMinusInelasticProcess* theInelasticProcess = 
                            new G4XiMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi0") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiZeroInelasticProcess* theInelasticProcess = 
                            new G4AntiXiZeroInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_xi-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiXiMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiXiMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "deuteron") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4DeuteronInelasticProcess* theInelasticProcess = 
                            new G4DeuteronInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "triton") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4TritonInelasticProcess* theInelasticProcess = 
                            new G4TritonInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "alpha") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AlphaInelasticProcess* theInelasticProcess = 
                            new G4AlphaInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theIonBC);
         theInelasticProcess->RegisterMe(theFTFPforIon);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4OmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4OmegaMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theBertini);
         theInelasticProcess->RegisterMe(theFTFPwtBERT);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
      else if (particleName == "anti_omega-") {
         pmanager->AddDiscreteProcess(theElasticProcess);
         G4AntiOmegaMinusInelasticProcess* theInelasticProcess = 
                            new G4AntiOmegaMinusInelasticProcess("inelastic");
         theInelasticProcess->RegisterMe(theFTFP);
         pmanager->AddDiscreteProcess(theInelasticProcess);
      }
   }
}

void Tst11PhysicsList::ConstructLeptHad()
{;}

#include "G4Decay.hh"
void Tst11PhysicsList::ConstructGeneral()
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

void Tst11PhysicsList::SetCuts()
{
  if (verboseLevel >0){
    G4cout << "Tst11PhysicsList::SetCuts:";
    G4cout << "CutLength : " << defaultCutValue/mm << " (mm)" << G4endl;
  }  
 //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault();   
}
