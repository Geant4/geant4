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
// --------------------------------------------------------------
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//
// This physics list is taken from the underground_physics example with small
// modifications.  It is an example of a "flat" physics list with no dependence
// on builders.  The physics covered would be suitable for a low background
// experiment including the neutron_hp package
//
//
//
// PhysicsList program
//
// Modified:
//
// 14-02-03 Fix bugs in msc and hIon instanciation + cut per region
// 16-08-10 Remove inclusion of obsolete class of G4ParticleWithCuts 
// 20-10-10 Migrate LowEnergy process to Livermore models, LP
// 28-03-13 Replace LEP/HEP with FTFP+BERT (A.R.)
// 15-01-20 Updated cross sections (A.R.)
// --------------------------------------------------------------

#include <iomanip>  

#include "globals.hh"
#include <CLHEP/Units/SystemOfUnits.h>
#include "G4ios.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4ProductionCutsTable.hh"

#include "G4UserLimits.hh"
#include "G4WarnPLStatus.hh"

// Builder for all stopping processes
#include "G4StoppingPhysics.hh"

#include "G4HadronicParameters.hh"
#include "G4ShortLivedConstructor.hh"
#include "LBE.hh"

// Constructor /////////////////////////////////////////////////////////////
LBE::LBE(G4int ver)
{
  if(ver > 0) {
    G4cout << "You are using the simulation engine: LBE"<<G4endl;
    G4cout <<G4endl;
  }
  defaultCutValue     = 1.0*CLHEP::micrometer; //
  cutForGamma         = defaultCutValue;
  cutForElectron      = 1.0*CLHEP::micrometer;
  cutForPositron      = defaultCutValue;
  //not used:
  // cutForProton        = defaultCutValue;
  // cutForAlpha         = 1.0*CLHEP::nanometer;
  // cutForGenericIon    = 1.0*CLHEP::nanometer;

  stoppingPhysics = new G4StoppingPhysics;

  VerboseLevel = ver;
  OpVerbLevel = 0;

  SetVerboseLevel(VerboseLevel);
}


// Destructor //////////////////////////////////////////////////////////////
LBE::~LBE() 
{
  delete stoppingPhysics;
}


// Construct Particles /////////////////////////////////////////////////////
 void LBE::ConstructParticle() 
{

  // In this method, static member functions should be called
  // for all particles which you want to use.
  // This ensures that objects of these particle types will be
  // created in the program. 

  ConstructMyBosons();
  ConstructMyLeptons();
  ConstructMyMesons();
  ConstructMyBaryons();
  ConstructMyIons();
  ConstructMyShortLiveds();
  stoppingPhysics->ConstructParticle();	// Anything not included above
}


// construct Bosons://///////////////////////////////////////////////////
 void LBE::ConstructMyBosons()
{
  // pseudo-particles
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();

  //OpticalPhotons
  G4OpticalPhoton::OpticalPhotonDefinition();
}


// construct Leptons://///////////////////////////////////////////////////
 void LBE::ConstructMyLeptons()
{
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"


// construct Mesons://///////////////////////////////////////////////////
 void LBE::ConstructMyMesons()
{
 //  mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

}


// construct Baryons://///////////////////////////////////////////////////
 void LBE::ConstructMyBaryons()
{
 //  baryons
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

}


// construct Ions://///////////////////////////////////////////////////
 void LBE::ConstructMyIons()
{
 //  ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();

}

// construct Shortliveds://///////////////////////////////////////////////////
 void LBE::ConstructMyShortLiveds()
{
  // ShortLiveds
  G4ShortLivedConstructor pShortLivedConstructor;
  pShortLivedConstructor.ConstructParticle();
}




// Construct Processes //////////////////////////////////////////////////////
 void LBE::ConstructProcess() 
{
  AddTransportation();
  ConstructEM();
  ConstructOp();
  ConstructHad();
  ConstructGeneral();
}


// Transportation ///////////////////////////////////////////////////////////
#include "G4MaxTimeCuts.hh"
#include "G4MinEkineCuts.hh"

 void LBE::AddTransportation() {
  
  G4VUserPhysicsList::AddTransportation();
  
  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*(myParticleIterator))() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // time cuts for ONLY neutrons:
    if(particleName == "neutron") 
    pmanager->AddDiscreteProcess(new G4MaxTimeCuts());
    // Energy cuts to kill charged (embedded in method) particles:
    pmanager->AddDiscreteProcess(new G4MinEkineCuts());
  }		      
}


// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles

#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"

// gamma. Use Livermore models
#include "G4PhotoElectricEffect.hh"
#include "G4LivermorePhotoElectricModel.hh"
#include "G4ComptonScattering.hh"
#include "G4LivermoreComptonModel.hh"
#include "G4GammaConversion.hh"
#include "G4LivermoreGammaConversionModel.hh"
#include "G4RayleighScattering.hh" 
#include "G4LivermoreRayleighModel.hh"


// e-
#include "G4eMultipleScattering.hh"
#include "G4UniversalFluctuation.hh"
#include "G4UrbanMscModel.hh"

#include "G4eIonisation.hh"
#include "G4LivermoreIonisationModel.hh"

#include "G4eBremsstrahlung.hh"
#include "G4LivermoreBremsstrahlungModel.hh"

// e+
#include "G4eplusAnnihilation.hh"


// alpha and GenericIon and deuterons, triton, He3:
#include "G4ionIonisation.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
//
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
#include "G4EnergyLossTables.hh"

//muon:
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCapture.hh"

//OTHERS:
//#include "G4hIonisation.hh" // standard hadron ionisation

 void LBE::ConstructEM() {

 // models & processes:
 // Use Livermore models up to 20 MeV, and standard 
 // models for higher energy
 G4double LivermoreHighEnergyLimit = 20*CLHEP::MeV;
 //
  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*(myParticleIterator))() ){
    G4ParticleDefinition* particle = myParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();
    
    if (particleName == "gamma") 
      {
      G4PhotoElectricEffect* thePhotoElectricEffect = new G4PhotoElectricEffect();
      G4LivermorePhotoElectricModel* theLivermorePhotoElectricModel = 
	new G4LivermorePhotoElectricModel();
      theLivermorePhotoElectricModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      thePhotoElectricEffect->AddEmModel(0, theLivermorePhotoElectricModel);
      pmanager->AddDiscreteProcess(thePhotoElectricEffect);

      G4ComptonScattering* theComptonScattering = new G4ComptonScattering();
      G4LivermoreComptonModel* theLivermoreComptonModel = 
	new G4LivermoreComptonModel();
      theLivermoreComptonModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theComptonScattering->AddEmModel(0, theLivermoreComptonModel);
      pmanager->AddDiscreteProcess(theComptonScattering);

      G4GammaConversion* theGammaConversion = new G4GammaConversion();
      G4LivermoreGammaConversionModel* theLivermoreGammaConversionModel = 
	new G4LivermoreGammaConversionModel();
      theLivermoreGammaConversionModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theGammaConversion->AddEmModel(0, theLivermoreGammaConversionModel);
      pmanager->AddDiscreteProcess(theGammaConversion);

      G4RayleighScattering* theRayleigh = new G4RayleighScattering();
      G4LivermoreRayleighModel* theRayleighModel = new G4LivermoreRayleighModel();
      theRayleighModel->SetHighEnergyLimit(LivermoreHighEnergyLimit);
      theRayleigh->AddEmModel(0, theRayleighModel);
      pmanager->AddDiscreteProcess(theRayleigh);

      } 
    else if (particleName == "e-") 
      {
       //electron
       // process ordering: AddProcess(name, at rest, along step, post step)
       // -1 = not implemented, then ordering
       	G4eMultipleScattering* msc = new G4eMultipleScattering();     
        //msc->AddEmModel(0, new G4UrbanMscModel());
        msc->SetStepLimitType(fUseDistanceToBoundary);
        pmanager->AddProcess(msc,                   -1, 1, 1);
      
       // Ionisation
       G4eIonisation* eIoni = new G4eIonisation();
       G4LivermoreIonisationModel* theIoniLivermore = new
        G4LivermoreIonisationModel();
       theIoniLivermore->SetHighEnergyLimit(1*CLHEP::MeV); 
       eIoni->AddEmModel(0, theIoniLivermore, new G4UniversalFluctuation() );
       eIoni->SetStepFunction(0.2, 100*CLHEP::um); //     
       pmanager->AddProcess(eIoni,                 -1, 2, 2);
      
       // Bremsstrahlung
       G4eBremsstrahlung* eBrem = new G4eBremsstrahlung();
       G4LivermoreBremsstrahlungModel* theBremLivermore = new
         G4LivermoreBremsstrahlungModel();
       theBremLivermore->SetHighEnergyLimit(LivermoreHighEnergyLimit);
       eBrem->AddEmModel(0, theBremLivermore);
       pmanager->AddProcess(eBrem, -1,-3, 3);	
      } 
    else if (particleName == "e+") 
      {
	//positron
      G4eMultipleScattering* msc = new G4eMultipleScattering();
      //msc->AddEmModel(0, new G4UrbanMscModel());      
      msc->SetStepLimitType(fUseDistanceToBoundary);
      pmanager->AddProcess(msc,                   -1, 1, 1);
      G4eIonisation* eIoni = new G4eIonisation();
      eIoni->SetStepFunction(0.2, 100*CLHEP::um);      
      pmanager->AddProcess(eIoni,                 -1, 2, 2);
      pmanager->AddProcess(new G4eBremsstrahlung, -1,-3, 3);      
      pmanager->AddProcess(new G4eplusAnnihilation,0,-1, 4);
      } 
    else if( particleName == "mu+" || 
	     particleName == "mu-"    ) 
      {
	//muon  
        G4MuMultipleScattering* aMultipleScattering = new G4MuMultipleScattering();
	pmanager->AddProcess(aMultipleScattering,           -1, 1, 1);
	pmanager->AddProcess(new G4MuIonisation(),          -1, 2, 2);
	pmanager->AddProcess(new G4MuBremsstrahlung(),      -1,-1, 3);
	pmanager->AddProcess(new G4MuPairProduction(),      -1,-1, 4);
	if( particleName == "mu-" )
	  pmanager->AddProcess(new G4MuonMinusCapture(), 0,-1,-1);
      } 
    else if (particleName == "GenericIon")
    {
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetEmModel(new G4IonParametrisedLossModel());
      ionIoni->SetStepFunction(0.1, 10*CLHEP::um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);	
    }
    else if (particleName == "alpha" || particleName == "He3")
    {
      //MSC, ion-Ionisation, Nuclear Stopping
      pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);

      G4ionIonisation* ionIoni = new G4ionIonisation();
      ionIoni->SetStepFunction(0.1, 20*CLHEP::um);
      pmanager->AddProcess(ionIoni,                   -1, 2, 2);
      pmanager->AddProcess(new G4NuclearStopping(),   -1, 3,-1);
    }
    else if (particleName == "proton"     ||	  
	     particleName == "deuteron"   ||
	     particleName == "triton"     ||
             particleName == "pi+" ||
             particleName == "pi-" ||
	     particleName == "kaon+" ||
             particleName == "kaon-") 
      {
       //MSC, h-ionisation, bremsstrahlung
       pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);      
       G4hIonisation* hIoni = new G4hIonisation();
       hIoni->SetStepFunction(0.2, 50*CLHEP::um);
       pmanager->AddProcess(hIoni,                     -1, 2, 2);      
       pmanager->AddProcess(new G4hBremsstrahlung,     -1,-3, 3);    
      } 
    else if ((!particle->IsShortLived()) &&
	     (charge != 0.0) && 
	     (particle->GetParticleName() != "chargedgeantino")) 
      {
	//all others charged particles except geantino
        pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
        pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
      }
    
  }
}


// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
//#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

 void LBE::ConstructOp() 
{
  // default scintillation process
  //Coverity report: check that the process is actually used, if not must delete
  G4bool theScintProcessDefNeverUsed = true;
  G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
  // theScintProcessDef->DumpPhysicsTable();
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetVerboseLevel(OpVerbLevel);

  // scintillation process for alpha:
  G4bool theScintProcessAlphaNeverUsed = true;
  G4Scintillation* theScintProcessAlpha = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessAlpha->SetTrackSecondariesFirst(true);
  theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

  // scintillation process for heavy nuclei
  G4bool theScintProcessNucNeverUsed = true;  
  G4Scintillation* theScintProcessNuc = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessNuc->SetTrackSecondariesFirst(true);
  theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

  // optical processes
  G4bool theAbsorptionProcessNeverUsed = true;
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  //  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4bool theBoundaryProcessNeverUsed = true;
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  //  theAbsorptionProcess->DumpPhysicsTable();
  //  theRayleighScatteringProcess->DumpPhysicsTable();
  theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
  // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
  theBoundaryProcess->SetVerboseLevel(OpVerbLevel);

  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*(myParticleIterator))() )
    {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (theScintProcessDef->IsApplicable(*particle)) {
	//      if(particle->GetPDGMass() > 5.0*CLHEP::GeV) 
	if(particle->GetParticleName() == "GenericIon") {
	  pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxPostStep);
	  theScintProcessNucNeverUsed = false;
	}	  
	else if(particle->GetParticleName() == "alpha") {
	  pmanager->AddProcess(theScintProcessAlpha);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxPostStep);
	  theScintProcessAlphaNeverUsed = false;
	}
	else {
	  pmanager->AddProcess(theScintProcessDef);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	  theScintProcessDefNeverUsed = false;
	}	  
      }
      
      if (particleName == "opticalphoton") {
	pmanager->AddDiscreteProcess(theAbsorptionProcess);
	theAbsorptionProcessNeverUsed = false;
	//	pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	theBoundaryProcessNeverUsed = false;
	pmanager->AddDiscreteProcess(theBoundaryProcess);
      }
    }
    if ( theScintProcessDefNeverUsed ) delete theScintProcessDef;
    if ( theScintProcessAlphaNeverUsed ) delete theScintProcessAlpha;
    if ( theScintProcessNucNeverUsed ) delete theScintProcessNuc;
    if ( theBoundaryProcessNeverUsed ) delete theBoundaryProcess;
    if ( theAbsorptionProcessNeverUsed ) delete theAbsorptionProcess;
}


// Hadronic processes ////////////////////////////////////////////////////////

// Elastic processes:
#include "G4HadronElasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4HadronElastic.hh"
#include "G4ChipsElasticModel.hh"
#include "G4ElasticHadrNucleusHE.hh"
#include "G4AntiNuclElastic.hh"
#include "G4BGGPionElasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4ChipsProtonElasticXS.hh"
#include "G4ChipsNeutronElasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"  
#include "G4ChipsKaonMinusElasticXS.hh"
#include "G4ChipsKaonPlusElasticXS.hh"
#include "G4ChipsKaonZeroElasticXS.hh"
#include "G4BGGNucleonElasticXS.hh"
#include "G4CrossSectionElastic.hh"

// Inelastic processes:
#include "G4HadronInelasticProcess.hh"

// FTFP + BERT model
#include "G4TheoFSGenerator.hh"
#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4GeneratorPrecompoundInterface.hh"
#include "G4FTFModel.hh"
#include "G4LundStringFragmentation.hh"
#include "G4ExcitedStringDecay.hh"
#include "G4CascadeInterface.hh"
#include "G4CrossSectionInelastic.hh"
#include "G4BGGPionInelasticXS.hh"
#include "G4ChipsKaonMinusInelasticXS.hh"
#include "G4ChipsKaonPlusInelasticXS.hh"
#include "G4ChipsKaonZeroInelasticXS.hh"
#include "G4CrossSectionDataSetRegistry.hh"
#include "G4BGGNucleonInelasticXS.hh"
#include "G4ComponentAntiNuclNuclearXS.hh"
#include "G4ComponentGGNuclNuclXsc.hh"

// Neutron high-precision models: <20 MeV
#include "G4ParticleHPElastic.hh"
#include "G4ParticleHPElasticData.hh"
#include "G4ParticleHPCapture.hh"
#include "G4ParticleHPCaptureData.hh"
#include "G4ParticleHPInelastic.hh"
#include "G4ParticleHPInelasticData.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4NeutronRadCapture.hh"

// Binary light ion cascade for alpha, deuteron and triton
#include "G4BinaryLightIonReaction.hh"

// ConstructHad()
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
// F.W.Jones  09-JUL-1998
 void LBE::ConstructHad() 
{
  // Elastic scattering
  G4HadronElastic* elastic_lhep0 = new G4HadronElastic();
  G4ChipsElasticModel* elastic_chip = new G4ChipsElasticModel();
  G4ElasticHadrNucleusHE* elastic_he = new G4ElasticHadrNucleusHE(); 

  const G4double elastic_elimitAntiNuc = 100.0*CLHEP::MeV;
  G4AntiNuclElastic* elastic_anuc = new G4AntiNuclElastic();
  elastic_anuc->SetMinEnergy( elastic_elimitAntiNuc );
  G4CrossSectionElastic* elastic_anucxs = new G4CrossSectionElastic( elastic_anuc->GetComponentCrossSection() );
  G4HadronElastic* elastic_lhep2 = new G4HadronElastic();
  elastic_lhep2->SetMaxEnergy( elastic_elimitAntiNuc );

  // Inelastic scattering
  const G4double theFTFMin0 =    0.0*CLHEP::GeV;
  const G4double theFTFMin1 =    4.0*CLHEP::GeV;
  const G4double theFTFMax = G4HadronicParameters::Instance()->GetMaxEnergy();
  const G4double theBERTMin0 =   0.0*CLHEP::GeV;
  const G4double theBERTMin1 =  19.0*CLHEP::MeV;
  const G4double theBERTMax =    5.0*CLHEP::GeV;
  const G4double theHPMin =      0.0*CLHEP::GeV;
  const G4double theHPMax =     20.0*CLHEP::MeV;
  const G4double theIonBCMin =   0.0*CLHEP::GeV;
  const G4double theIonBCMax =   5.0*CLHEP::GeV;


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

  G4CascadeInterface * theBERTModel0 = new G4CascadeInterface;
  theBERTModel0->SetMinEnergy( theBERTMin0 );
  theBERTModel0->SetMaxEnergy( theBERTMax );

  G4CascadeInterface * theBERTModel1 = new G4CascadeInterface;
  theBERTModel1->SetMinEnergy( theBERTMin1 );
  theBERTModel1->SetMaxEnergy( theBERTMax );

  // Binary Cascade
  G4BinaryLightIonReaction * theIonBC = new G4BinaryLightIonReaction( thePreEquilib );
  theIonBC->SetMinEnergy( theIonBCMin );
  theIonBC->SetMaxEnergy( theIonBCMax );

  G4VCrossSectionDataSet * theAntiNucleonData = new G4CrossSectionInelastic( new G4ComponentAntiNuclNuclearXS );
  G4ComponentGGNuclNuclXsc * ggNuclNuclXsec = new G4ComponentGGNuclNuclXsc();
  G4VCrossSectionDataSet * theGGNuclNuclData = new G4CrossSectionInelastic(ggNuclNuclXsec);

  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while ((*(myParticleIterator))()) 
    {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4PionPlus::Definition() );
          theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( G4PionPlus::Definition() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	} 

      else if (particleName == "pi-") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( new G4BGGPionElasticXS( particle ) );
          theElasticProcess->RegisterMe( elastic_he );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4PionMinus::Definition() );	  
          theInelasticProcess->AddDataSet( new G4BGGPionInelasticXS( G4PionMinus::Definition() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon+") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4KaonPlus::Definition() );
	  
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonPlusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "kaon0S") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4KaonZeroShort::Definition() );
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon0L") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
          //G4KaonZeroLInelasticProcess* theInelasticProcess = new G4KaonZeroLInelasticProcess("inelastic");
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4KaonZeroLong::Definition() );
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonZeroInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 ); 
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "kaon-") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusElasticXS::Default_Name()));
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4KaonMinus::Definition() );
          theInelasticProcess->AddDataSet( G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsKaonMinusInelasticXS::Default_Name()));
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "proton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsProtonElasticXS::Default_Name()));
          theElasticProcess->AddDataSet( new G4BGGNucleonElasticXS( G4Proton::Proton() ) );
          theElasticProcess->RegisterMe( elastic_chip );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4Proton::Definition() );
          theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Proton::Proton() ) );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
          theInelasticProcess->RegisterMe( theBERTModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "anti_proton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4AntiProton::Definition() );
          theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
        theElasticProcess->AddDataSet(G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4ChipsNeutronElasticXS::Default_Name()));
        G4HadronElastic* elastic_neutronChipsModel = new G4ChipsElasticModel();
	elastic_neutronChipsModel->SetMinEnergy( 19.0*CLHEP::MeV );
        theElasticProcess->RegisterMe( elastic_neutronChipsModel );
	G4ParticleHPElastic * theElasticNeutronHP = new G4ParticleHPElastic;
        theElasticNeutronHP->SetMinEnergy( theHPMin );
        theElasticNeutronHP->SetMaxEnergy( theHPMax );
	theElasticProcess->RegisterMe( theElasticNeutronHP );
	theElasticProcess->AddDataSet( new G4ParticleHPElasticData );
	pmanager->AddDiscreteProcess( theElasticProcess );
	// inelastic scattering
	G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4Neutron::Definition() );
        theInelasticProcess->AddDataSet( new G4BGGNucleonInelasticXS( G4Neutron::Neutron() ) );
	theInelasticProcess->RegisterMe( theFTFModel1 );
        theInelasticProcess->RegisterMe( theBERTModel1 );
	G4ParticleHPInelastic * theNeutronInelasticHPModel = new G4ParticleHPInelastic;
        theNeutronInelasticHPModel->SetMinEnergy( theHPMin );
        theNeutronInelasticHPModel->SetMaxEnergy( theHPMax );
	theInelasticProcess->RegisterMe( theNeutronInelasticHPModel );
	theInelasticProcess->AddDataSet( new G4ParticleHPInelasticData );
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4NeutronCaptureProcess* theCaptureProcess = new G4NeutronCaptureProcess;
	G4ParticleHPCapture * theNeutronCaptureHPModel = new G4ParticleHPCapture;
        theNeutronCaptureHPModel->SetMinEnergy( theHPMin );
        theNeutronCaptureHPModel->SetMaxEnergy( theHPMax );
	G4NeutronRadCapture* theNeutronRadCapture = new G4NeutronRadCapture(); 
 	theNeutronRadCapture->SetMinEnergy(theHPMax*0.99); 
	theCaptureProcess->RegisterMe( theNeutronCaptureHPModel );
	theCaptureProcess->RegisterMe( theNeutronRadCapture);
	theCaptureProcess->AddDataSet( new G4ParticleHPCaptureData );
	theCaptureProcess->AddDataSet((G4NeutronCaptureXS*)G4CrossSectionDataSetRegistry::Instance()->GetCrossSectionDataSet(G4NeutronCaptureXS::Default_Name()));
	pmanager->AddDiscreteProcess(theCaptureProcess);
      }
      else if (particleName == "anti_neutron") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->AddDataSet( elastic_anucxs );
          theElasticProcess->RegisterMe( elastic_lhep2 );
          theElasticProcess->RegisterMe( elastic_anuc );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4AntiNeutron::Definition() );
          theInelasticProcess->AddDataSet( theAntiNucleonData );
	  theInelasticProcess->RegisterMe( theFTFModel0 );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "deuteron") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4Deuteron::Definition() );	  
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
      
      else if (particleName == "triton") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4Triton::Definition() );	  
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}

      else if (particleName == "alpha") 
	{
          // Elastic scattering
          G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
          theElasticProcess->RegisterMe( elastic_lhep0 );
          theElasticProcess->AddDataSet( theGGNuclNuclData );
	  pmanager->AddDiscreteProcess( theElasticProcess );
          // Inelastic scattering
	  G4HadronInelasticProcess* theInelasticProcess = new G4HadronInelasticProcess( "inelastic", G4Alpha::Definition() );	  
          theInelasticProcess->AddDataSet( theGGNuclNuclData );
	  theInelasticProcess->RegisterMe( theFTFModel1 );
	  theInelasticProcess->RegisterMe( theIonBC );
	  pmanager->AddDiscreteProcess( theInelasticProcess );
	}
    }	// while ((*(myParticleIterator))()) 

  // Add stopping processes with builder
  stoppingPhysics->ConstructProcess();
}


// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

#include "G4LossTableManager.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4NuclearLevelData.hh"
#include "G4NuclideTable.hh"

 void LBE::ConstructGeneral() {

  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  G4bool theDecayProcessNeverUsed = true; //Check if theDecayProcess will be used
  auto myParticleIterator=G4ParticleTable::GetParticleTable()->GetIterator();
  myParticleIterator->reset();
  while( (*(myParticleIterator))() )
    {
      G4ParticleDefinition* particle = myParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      
      if (theDecayProcess->IsApplicable(*particle) && !particle->IsShortLived()) 
	{ 
	  theDecayProcessNeverUsed = false;
	  pmanager ->AddProcess(theDecayProcess);
	  // set ordering for PostStepDoIt and AtRestDoIt
	  pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	  pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
	}
    }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable = 
    G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay* theRadioactiveDecay = new G4RadioactiveDecay();
 
  //Fix for activation of RadioactiveDecay, based on G4RadioactiveDecayPhysics 
  G4EmParameters* param = G4EmParameters::Instance();
  param->SetAugerCascade(true);
  param->AddPhysics("world","G4RadioactiveDecay");

  G4DeexPrecoParameters* deex = G4NuclearLevelData::GetInstance()->GetParameters();
  deex->SetStoreAllLevels(true);
  deex->SetMaxLifeTime(G4NuclideTable::GetInstance()->GetThresholdOfHalfLife()
                       /std::log(2.));

  G4LossTableManager* man = G4LossTableManager::Instance();
  G4VAtomDeexcitation* ad = man->AtomDeexcitation();
  if(!ad) {
    ad = new G4UAtomicDeexcitation();
    man->SetAtomDeexcitation(ad);
    ad->InitialiseAtomicDeexcitation();
  }

  for (G4int i=0; i<theIonTable->Entries(); i++) 
    {
      G4String particleName = theIonTable->GetParticle(i)->GetParticleName();
      G4String particleType = theIonTable->GetParticle(i)->GetParticleType();
      
      if (particleName == "GenericIon") 
	{
	  G4ProcessManager* pmanager = 
	    theIonTable->GetParticle(i)->GetProcessManager();
	  pmanager->SetVerboseLevel(VerboseLevel);
	  pmanager ->AddProcess(theRadioactiveDecay);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxPostStep);
	  pmanager ->SetProcessOrdering(theRadioactiveDecay, idxAtRest);
	} 
    }
    //If we actually never used the process, delete it
    //From Coverity report
    if ( theDecayProcessNeverUsed ) delete theDecayProcess;
}

// Cuts /////////////////////////////////////////////////////////////////////
void LBE::SetCuts() 
{
  
  if (verboseLevel >1)
    G4cout << "LBE::SetCuts:";
  
  if (verboseLevel>0){
    G4cout << "LBE::SetCuts:";
    G4cout << "CutLength : " 
	   << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  //special for low energy physics
  G4double lowlimit=250*CLHEP::eV;  
  G4ProductionCutsTable * aPCTable = G4ProductionCutsTable::GetProductionCutsTable();
  aPCTable->SetEnergyRange(lowlimit,100*CLHEP::GeV);
   
  // set cut values for gamma at first and for e- second and next for e+,
  // because some processes for e+/e- need cut values for gamma 
  SetCutValue(cutForGamma, "gamma");
  SetCutValue(cutForElectron, "e-");
  SetCutValue(cutForPositron, "e+");
  
  //  SetCutValue(cutForProton, "proton");
  //  SetCutValue(cutForProton, "anti_proton");
  //  SetCutValue(cutForAlpha,  "alpha");
  //  SetCutValue(cutForGenericIon,  "GenericIon");
  
  //  SetCutValueForOthers(defaultCutValue);
  
  if (verboseLevel>0) DumpCutValuesTable();
}

