//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PhysicsList program
//
// Modified:
//
// 14-02-03 Fix bugs in msc and hIon instanciation + cut per region
//
// --------------------------------------------------------------

#include "DMXPhysicsList.hh"

#include "globals.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleWithCuts.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"

#include "G4ios.hh"
#include "g4std/iomanip"  

#include "G4UserLimits.hh"


// Constructor /////////////////////////////////////////////////////////////
DMXPhysicsList::DMXPhysicsList() : G4VUserPhysicsList() 
{

  defaultCutValue     = 1.0*micrometer; //
  cutForGamma         = defaultCutValue;
  cutForElectron      = 1.0*nanometer;
  cutForPositron      = defaultCutValue;
  cutForProton        = defaultCutValue;
  cutForAlpha         = 1.0*nanometer;
  cutForGenericIon    = 1.0*nanometer;

  VerboseLevel = 1;
  OpVerbLevel = 0;

  SetVerboseLevel(VerboseLevel);
}


// Destructor //////////////////////////////////////////////////////////////
DMXPhysicsList::~DMXPhysicsList() 
{;}


// Construct Particles /////////////////////////////////////////////////////
void DMXPhysicsList::ConstructParticle() 
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

}


// construct Bosons://///////////////////////////////////////////////////
void DMXPhysicsList::ConstructMyBosons()
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
void DMXPhysicsList::ConstructMyLeptons()
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


// construct Mesons://///////////////////////////////////////////////////
void DMXPhysicsList::ConstructMyMesons()
{
 //  mesons
  G4PionPlus::PionPlusDefinition();
  G4PionMinus::PionMinusDefinition();
  G4PionZero::PionZeroDefinition();
  G4KaonPlus::KaonPlusDefinition();
  G4KaonMinus::KaonMinusDefinition();
  G4Eta::EtaDefinition();
  G4EtaPrime::EtaPrimeDefinition();
  G4KaonZero::KaonZeroDefinition();
  G4AntiKaonZero::AntiKaonZeroDefinition();
  G4KaonZeroLong::KaonZeroLongDefinition();
  G4KaonZeroShort::KaonZeroShortDefinition();
}


// construct Baryons://///////////////////////////////////////////////////
void DMXPhysicsList::ConstructMyBaryons()
{
//  barions
  G4Proton::ProtonDefinition();
  G4AntiProton::AntiProtonDefinition();
  G4Neutron::NeutronDefinition();
  G4AntiNeutron::AntiNeutronDefinition();
}


// construct Ions://///////////////////////////////////////////////////
void DMXPhysicsList::ConstructMyIons()
{
//  Ions
  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4He3::He3Definition();
  G4Alpha::AlphaDefinition();
  G4GenericIon::GenericIonDefinition();
}

// construct Shortliveds://///////////////////////////////////////////////////
void DMXPhysicsList::ConstructMyShortLiveds()
{
  // ShortLiveds
  ;
}




// Construct Processes //////////////////////////////////////////////////////
void DMXPhysicsList::ConstructProcess() 
{

  AddTransportation();

  ConstructEM();

  ConstructOp();

  ConstructHad();

  ConstructGeneral();

}


// Transportation ///////////////////////////////////////////////////////////
#include "DMXMaxTimeCuts.hh"
#include "DMXMinEkineCuts.hh"

void DMXPhysicsList::AddTransportation() {
  
  G4VUserPhysicsList::AddTransportation();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // time cuts for ONLY neutrons:
    if(particleName == "neutron") 
    pmanager->AddDiscreteProcess(new DMXMaxTimeCuts());
    // Energy cuts to kill charged (embedded in method) particles:
    pmanager->AddDiscreteProcess(new DMXMinEkineCuts());
  }		      
}


// Electromagnetic Processes ////////////////////////////////////////////////
// all charged particles
#include "G4MultipleScattering.hh"

// gamma
#include "G4LowEnergyRayleigh.hh" 
#include "G4LowEnergyPhotoElectric.hh"
#include "G4LowEnergyCompton.hh"  
#include "G4LowEnergyGammaConversion.hh" 


// e-
#include "G4LowEnergyIonisation.hh" 
#include "G4LowEnergyBremsstrahlung.hh" 

// e+
#include "G4eIonisation.hh" 
#include "G4eBremsstrahlung.hh" 
#include "G4eplusAnnihilation.hh"


// alpha and GenericIon and deuterons, triton, He3:
#include "G4hLowEnergyIonisation.hh"
#include "G4EnergyLossTables.hh"
// hLowEnergyIonisation uses Ziegler 1988 as the default


//muon:
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4MuonMinusCaptureAtRest.hh"

//OTHERS:
//#include "G4hIonisation.hh" // standard hadron ionisation

void DMXPhysicsList::ConstructEM() {

// processes:

  G4LowEnergyPhotoElectric* lowePhot = new G4LowEnergyPhotoElectric();
  G4LowEnergyIonisation* loweIon  = new G4LowEnergyIonisation();
  G4LowEnergyBremsstrahlung* loweBrem = new G4LowEnergyBremsstrahlung();

  // note LowEIon uses proton as basis for its data-base, therefore
  // cannot specify different LowEnergyIonisation models for different
  // particles, but can change model globally for Ion, Alpha and Proton.


  //fluorescence apply specific cut for fluorescence from photons, electrons
  //and bremsstrahlung photons:
  G4double fluorcut = 250*eV;
  lowePhot->SetCutForLowEnSecPhotons(fluorcut);
  loweIon->SetCutForLowEnSecPhotons(fluorcut);
  loweBrem->SetCutForLowEnSecPhotons(fluorcut);
  
  // setting tables explicitly for electronic stopping power
  //  ahadronLowEIon->SetElectronicStoppingPowerModel
  //  (G4GenericIon::GenericIonDefinition(), "ICRU_R49p") ;
  //  ahadronLowEIon->SetElectronicStoppingPowerModel
  //  (G4Proton::ProtonDefinition(), "ICRU_R49p") ;

  // Switch off the Barkas and Bloch corrections
  //  ahadronLowEIon->SetBarkasOff();


  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition* particle = theParticleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    G4String particleType = particle->GetParticleType();
    G4double charge = particle->GetPDGCharge();
    
    if (particleName == "gamma") 
      {
	//gamma
	pmanager->AddDiscreteProcess(new G4LowEnergyRayleigh());
	pmanager->AddDiscreteProcess(lowePhot);
	pmanager->AddDiscreteProcess(new G4LowEnergyCompton());
	pmanager->AddDiscreteProcess(new G4LowEnergyGammaConversion());
      } 
    else if (particleName == "e-") 
      {
	//electron
	// process ordering: AddProcess(name, at rest, along step, post step)
	// -1 = not implemented, then ordering
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
	pmanager->AddProcess(aMultipleScattering,     -1, 1, 1);
	pmanager->AddProcess(loweIon,                 -1, 2, 2);
	pmanager->AddProcess(loweBrem,                -1,-1, 3);
      } 
    else if (particleName == "e+") 
      {
	//positron
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
	pmanager->AddProcess(aMultipleScattering,     -1, 1, 1);
	pmanager->AddProcess(new G4eIonisation(),     -1, 2, 2);
	pmanager->AddProcess(new G4eBremsstrahlung(), -1,-1, 3);
	pmanager->AddProcess(new G4eplusAnnihilation(),0,-1, 4);      
      } 
    else if( particleName == "mu+" || 
	     particleName == "mu-"    ) 
      {
	//muon  
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
	pmanager->AddProcess(aMultipleScattering,           -1, 1, 1);
	pmanager->AddProcess(new G4MuIonisation(),          -1, 2, 2);
	pmanager->AddProcess(new G4MuBremsstrahlung(),      -1,-1, 3);
	pmanager->AddProcess(new G4MuPairProduction(),      -1,-1, 4);
	if( particleName == "mu-" )
	  pmanager->AddProcess(new G4MuonMinusCaptureAtRest(), 0,-1,-1);
      } 
    else if (particleName == "proton"     ||
	     particleName == "alpha"      ||
	     particleName == "deuteron"   ||
	     particleName == "triton"     ||
	     particleName == "He3"        ||
	     particleName == "GenericIon" || 
	     (particleType == "nucleus" && charge != 0)) 
      {
	// OBJECT may be dynamically created as either a GenericIon or nucleus
	// G4Nucleus exists and therefore has particle type nucleus
	// genericIon:
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
        G4hLowEnergyIonisation* ahadronLowEIon = new G4hLowEnergyIonisation();
	pmanager->AddProcess(aMultipleScattering,-1,1,1);
	pmanager->AddProcess(ahadronLowEIon,-1,2,2); 
        // ahadronLowEIon->SetNuclearStoppingOff() ;
	//        ahadronLowEIon->SetNuclearStoppingPowerModel("ICRU_R49") ;
	//        ahadronLowEIon->SetNuclearStoppingOn() ;
  
        //fluorescence switch off for hadrons (for now) PIXE:
        ahadronLowEIon->SetFluorescence(false);
      } 
    else if ((!particle->IsShortLived()) &&
	     (charge != 0.0) && 
	     (particle->GetParticleName() != "chargedgeantino")) 
      {
	//all others charged particles except geantino
        G4MultipleScattering* aMultipleScattering = new G4MultipleScattering();
        G4hLowEnergyIonisation* ahadronLowEIon = new G4hLowEnergyIonisation();
	pmanager->AddProcess(aMultipleScattering,-1,1,1);
	pmanager->AddProcess(ahadronLowEIon,       -1,2,2);      
	//      pmanager->AddProcess(new G4hIonisation(),       -1,2,2);      
      }
    
  }
}


// Optical Processes ////////////////////////////////////////////////////////
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
//#include "G4OpRayleigh.hh"
#include "G4OpBoundaryProcess.hh"

void DMXPhysicsList::ConstructOp() 
{
  // default scintillation process
  G4Scintillation* theScintProcessDef = new G4Scintillation("Scintillation");
  // theScintProcessDef->DumpPhysicsTable();
  theScintProcessDef->SetTrackSecondariesFirst(true);
  theScintProcessDef->SetScintillationYieldFactor(1.0); //
  theScintProcessDef->SetScintillationExcitationRatio(0.0); //
  theScintProcessDef->SetVerboseLevel(OpVerbLevel);

  // scintillation process for alpha:
  G4Scintillation* theScintProcessAlpha = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessAlpha->SetTrackSecondariesFirst(true);
  theScintProcessAlpha->SetScintillationYieldFactor(1.1);
  theScintProcessAlpha->SetScintillationExcitationRatio(1.0);
  theScintProcessAlpha->SetVerboseLevel(OpVerbLevel);

  // scintillation process for heavy nuclei
  G4Scintillation* theScintProcessNuc = new G4Scintillation("Scintillation");
  // theScintProcessNuc->DumpPhysicsTable();
  theScintProcessNuc->SetTrackSecondariesFirst(true);
  theScintProcessNuc->SetScintillationYieldFactor(0.2);
  theScintProcessNuc->SetScintillationExcitationRatio(1.0);
  theScintProcessNuc->SetVerboseLevel(OpVerbLevel);

  // optical processes
  G4OpAbsorption* theAbsorptionProcess = new G4OpAbsorption();
  //  G4OpRayleigh* theRayleighScatteringProcess = new G4OpRayleigh();
  G4OpBoundaryProcess* theBoundaryProcess = new G4OpBoundaryProcess();
  //  theAbsorptionProcess->DumpPhysicsTable();
  //  theRayleighScatteringProcess->DumpPhysicsTable();
  theAbsorptionProcess->SetVerboseLevel(OpVerbLevel);
  // theRayleighScatteringProcess->SetVerboseLevel(OpVerbLevel);
  theBoundaryProcess->SetVerboseLevel(OpVerbLevel);
  G4OpticalSurfaceModel themodel = unified;
  theBoundaryProcess->SetModel(themodel);

  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
      if (theScintProcessDef->IsApplicable(*particle)) {
	//      if(particle->GetPDGMass() > 5.0*GeV) 
	if(particle->GetParticleName() == "GenericIon") {
	  pmanager->AddProcess(theScintProcessNuc); // AtRestDiscrete
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessNuc,idxPostStep);
	}	  
	else if(particle->GetParticleName() == "alpha") {
	  pmanager->AddProcess(theScintProcessAlpha);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessAlpha,idxPostStep);
	}
	else {
	  pmanager->AddProcess(theScintProcessDef);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxAtRest);
	  pmanager->SetProcessOrderingToLast(theScintProcessDef,idxPostStep);
	}	  
      }
      
      if (particleName == "opticalphoton") {
	pmanager->AddDiscreteProcess(theAbsorptionProcess);
	//	pmanager->AddDiscreteProcess(theRayleighScatteringProcess);
	pmanager->AddDiscreteProcess(theBoundaryProcess);
      }
    }
}


// Hadronic processes ////////////////////////////////////////////////////////

// Elastic processes:
#include "G4HadronElasticProcess.hh"

// Inelastic processes:
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

// Low-energy Models: < 20GeV
#include "G4LElastic.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LEKaonPlusInelastic.hh"
#include "G4LEKaonZeroSInelastic.hh"
#include "G4LEKaonZeroLInelastic.hh"
#include "G4LEKaonMinusInelastic.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LEAntiNeutronInelastic.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

// High-energy Models: >20 GeV
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4HEKaonPlusInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonZeroInelastic.hh"
#include "G4HEKaonMinusInelastic.hh"
#include "G4HEProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4HENeutronInelastic.hh"
#include "G4HEAntiNeutronInelastic.hh"

// Neutron high-precision models: <20 MeV
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronHPCapture.hh"
#include "G4NeutronHPCaptureData.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
#include "G4LCapture.hh"

// Stopping processes
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4KaonMinusAbsorptionAtRest.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4AntiNeutronAnnihilationAtRest.hh"


// ConstructHad()
// Makes discrete physics processes for the hadrons, at present limited
// to those particles with GHEISHA interactions (INTRC > 0).
// The processes are: Elastic scattering and Inelastic scattering.
// F.W.Jones  09-JUL-1998
void DMXPhysicsList::ConstructHad() 
{
  G4HadronElasticProcess* theElasticProcess = new G4HadronElasticProcess;
  G4LElastic* theElasticModel = new G4LElastic;
  theElasticProcess->RegisterMe(theElasticModel);
  
  theParticleIterator->reset();
  while ((*theParticleIterator)()) 
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();

      if (particleName == "pi+") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4PionPlusInelasticProcess* theInelasticProcess = 
	    new G4PionPlusInelasticProcess("inelastic");
	  G4LEPionPlusInelastic* theLEInelasticModel = 
	    new G4LEPionPlusInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEPionPlusInelastic* theHEInelasticModel = 
	    new G4HEPionPlusInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	} 

      else if (particleName == "pi-") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4PionMinusInelasticProcess* theInelasticProcess = 
	    new G4PionMinusInelasticProcess("inelastic");
	  G4LEPionMinusInelastic* theLEInelasticModel = 
	    new G4LEPionMinusInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEPionMinusInelastic* theHEInelasticModel = 
	    new G4HEPionMinusInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	  G4String prcNam;
	  pmanager->AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);
	}
      
      else if (particleName == "kaon+") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4KaonPlusInelasticProcess* theInelasticProcess = 
	    new G4KaonPlusInelasticProcess("inelastic");
	  G4LEKaonPlusInelastic* theLEInelasticModel = 
	    new G4LEKaonPlusInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEKaonPlusInelastic* theHEInelasticModel = 
	    new G4HEKaonPlusInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}
      
      else if (particleName == "kaon0S") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4KaonZeroSInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroSInelasticProcess("inelastic");
	  G4LEKaonZeroSInelastic* theLEInelasticModel = 
	    new G4LEKaonZeroSInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEKaonZeroInelastic* theHEInelasticModel = 
	    new G4HEKaonZeroInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "kaon0L") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4KaonZeroLInelasticProcess* theInelasticProcess = 
	    new G4KaonZeroLInelasticProcess("inelastic");
	  G4LEKaonZeroLInelastic* theLEInelasticModel = 
	    new G4LEKaonZeroLInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEKaonZeroInelastic* theHEInelasticModel = 
	    new G4HEKaonZeroInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "kaon-") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4KaonMinusInelasticProcess* theInelasticProcess = 
	    new G4KaonMinusInelasticProcess("inelastic");
	  G4LEKaonMinusInelastic* theLEInelasticModel = 
	    new G4LEKaonMinusInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEKaonMinusInelastic* theHEInelasticModel = 
	    new G4HEKaonMinusInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	  pmanager->AddRestProcess(new G4KaonMinusAbsorptionAtRest, ordDefault);
	}

      else if (particleName == "proton") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4ProtonInelasticProcess* theInelasticProcess = 
	    new G4ProtonInelasticProcess("inelastic");
	  G4LEProtonInelastic* theLEInelasticModel = new G4LEProtonInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEProtonInelastic* theHEInelasticModel = new G4HEProtonInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "anti_proton") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4AntiProtonInelasticProcess* theInelasticProcess = 
	    new G4AntiProtonInelasticProcess("inelastic");
	  G4LEAntiProtonInelastic* theLEInelasticModel = 
	    new G4LEAntiProtonInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEAntiProtonInelastic* theHEInelasticModel = 
	    new G4HEAntiProtonInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "neutron") {
	// elastic scattering
	G4HadronElasticProcess* theNeutronElasticProcess = 
	  new G4HadronElasticProcess;
	G4LElastic* theElasticModel1 = new G4LElastic;
	G4NeutronHPElastic * theElasticNeutron = new G4NeutronHPElastic;
	theNeutronElasticProcess->RegisterMe(theElasticModel1);
	theElasticModel1->SetMinEnergy(19*MeV);
	theNeutronElasticProcess->RegisterMe(theElasticNeutron);
	G4CrossSectionDataStore * theStore = 
	  ((G4HadronElasticProcess*)theNeutronElasticProcess)
	  ->GetCrossSectionDataStore();
	G4NeutronHPElasticData * theNeutronData = new G4NeutronHPElasticData;
	theStore->AddDataSet(theNeutronData);
	pmanager->AddDiscreteProcess(theNeutronElasticProcess);
	// inelastic scattering
	G4NeutronInelasticProcess* theInelasticProcess =
	  new G4NeutronInelasticProcess("inelastic");
	G4LENeutronInelastic* theInelasticModel = new G4LENeutronInelastic;
	theInelasticModel->SetMinEnergy(19*MeV);
	theInelasticProcess->RegisterMe(theInelasticModel);
	G4NeutronHPInelastic * theLENeutronInelasticModel =
	  new G4NeutronHPInelastic;
	theInelasticProcess->RegisterMe(theLENeutronInelasticModel);
	G4CrossSectionDataStore * theStore1 = 
	  ((G4HadronInelasticProcess*)theInelasticProcess)
	  ->GetCrossSectionDataStore();
	G4NeutronHPInelasticData * theNeutronData1 = 
	  new G4NeutronHPInelasticData;
	theStore1->AddDataSet(theNeutronData1);
	pmanager->AddDiscreteProcess(theInelasticProcess);
	// capture
	G4HadronCaptureProcess* theCaptureProcess =
	  new G4HadronCaptureProcess;
	G4LCapture* theCaptureModel = new G4LCapture;
	theCaptureModel->SetMinEnergy(19*MeV);
	theCaptureProcess->RegisterMe(theCaptureModel);
	G4NeutronHPCapture * theLENeutronCaptureModel = new G4NeutronHPCapture;
	theCaptureProcess->RegisterMe(theLENeutronCaptureModel);
	G4CrossSectionDataStore * theStore3 = 
	  ((G4HadronCaptureProcess*)theCaptureProcess)->
	  GetCrossSectionDataStore();
	G4NeutronHPCaptureData * theNeutronData3 = new G4NeutronHPCaptureData;
	theStore3->AddDataSet(theNeutronData3);
	pmanager->AddDiscreteProcess(theCaptureProcess);
	//  G4ProcessManager* pmanager = G4Neutron::Neutron->GetProcessManager();
	//  pmanager->AddProcess(new G4UserSpecialCuts(),-1,-1,1);
      }
      else if (particleName == "anti_neutron") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4AntiNeutronInelasticProcess* theInelasticProcess = 
	    new G4AntiNeutronInelasticProcess("inelastic");
	  G4LEAntiNeutronInelastic* theLEInelasticModel = 
	    new G4LEAntiNeutronInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  G4HEAntiNeutronInelastic* theHEInelasticModel = 
	    new G4HEAntiNeutronInelastic;
	  theInelasticProcess->RegisterMe(theHEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "deuteron") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4DeuteronInelasticProcess* theInelasticProcess = 
	    new G4DeuteronInelasticProcess("inelastic");
	  G4LEDeuteronInelastic* theLEInelasticModel = 
	    new G4LEDeuteronInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}
      
      else if (particleName == "triton") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4TritonInelasticProcess* theInelasticProcess = 
	    new G4TritonInelasticProcess("inelastic");
	  G4LETritonInelastic* theLEInelasticModel = 
	    new G4LETritonInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

      else if (particleName == "alpha") 
	{
	  pmanager->AddDiscreteProcess(theElasticProcess);
	  G4AlphaInelasticProcess* theInelasticProcess = 
	    new G4AlphaInelasticProcess("inelastic");
	  G4LEAlphaInelastic* theLEInelasticModel = 
	    new G4LEAlphaInelastic;
	  theInelasticProcess->RegisterMe(theLEInelasticModel);
	  pmanager->AddDiscreteProcess(theInelasticProcess);
	}

    }
}


// Decays ///////////////////////////////////////////////////////////////////
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4Ions.hh"

void DMXPhysicsList::ConstructGeneral() {

  // Add Decay Process
  G4Decay* theDecayProcess = new G4Decay();
  theParticleIterator->reset();
  while( (*theParticleIterator)() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      
      if (theDecayProcess->IsApplicable(*particle)) 
	{ 
	  pmanager ->AddProcess(theDecayProcess);
	  // set ordering for PostStepDoIt and AtRestDoIt
	  pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
	  pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
	}
    }

  // Declare radioactive decay to the GenericIon in the IonTable.
  const G4IonTable *theIonTable = 
    G4ParticleTable::GetParticleTable()->GetIonTable();
  G4RadioactiveDecay *theRadioactiveDecay = new G4RadioactiveDecay();

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
}

// Cuts /////////////////////////////////////////////////////////////////////
void DMXPhysicsList::SetCuts() 
{
  
  if (verboseLevel >1)
    G4cout << "DMXPhysicsList::SetCuts:";
  
  if (verboseLevel>0){
    G4cout << "DMXPhysicsList::SetCuts:";
    G4cout << "CutLength : " 
	   << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }

  //special for low energy physics
  G4double lowlimit=250*eV;  
  G4Gamma   ::SetEnergyRange(lowlimit,100*GeV);
  G4Electron::SetEnergyRange(lowlimit,100*GeV);
  G4Positron::SetEnergyRange(lowlimit,100*GeV);
   
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

