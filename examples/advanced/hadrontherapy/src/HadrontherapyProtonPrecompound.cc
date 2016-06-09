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
// $Id: HadrontherapyProtonPrecompound.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the INFN, Catania, Italy
// (b) INFN Section of Genova, Genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyProtonPrecompound.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LElastic.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4PiNuclearCrossSection.hh" 
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4LETritonInelastic.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4HadronCaptureProcess.hh"
//
// PRECOMPOUND PHYSICS LIST
//
// PRECOMPOUND + EVAPORATION(DEFAULT EVAPORATION) NO FERMI BREAK-UP FOR PROTONS, NEUTRONS AND PIONS
// 
// LEP MODEL UP TO 100 MEV AND BINARY ION MODEL BETWEEN 80 MEV AND 40. GEV 
// FOR  DEUTERON, TRITON, HE3, ALPHA
// 
//FISSION AND HADRON CAPTURE FOR NEUTRONS BETWEEN 0. MEV AND 100. TEV
//
HadrontherapyProtonPrecompound::HadrontherapyProtonPrecompound(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  // Inelastic process, energy limits 

  // Protons, neutrons and pions
  // Energy limit of the precompound model
  precompoundLowLimit = 0.*MeV;
  precompoundHighLimit = 300.*MeV;

  // Energy limit of the neutron fission and capture
  neutronLowLimit = 0.*TeV;
  neutronHighLimit = 100.*TeV;
 
  // Ions
  // Energy limit of the binary ion model
  binaryLightIonLowLimit = 80.*MeV;
  binaryLightIonHighLimit = 40.*GeV;

  // Energy limit of the LEP model for ions
  LEPHighLimit = 100.*MeV;
}

HadrontherapyProtonPrecompound::~HadrontherapyProtonPrecompound()
{}

void HadrontherapyProtonPrecompound::ConstructProcess()
{

  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
 
  // Elastic scattering
  // for protons, neutrons, pions, ions
  G4LElastic* elastic_Model = new G4LElastic();
  G4HadronElasticProcess* elastic = new G4HadronElasticProcess();
  elastic -> RegisterMe(elastic_Model);
 
  // Activation of Precompound model + evaporation (default evaporation)
  // The Fermi break-up is not active (default condition) 
  // for protons, pions and neutrons
 
  G4PreCompoundModel* preEquilibrium = new G4PreCompoundModel(&theHandler);
  preEquilibrium -> SetMinEnergy(precompoundLowLimit);
  preEquilibrium -> SetMaxEnergy(precompoundHighLimit);

  //----------//
  // Protons  //
  //----------//

  // Inelastic process for protons
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  // Activation of the Precompound model for the proton 
  // inelastic process
  //G4ProtonInelasticProcess protonInelasticScattering;
  protonInelasticScattering.RegisterMe(preEquilibrium);  
  
  // Set the cross section for the proton inelastic process
  //G4ProtonInelasticCrossSection  protonInelasticCrossSection;
  protonInelasticScattering.AddDataSet(&protonInelasticCrossSection);

  // Set the proton hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&protonInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 
 

  //----------//
  // Neutrons //
  //----------//

  // Inelastic process for neutrons
  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
  
  // Activation of the Precompound model for the neutrons
  // inelastic process
  
  neutronInelasticScattering.RegisterMe(preEquilibrium);

  // Set the cross section for the neutron inelastic process
  
  neutronInelasticScattering.AddDataSet(&neutronInelasticCrossSection);

  // Set the neutron hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&neutronInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 

  //----------//
  // Pions    //
  //----------//

  // Inelastic process for pions
 
  // Pions plus
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
 
  // Activation of the Precompound model for the pion plus
  // inelastic process
  
  pionPlusInelasticScattering.RegisterMe(preEquilibrium);

  // Set the cross section for the pion inelastic process
  
  pionPlusInelasticScattering.AddDataSet(&pionPlusInelasticCrossSection);

  // Set the pion plus hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&pionPlusInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 

  // Pions minus
  particle = G4PionMinus::PionMinus(); 
  pmanager = particle -> GetProcessManager();
 
  // Activation of the Precompound model for the pion minus
  // inelastic process
 
  pionMinusInelasticScattering.RegisterMe(preEquilibrium);

  // Set the cross section for the pion minus inelastic process
  
  pionMinusInelasticScattering.AddDataSet(&pionMinusInelasticCrossSection);

  // Set the pion minus hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&pionMinusInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 
  pmanager -> AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault); 

  //------------------//
  //d, t, 3He, alpha //
  //------------------//

  //Inelastic scattering for d, t, 3He, alpha
  
  //Activation of the Low Energy Parameterised (LEP) model and binary ion model
  // for the inelastic scattering 

  G4BinaryLightIonReaction* binaryIonCascade = new G4BinaryLightIonReaction();
  binaryIonCascade -> SetMinEnergy(binaryLightIonLowLimit);
  binaryIonCascade -> SetMaxEnergy(binaryLightIonHighLimit);

  // Set the cross section for the inelastic scattering
  TripathiCrossSection = new G4TripathiCrossSection;
  ShenCrossSection = new G4IonsShenCrossSection;

  // Deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // Activate the LEP model for the d inelastic scattering 
  G4LEDeuteronInelastic* deuteronLEPModel = new G4LEDeuteronInelastic;
  deuteronLEPModel -> SetMaxEnergy(LEPHighLimit);

  // Activate the binary ion model model for the d inelastic scattering 
 
  deuteronInelasticScattering.RegisterMe(deuteronLEPModel);
  deuteronInelasticScattering.RegisterMe(binaryIonCascade);

  // Definition of the the deuteron inelastic scattering cross section 
  deuteronInelasticScattering.AddDataSet(TripathiCrossSection);
  deuteronInelasticScattering.AddDataSet(ShenCrossSection);

  // Set the deuteron hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&deuteronInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 

  // Triton
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();

  // Activate the LEP model for the t inelastic scattering 
  G4LETritonInelastic* tritonLEPModel = new G4LETritonInelastic;
  tritonLEPModel -> SetMaxEnergy(LEPHighLimit);

  // Activate the binary ion model model for the t inelastic scattering 
 
  tritonInelasticScattering.RegisterMe(tritonLEPModel);
  tritonInelasticScattering.RegisterMe(binaryIonCascade);

  // Definition of the the triton inelastic scattering cross section 
  tritonInelasticScattering.AddDataSet(TripathiCrossSection);
  tritonInelasticScattering.AddDataSet(ShenCrossSection);

  // Set the triton hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&tritonInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();

  // Activate the LEP model for the alpha inelastic scattering 
  G4LEAlphaInelastic* alphaLEPModel = new G4LEAlphaInelastic;
  alphaLEPModel -> SetMaxEnergy(LEPHighLimit);

  // Activate the binary ion model model for the alpha inelastic scattering 
  alphaInelasticScattering.RegisterMe(alphaLEPModel);
  alphaInelasticScattering.RegisterMe(binaryIonCascade);

  // Definition of the alpha inelastic scattering cross section 
  alphaInelasticScattering.AddDataSet(TripathiCrossSection);
  alphaInelasticScattering.AddDataSet(ShenCrossSection);

  // Set the alpha hadronic processes: inelastic and elastic scattering
  pmanager -> AddDiscreteProcess(&alphaInelasticScattering);
  pmanager -> AddDiscreteProcess(elastic); 



/////////////////////////////////////////////////////////////////////////////
// NEW HE3
/////////////////////////////////////////////////////////////////////////////

  //TripathiCrossSection = new G4TripathiCrossSection;
  //ShenCrossSection = new G4IonsShenCrossSection;

  pmanager->AddDiscreteProcess(elastic);
  He3InelasticProcess = new G4HadronInelasticProcess
    ("He3Inelastic", G4He3::He3());
  
  He3InelasticProcess->AddDataSet(TripathiCrossSection);
  He3InelasticProcess->AddDataSet(ShenCrossSection);
  He3InelasticProcess->RegisterMe(binaryIonCascade);
  pmanager->AddDiscreteProcess(He3InelasticProcess);


 //  // He3
//   particle = G4He3::He3();
//   pmanager = particle -> GetProcessManager();

//   // Only Binary Ion model activated
//   // LEP model is not availaboe for 3He
  
//   G4HadronInelasticProcess*  He3InelasticScattering = 
//     new G4HadronInelasticProcess("Inelastic scattering", particle);

//   // Change the binary ion model limits because the LEP is not available
//   // 
//   binaryIonCascade -> SetMinEnergy(0.*MeV);
//   binaryIonCascade -> SetMaxEnergy(binaryLightIonHighLimit);

//   // Definition of the 3He inelastic scattering cross section 
//   He3InelasticScattering -> AddDataSet(TripathiCrossSection);
//   He3InelasticScattering -> AddDataSet(aShen);
 
//   He3InelasticScattering -> RegisterMe(binaryIonCascade);

//   // Set the He3 hadronic processes: inelastic and elastic scattering
//   pmanager -> AddDiscreteProcess(He3InelasticScattering);
//   pmanager -> AddDiscreteProcess(elastic); 
  
  // Other neutron processes
 
  //Hadron Capture 
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  G4LCapture* capture_model = new G4LCapture();
  capture_model -> SetMinEnergy(neutronLowLimit);
  capture_model -> SetMaxEnergy(neutronHighLimit);
  neutronCapture -> RegisterMe(capture_model);
  pmanager -> AddDiscreteProcess(neutronCapture);

  //Fission
  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  G4LFission* fission_model = new G4LFission();
  fission_model -> SetMinEnergy(neutronLowLimit);
  fission_model -> SetMaxEnergy(neutronHighLimit);
  fission -> RegisterMe(fission_model); 
  pmanager -> AddDiscreteProcess(fission);     
  
}



