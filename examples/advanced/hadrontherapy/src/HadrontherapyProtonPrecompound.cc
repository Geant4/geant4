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
// PRECOMPOUND + EVAPORATION(DEFAULT EVAPORATION) FOR PROTONS, NEUTRONS AND PIONS
// 
// LEP MODEL UP TO 100 MEV AND BINARY ION MODEL BETWEEN 80 MEV AND 40. GEV 
// FOR  DEUTERON, TRITON, HE3, ALPHA
// 
//FISSION AND HADRON CAPTURE FOR NEUTRONS BETWEEN 0. MEV AND 100. TEV
//
HadrontherapyProtonPrecompound::HadrontherapyProtonPrecompound(const G4String& name): 
  G4VPhysicsConstructor(name)
{

  G4cout<<"****** Proton Precompound Physics List is active !!!!!! ******"<<G4endl;
  // Inelastic hadronic process: energy limits 

  // Protons, neutrons and pions
  // Energy limits of the precompound model
  precompoundLowLimit = 0.*MeV;
  precompoundHighLimit = 300.*MeV;

  // Energy limit of the neutron fission and capture
  neutronLowLimit = 0.*TeV;
  neutronHighLimit = 100.*TeV;
 
  // Ions 

  // Energy limit of the LEP model for ions
  LEPHighLimit = 100.*MeV;
  // Energy limit of the binary ion model
  binaryLightIonLowLimit = 80.*MeV;
  binaryLightIonHighLimit = 40.*GeV;
}

HadrontherapyProtonPrecompound::~HadrontherapyProtonPrecompound()
{}

void HadrontherapyProtonPrecompound::ConstructProcess()
{
 G4ParticleDefinition* particle = 0;
 G4ProcessManager* pmanager = 0;

  //  ELASTIC SCATTERING 
  // FOR PROTON, NEUTRON, IONS 
  G4LElastic* elasticScattering_model = new G4LElastic();
  G4HadronElasticProcess* elastic_scattering = new G4HadronElasticProcess();
  elastic_scattering -> RegisterMe(elasticScattering_model);
 
  // INELASTIC SCATTERING
  // Binary Cascade
  G4BinaryLightIonReaction* ionBinaryCascade = new G4BinaryLightIonReaction();
  // Set the min and max energy for the Binary Cascade model
  ionBinaryCascade -> SetMinEnergy(binaryLightIonLowLimit);
  ionBinaryCascade -> SetMaxEnergy(binaryLightIonHighLimit);
  
  // TRIPATHI CROSS SECTION
  // Implementation of formulas in analogy to NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* TripathiCrossSection = new G4TripathiCrossSection;
  
  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

//------------------------------------------//
// Activate the hadronic physics processes  //
//------------------------------------------//

  // Proton PRECOMPOUND + EVAPORATION(DEFAULT EVAPORATION)
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  G4PreCompoundModel* preequilibriumModel = new G4PreCompoundModel(&theHandler);
   
  // Set the minimum and maximum energy value of the pre-equilibrium model
  preequilibriumModel -> SetMinEnergy(precompoundLowLimit);
  preequilibriumModel -> SetMaxEnergy(precompoundHighLimit);
 
  // Model Registration 
  protonInelasticProcess.RegisterMe(preequilibriumModel);
  // Activate the cross-sections for proton nuclear scattering 
  protonInelasticProcess.AddDataSet(&protonInelasticCrossSection);
   
  // Active the proton inelastic scattering using the precompound model
  pmanager -> AddDiscreteProcess(&protonInelasticProcess);
  // Active the Hadron Elastic Process 
  pmanager -> AddDiscreteProcess(elastic_scattering); 

 // Neutron
  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
		  
  // Register the Precompound model
  neutronInelasticProcess.RegisterMe(preequilibriumModel);
  // Active the Cross-sections for neutron nuclear scattering
  neutronInelasticProcess.AddDataSet(&neutronInelasticCrossSection);
  // Active the neutron inelastic process
  pmanager -> AddDiscreteProcess(&neutronInelasticProcess);
  // Active the neutron elastic scattering
  pmanager -> AddDiscreteProcess(elastic_scattering); 
		  
  // Pions plus
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
 
  // Define the inelastic process for pions plus
  G4PionPlusInelasticProcess* pionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");
  pionPlusInelasticProcess -> RegisterMe(preequilibriumModel);
  // Active the inelastic process for pions plus
  pmanager -> AddDiscreteProcess(pionPlusInelasticProcess);
  pmanager -> AddDiscreteProcess(elastic_scattering);

  // Pion Minus
  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();

  // Define the inelastic process for pions minus
  G4PionMinusInelasticProcess* pionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
  // Register the inelastic model for pion minus  
  pionMinusInelasticProcess -> RegisterMe(preequilibriumModel);
  // Active the inelastic process for pion minus
  pmanager -> AddDiscreteProcess(pionMinusInelasticProcess); 
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // Final state production model for Deuteron inelastic scattering below 100 MeV
  G4LEDeuteronInelastic* deuteronLEModel = new G4LEDeuteronInelastic;
  // Set the maximum energy for LEP model
  deuteronLEModel -> SetMaxEnergy(LEPHighLimit);

  // Active the Tripathi and aShen Cross Section
  deuteronInelasticProcess.AddDataSet(TripathiCrossSection);
  deuteronInelasticProcess.AddDataSet(aShen);

  // Register the deuteron inelastic scattering models
  deuteronInelasticProcess.RegisterMe(deuteronLEModel);
  deuteronInelasticProcess.RegisterMe(ionBinaryCascade);
  
  // Active the deuteron inelastic scattering using the deuteron inelastic and binary cascade model
  pmanager -> AddDiscreteProcess(&deuteronInelasticProcess);
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // triton
  particle = G4Triton::Triton();
  pmanager = particle->GetProcessManager();
  
  // Final state production model for Triton inelastic scattering below 100 MeV
  G4LETritonInelastic* tritonLEModel = new G4LETritonInelastic;
  // Set the maximum energy for LEP model
  tritonLEModel -> SetMaxEnergy(LEPHighLimit);

  // Active the Tripathi and aShen Cross Section
  tritonInelasticProcess.AddDataSet(TripathiCrossSection);
  tritonInelasticProcess.AddDataSet(aShen);

  // Register the triton inelastic scattering models
  tritonInelasticProcess.RegisterMe(tritonLEModel);
  tritonInelasticProcess.RegisterMe(ionBinaryCascade);

  // Active the triton inelastic scattering process
  pmanager -> AddDiscreteProcess(&tritonInelasticProcess);
  // Active the triton elastic scattering
  pmanager -> AddDiscreteProcess(elastic_scattering);

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
		  
  // Final state production model for Alpha inelastic scattering below 20 GeV
  G4LEAlphaInelastic* alphaLEModel = new G4LEAlphaInelastic;
  // Set the maximum energy for LEP model
  alphaLEModel -> SetMaxEnergy(LEPHighLimit);
 
 // Register the alpha inelastic scattering models
  alphaInelasticProcess.AddDataSet(TripathiCrossSection);
  alphaInelasticProcess.AddDataSet(aShen);
 
  // Register the Alpha Inelastic and Binary Cascade Model
  alphaInelasticProcess.RegisterMe(alphaLEModel);
  alphaInelasticProcess.RegisterMe(ionBinaryCascade);
 
  // Active the alpha inelastic scattering using the alpha inelastic and binary cascade model
  pmanager -> AddDiscreteProcess(&alphaInelasticProcess);
  // Active the alpha elastic scattering
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // He3
  particle = G4He3::He3();
  pmanager = particle -> GetProcessManager();

   // Inelastic Scattering for ions
  G4HadronInelasticProcess* He3InelasticProcess = new G4HadronInelasticProcess("He3Inelastic",particle);

   G4BinaryLightIonReaction * theGenericIonBC= new G4BinaryLightIonReaction;
   // Active the Tripathi and aShen Cross Section
   He3InelasticProcess -> AddDataSet(TripathiCrossSection);
   He3InelasticProcess -> AddDataSet(aShen);

   // Register the He3 Binary Cascade Model
   He3InelasticProcess -> RegisterMe(theGenericIonBC);
   // Active the Inelastic Process for He3
   pmanager -> AddDiscreteProcess(He3InelasticProcess);
   // Active the Hadron Elastic Process
   pmanager -> AddDiscreteProcess(elastic_scattering); 
  
  //HADRON CAPTURE
  // Process for capture of neutral hadrons
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  // Final state production model for capture of neutral hadrons in nuclei
  G4LCapture* capture_model = new G4LCapture();
  // Set the energy range for the capture model
  capture_model -> SetMinEnergy(neutronLowLimit);
  capture_model -> SetMaxEnergy(neutronHighLimit);
  // Register the capture model
  neutronCapture -> RegisterMe(capture_model);
  // Active the neutron capture process
  pmanager -> AddDiscreteProcess(neutronCapture);

  //FISSION
  // Process for induced fission
  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  //Final state production model for induced fission
  G4LFission* fission_model = new G4LFission();
  // Set the energy range for the fission model
  fission_model -> SetMinEnergy(neutronLowLimit);
  fission_model -> SetMaxEnergy(neutronHighLimit);
  // Register the fission model
  fission -> RegisterMe(fission_model); 
  // Active the fission process
  pmanager -> AddDiscreteProcess(fission);  
}



