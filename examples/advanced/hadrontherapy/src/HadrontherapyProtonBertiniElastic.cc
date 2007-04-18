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
// $Id: HadrontherapyProtonBertiniElastic.cc; May 2005
// ----------------------------------------------------------------------------
//                 GEANT 4 - Hadrontherapy example
// ----------------------------------------------------------------------------
// Code developed by:
//
// G.A.P. Cirrone(a)*, F. Di Rosa(a), S. Guatelli(b), G. Russo(a)
// 
// (a) Laboratori Nazionali del Sud 
//     of the National Institute for Nuclear Physics, Catania, Italy
// (b) National Institute for Nuclear Physics Section of Genova, genova, Italy
// 
// * cirrone@lns.infn.it
// ----------------------------------------------------------------------------

#include "HadrontherapyProtonBertiniElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4LElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"
#include "G4CascadeElasticInterface.hh"

//
// BERTINI PHYSICS LIST
//
// BERTINI FOR PROTONS, NEUTRONS AND PIONS
// 
// LEP MODEL UP TO 100 MEV AND BINARY ION MODEL BETWEEN 80 MEV AND 40. GEV 
// FOR  DEUTERON, TRITON, HE3, ALPHA
// 
// FISSION AND HADRON CAPTURE FOR NEUTRONS BETWEEN 0. MEV AND 100. TEV
//
HadrontherapyProtonBertiniElastic::HadrontherapyProtonBertiniElastic(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout << "The Bertini model (for elastic and inelastic scattering) is set for protons, neutrons and pions !!!!" << G4endl;

  // Inelastic process, energy limits 

  //
  // The Bertini model is set for protons, neutrons and pions
  // This model contains a pre-equilibrium model and a de-excitation model
  // Energy limit of the Bertini model
  bertiniLowEnergyLimit = 0.* MeV;
  bertiniHighEnergyLimit = 300.*MeV;

  // Energy limit of the neutron fission and capture
  neutronLowEnergyLimit = 0.* TeV;
  neutronHighEnergyLimit = 100.* TeV;
 
  // Ions:
  // The inelastic scattering is modelled with  LEP model up to 100 MeV,
  // then Binary Ion Model 
  // Energy limit of the LEP model for ions
  LEPHighEnergyLimit = 100.*MeV;

  // Energy limit of the binary ion model
  binaryLightIonLowEnergyLimit = 80.* MeV;
  binaryLightIonHighEnergyLimit = 40.*GeV;
}

HadrontherapyProtonBertiniElastic::~HadrontherapyProtonBertiniElastic()
{}

void HadrontherapyProtonBertiniElastic::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

  // BERTINI ELASTIC SCATTERING 
  // FOR PROTON, NEUTRON, PIONS
 
  // eliminate "the"
  G4CascadeElasticInterface* theBertiniElasticModel = new G4CascadeElasticInterface;
  // G4HadronElasticProcess* hadronElasticProcess = new G4HadronElasticProcess();
  G4HadronElasticProcess* bertiniElasticModel = new G4HadronElasticProcess();
  bertiniElasticModel -> RegisterMe(theBertiniElasticModel);

  // LOW ENERGY ELASTIC SCATTERING 
  // FOR  IONS
  G4LElastic* elasticLEmodel = new G4LElastic();
  // G4HadronElasticProcess* ionElasticProcess = new G4HadronElasticProcess();
  G4HadronElasticProcess* elasticScattering = new G4HadronElasticProcess();
  elasticScattering -> RegisterMe(elasticLEmodel);
  
  // INELASTIC SCATTERING
  // Bertini Model for protons, pions and neutrons
  G4CascadeInterface* theBertiniModel = new G4CascadeInterface;
  
  // Set the min and max energy for the Bertini Model
  theBertiniModel -> SetMinEnergy(bertiniLowEnergyLimit);
  theBertiniModel -> SetMaxEnergy(bertiniHighEnergyLimit);

  // Binary Cascade for deuteron, triton, alpha particle, He3
  G4BinaryLightIonReaction* theBinaryCascade = new G4BinaryLightIonReaction();
  // Set the min and max energy for the Binary Cascade
  theBinaryCascade -> SetMinEnergy(binaryLightIonLowEnergyLimit);
  theBinaryCascade -> SetMaxEnergy(binaryLightIonHighEnergyLimit);

  // TRIPATHI CROSS SECTION
  // Implementation of formulas in analogy to NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* tripathiCrossSection = new G4TripathiCrossSection;
  
  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

  //--------------------------------------------------------------------------------------
  // Proton BERTINI MODEL
  particle = G4Proton::Proton();
  processManager = particle -> GetProcessManager();
  
  // Model Registration 
  theProtonInelasticProcess.RegisterMe(theBertiniModel);
  // Active the Cross-sections for proton nuclear scattering up to 20 GeV
  theProtonInelasticProcess.AddDataSet(&theProtonCrossSection);

  // Active the proton inelastic scattering 
  processManager -> AddDiscreteProcess(&theProtonInelasticProcess);
  // Active the Hadron Elastic Process 
  processManager -> AddDiscreteProcess(bertiniElasticModel); 

  // Deuteron 
  particle = G4Deuteron::Deuteron();
  processManager = particle -> GetProcessManager();

  // Final state production model for deuteron inelastic scattering below 100 MeV
  G4LEDeuteronInelastic* theDeuteronLEInelasticModel = new G4LEDeuteronInelastic;
  // Set the maximum energy for LEP model
  theDeuteronLEInelasticModel -> SetMaxEnergy(LEPHighEnergyLimit);

  // Active the Tripathi and aShen Cross Section
  theDeuteronInelasticProcess.AddDataSet(tripathiCrossSection);
  theDeuteronInelasticProcess.AddDataSet(aShen);

  // Register the Parameterised Deuteron Inelastic Model and the Ion Binary Cascade Model
  theDeuteronInelasticProcess.RegisterMe(theDeuteronLEInelasticModel);
  theDeuteronInelasticProcess.RegisterMe(theBinaryCascade);

  // Active the deuteron elastic and inelastic scattering 
  processManager -> AddDiscreteProcess(&theDeuteronInelasticProcess);
  // Active the Hadron Elastic Process
  processManager -> AddDiscreteProcess(elasticScattering); 

 // triton
  particle = G4Triton::Triton();
  processManager = particle -> GetProcessManager();
  
  // Final state production model for Triton inelastic scattering below 100 MeV
  G4LETritonInelastic* theTritonLEInelasticModel = new G4LETritonInelastic;
  
  // Set the maximum energy for LEP model
  theTritonLEInelasticModel -> SetMaxEnergy(LEPHighEnergyLimit);
  
  // Active the Tripathi and aShen Cross Section
  theTritonInelasticProcess.AddDataSet(tripathiCrossSection);
  theTritonInelasticProcess.AddDataSet(aShen);
  
  // Register the Triton Inelastic and Binary Cascade Model
   theTritonInelasticProcess.RegisterMe(theTritonLEInelasticModel);
   theTritonInelasticProcess.RegisterMe(theBinaryCascade);
  
  // Active the triton inelastic scattering using the triton inelastic and binary cascade model
  processManager -> AddDiscreteProcess(&theTritonInelasticProcess);
  
  // Active the Hadron Elastic Process
  processManager -> AddDiscreteProcess(elasticScattering);

  // Alpha
  particle = G4Alpha::Alpha();
  processManager = particle -> GetProcessManager();
  // Final state production model for Alpha inelastic scattering below 20 GeV
  G4LEAlphaInelastic* theAlphaLEInelasticModel = new G4LEAlphaInelastic;
  // Set the maximum energy for LEP model
  theAlphaLEInelasticModel -> SetMaxEnergy(LEPHighEnergyLimit);
  // Register the Triton Inelastic and Binary Cascade Model
  theAlphaInelasticProcess.AddDataSet(tripathiCrossSection);
  theAlphaInelasticProcess.AddDataSet(aShen);
  // Register the Alpha Inelastic and Binary Cascade Model
  theAlphaInelasticProcess.RegisterMe(theAlphaLEInelasticModel);
  theAlphaInelasticProcess.RegisterMe(theBinaryCascade);
  // Active the alpha inelastic scattering using the alpha inelastic and binary cascade model
  processManager -> AddDiscreteProcess(&theAlphaInelasticProcess);
  // Active the Hadron Elastic Process
  processManager -> AddDiscreteProcess(elasticScattering); 

// Neutron processes
  particle = G4Neutron::Neutron();
  processManager = particle -> GetProcessManager();
  // Register the Precompound model
  theNeutronInelasticProcess.RegisterMe(theBertiniModel);
  // Active the Cross-sections for neutron nuclear scattering from 14 MeV up to 20 GeV
  theNeutronInelasticProcess.AddDataSet(&theNeutronCrossSection);
  // Active the neutron inelastic process
  processManager -> AddDiscreteProcess(&theNeutronInelasticProcess);
  // Active the Hadron Elastic Process
  processManager -> AddDiscreteProcess(bertiniElasticModel);

   //HADRON CAPTURE
  // Process for capture of neutral hadrons
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  // Final state production model for capture of neutral hadrons in nuclei
  G4LCapture* captureModel = new G4LCapture();
  // Set the energy range for the capture model
  captureModel -> SetMinEnergy(neutronLowEnergyLimit);
  captureModel -> SetMaxEnergy(neutronHighEnergyLimit);
  // Register the capture model
  neutronCapture -> RegisterMe(captureModel);
  // Active the neutron capture process
  processManager -> AddDiscreteProcess(neutronCapture);

  //FISSION
  // Process for induced fission
  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  //Final state production model for induced fission
  G4LFission* fissionModel = new G4LFission();
  // Set the energy range for the fission model
  fissionModel -> SetMinEnergy(neutronLowEnergyLimit);
  fissionModel -> SetMaxEnergy(neutronHighEnergyLimit);
  // Register the fission model
  fission -> RegisterMe(fissionModel); 
  // Active the fission process
  processManager -> AddDiscreteProcess(fission); 

// Pions plus processes
  particle = G4PionPlus::PionPlus(); 
  processManager = particle -> GetProcessManager();

  // Define the inelastic process for pions plus
  G4PionPlusInelasticProcess* thePionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");

  // Register the Low Energy Inelastic Model for pions plus
  thePionPlusInelasticProcess -> RegisterMe(theBertiniModel);

  // Active the inelastic process for pions plus
  processManager -> AddDiscreteProcess(thePionPlusInelasticProcess);
  processManager -> AddDiscreteProcess(bertiniElasticModel);

  // Pion Minus processes
  particle = G4PionMinus::PionMinus();
  processManager = particle -> GetProcessManager();
  // Define the inelastic process for pions minus
  G4PionMinusInelasticProcess* thePionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
  // Register the inelastic model for pion minus  
  thePionMinusInelasticProcess -> RegisterMe(theBertiniModel);
  // Active the inelastic process for pion minus
  processManager -> AddDiscreteProcess(thePionMinusInelasticProcess); 
  // Active Absorption process for pion minus
  processManager -> AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);
  processManager -> AddDiscreteProcess(bertiniElasticModel);   
}



