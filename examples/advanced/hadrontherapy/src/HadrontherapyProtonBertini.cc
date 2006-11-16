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
// $Id: HadrontherapyProtonBertini.cc; May 2005
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

// Code review by M.G. Pia, 2 November 2006
// Further code review is needed
// ----------------------------------------------------------------------------

#include "HadrontherapyProtonBertini.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleTypes.hh"
#include "G4HadronElasticProcess.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4ExcitationHandler.hh" 
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LElastic.hh"
#include "G4CascadeInterface.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4CascadeElasticInterface.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"

// BERTINI PHYSICS LIST
//
// BERTINI FOR PROTONS, NEUTRONS AND PIONS
// 
// LEP MODEL UP TO 100 MEV AND BINARY ION MODEL BETWEEN 80 MEV AND 40. GEV 
// FOR  DEUTERON, TRITON, ALPHA
// 
// FISSION AND HADRON CAPTURE FOR NEUTRONS BETWEEN 0. MEV AND 100. TEV
//
HadrontherapyProtonBertini::HadrontherapyProtonBertini(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout << "The Bertini model (for inelastic scattering) is set for protons, neutrons and pions" << G4endl;

  // The Bertini model is set for protons, neutrons and pions
  // This model contains a pre-equilibrium model and a de-excitation model
 
  // Ions:
  // The inelastic scattering is modelled with  LEP model up to 100 MeV,
  // then Binary Ion Model 
}

HadrontherapyProtonBertini::~HadrontherapyProtonBertini()
{}

void HadrontherapyProtonBertini::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

  // Physics for proton, neutron, pion+ and pion-

  // Elastic scattering: Low Energy Parameterised model 
  G4LElastic* elasticModel = new G4LElastic();
  G4HadronElasticProcess* elasticScattering = new G4HadronElasticProcess();
  elasticScattering->RegisterMe(elasticModel);

  // Inelastic scattering: Bertini Inelastic model
  G4CascadeInterface* theBertiniModel = new G4CascadeInterface;
  
  // Energy limit of the Bertini model
  G4double bertiniLowEnergyLimit = 0.* MeV;
  G4double bertiniHighEnergyLimit = 300.*MeV;
 
  theBertiniModel->SetMinEnergy(bertiniLowEnergyLimit);
  theBertiniModel->SetMaxEnergy(bertiniHighEnergyLimit);

  //--------------------------------------------------------------------------------------
  // Proton processes
  particle = G4Proton::Proton();
  processManager = particle->GetProcessManager();
  
  // Model Registration 
  G4ProtonInelasticProcess* theProtonInelasticProcess = new G4ProtonInelasticProcess();
  theProtonInelasticProcess->RegisterMe(theBertiniModel);

  // Activate the cross-sections for proton nuclear scattering up to 20 GeV
  theProtonInelasticProcess->AddDataSet(&theProtonCrossSection);
 
  // Activate the proton inelastic scattering 
  processManager->AddDiscreteProcess(theProtonInelasticProcess);
  // Activate the elastic scattering 
  processManager->AddDiscreteProcess(elasticScattering); 

  //--------------------------------------------------------------------------------------
  // Pions plus processes
  particle = G4PionPlus::PionPlus(); 
  processManager = particle->GetProcessManager();

  // Define the inelastic process for pions plus
  G4PionPlusInelasticProcess* thePionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");
  // Register the Low Energy Inelastic Model for pions plus
  thePionPlusInelasticProcess->RegisterMe(theBertiniModel);
  // Activate the inelastic process for pions plus
  processManager->AddDiscreteProcess(thePionPlusInelasticProcess);
  // Activate the elastic process for pions plus
  processManager->AddDiscreteProcess(elasticScattering);

  //--------------------------------------------------------------------------------------
  // Pion Minus processes
  particle = G4PionMinus::PionMinus();
  processManager = particle->GetProcessManager();

  // Define the inelastic process for pions minus
  G4PionMinusInelasticProcess* thePionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
  // Register the inelastic model for pion minus  
  thePionMinusInelasticProcess->RegisterMe(theBertiniModel);
  // Activate the inelastic process for pion minus
  processManager->AddDiscreteProcess(thePionMinusInelasticProcess); 
  // Activate the elastic process for pion minus
  processManager->AddDiscreteProcess(elasticScattering);   

  //--------------------------------------------------------------------------------------
  // Neutron processes
  particle = G4Neutron::Neutron();
  processManager = particle->GetProcessManager();

  // Register the Bertini model
  G4NeutronInelasticProcess* theNeutronInelasticProcess = new G4NeutronInelasticProcess();
  theNeutronInelasticProcess->RegisterMe(theBertiniModel);

  // Activate the Cross-sections for neutron nuclear scattering from 14 MeV up to 20 GeV
  theNeutronInelasticProcess->AddDataSet(&theNeutronCrossSection);
  // Activate the neutron inelastic process
  processManager->AddDiscreteProcess(theNeutronInelasticProcess);
  // Activate the Hadron Elastic Process
  processManager->AddDiscreteProcess(elasticScattering);

  // Neutron capture process
 
  // Energy limits 
  G4double neutronLowEnergyLimit = 0. * MeV;
  G4double neutronHighEnergyLimit = 100. * TeV;

  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();
  // Final state production model for capture of neutral hadrons in nuclei
  G4LCapture* captureModel = new G4LCapture();
  // Set the energy range for the capture model
  captureModel->SetMinEnergy(neutronLowEnergyLimit);
  captureModel->SetMaxEnergy(neutronHighEnergyLimit);
  // Register the neutron capture model
  neutronCapture->RegisterMe(captureModel);
  // Activate the neutron capture process
  processManager->AddDiscreteProcess(neutronCapture);

  // Process for induced fission
  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  //Final state production model for induced fission
  G4LFission* fissionModel = new G4LFission();
  // Set the energy range for the fission model
  fissionModel->SetMinEnergy(neutronLowEnergyLimit);
  fissionModel->SetMaxEnergy(neutronHighEnergyLimit);
  // Register the fission model
  fission->RegisterMe(fissionModel); 
  // Activate the fission process
  processManager->AddDiscreteProcess(fission); 

  //--------------------------------------------------------------------------------------
  // Physics for ions

  // Energy limit of the LEP model for ions
  G4double LEPHighEnergyLimit = 100.* MeV;

  // Energy limit of the binary ion model
  G4double binaryLightIonLowEnergyLimit = 80.* MeV;
  G4double binaryLightIonHighEnergyLimit = 40.* GeV;

  // Cross section data sets

  // TRIPATHI CROSS SECTION
  // Implementation of formulas taken from NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* tripathiCrossSection = new G4TripathiCrossSection;
  
  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

  // Intra-nuclear transport: Binary Cascade Model
  // Binary Cascade for deuteron, triton, alpha particle
  G4BinaryLightIonReaction* theBinaryCascade = new G4BinaryLightIonReaction();
  // Set the min and max energy for the Binary Cascade
  theBinaryCascade->SetMinEnergy(binaryLightIonLowEnergyLimit);
  theBinaryCascade->SetMaxEnergy(binaryLightIonHighEnergyLimit);

  //--------------------------------------------------------------------------------------
  // Deuteron 
  particle = G4Deuteron::Deuteron();
  processManager = particle->GetProcessManager();

  // Final state production model for deuteron inelastic scattering below 100 MeV: Low Energy Parameterised model
  G4LEDeuteronInelastic* theDeuteronLEInelasticModel = new G4LEDeuteronInelastic;
  // Set the maximum energy for LEP model
  theDeuteronLEInelasticModel->SetMaxEnergy(LEPHighEnergyLimit);
 // G4DeuteronInelasticProcess theDeuteronInelasticProcess;
 
  // Activate the Tripathi and Shen Cross Section
  theDeuteronInelasticProcess.AddDataSet(tripathiCrossSection);
  theDeuteronInelasticProcess.AddDataSet(aShen);

  // Register the Parameterised Deuteron Inelastic Model and the Ion Binary Cascade Model
  theDeuteronInelasticProcess.RegisterMe(theDeuteronLEInelasticModel);
  theDeuteronInelasticProcess.RegisterMe(theBinaryCascade);

  // Activate the deuteron elastic and inelastic scattering 
  processManager->AddDiscreteProcess(&theDeuteronInelasticProcess);
  // Activate the Hadron Elastic Process
  processManager->AddDiscreteProcess(elasticScattering); 

  //--------------------------------------------------------------------------------------
  // Triton
  particle = G4Triton::Triton();
  processManager = particle->GetProcessManager();
  
  // Final state production model for Triton inelastic scattering below 100 MeV: Low Energy Parameterised model
  G4LETritonInelastic* theTritonLEInelasticModel = new G4LETritonInelastic;  
  // Set the maximum energy for LEP model
  theTritonLEInelasticModel->SetMaxEnergy(LEPHighEnergyLimit);

  // Activate the Tripathi and Shen Cross Section
  //G4TritonInelasticProcess theTritonInelasticProcess;
  theTritonInelasticProcess.AddDataSet(tripathiCrossSection);
  theTritonInelasticProcess.AddDataSet(aShen);

  // Register the Triton Inelastic and Binary Cascade Models
  theTritonInelasticProcess.RegisterMe(theTritonLEInelasticModel);
  theTritonInelasticProcess.RegisterMe(theBinaryCascade);

  // Activate the triton inelastic scattering using the parameterised Triton Inelastic and Binary Cascade models
  processManager->AddDiscreteProcess(&theTritonInelasticProcess);
  // Activate the Hadron Elastic Process
  processManager->AddDiscreteProcess(elasticScattering);

  //--------------------------------------------------------------------------------------
  // Alpha
  particle = G4Alpha::Alpha();
  processManager = particle->GetProcessManager();

  // Final state production model for Alpha inelastic scattering below 20 GeV: Low Energy Parameterised model
  G4LEAlphaInelastic* theAlphaLEInelasticModel = new G4LEAlphaInelastic;
  // Set the maximum energy for LEP model
  theAlphaLEInelasticModel->SetMaxEnergy(LEPHighEnergyLimit);

  //G4AlphaInelasticProcess theAlphaInelasticProcess;

  //  Activate the Tripathi and Shen Cross Section
  theAlphaInelasticProcess.AddDataSet(tripathiCrossSection);
  theAlphaInelasticProcess.AddDataSet(aShen);

  // Register the Alpha Inelastic and Binary Cascade Models
  theAlphaInelasticProcess.RegisterMe(theAlphaLEInelasticModel);
  theAlphaInelasticProcess.RegisterMe(theBinaryCascade);

  // Activate the alpha inelastic scattering using the parameterised Alpha Inelastic and Binary Cascade models
  processManager->AddDiscreteProcess(&theAlphaInelasticProcess);
  // Activate the Hadron Elastic Process
  processManager->AddDiscreteProcess(elasticScattering); 
}
