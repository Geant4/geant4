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
// ----------------------------------------------------------------------------

#include "HadrontherapyProtonBertini.hh"
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
HadrontherapyProtonBertini::HadrontherapyProtonBertini(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout << "The Bertini model is set for protons, neutrons and pions !!!!" << G4endl;
  // Inelastic process, energy limits 

  //
  // The Bertini model is set for protons, neutrons and pions
  // This model contains a pre-equilibrium model and a de-excitation model
  // Energy limit of the Bertini model
  bertiniLowLimit = 0.*MeV;
  bertiniHighLimit = 300.*MeV;

  // Energy limit of the neutron fission and capture
  neutronLowLimit = 0.*TeV;
  neutronHighLimit = 100.*TeV;
 
  // Ions:
  // The inelastic scattering is modelled with  LEP model up to 100 MeV,
  // then Binary Ion Model 
  // Energy limit of the LEP model for ions
  LEPHighLimit = 100.*MeV;
  // Energy limit of the binary ion model
  binaryLightIonLowLimit = 80.*MeV;
  binaryLightIonHighLimit = 40.*GeV;
}

HadrontherapyProtonBertini::~HadrontherapyProtonBertini()
{}

void HadrontherapyProtonBertini::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;

  // LOW ENERGY ELASTIC SCATTERING 
  // FOR PROTON, NEUTRON, IONS 
  G4LElastic* elastic_model = new G4LElastic();
  G4HadronElasticProcess* elastic_scattering = new G4HadronElasticProcess();
  elastic_scattering -> RegisterMe(elastic_model);
  
  // INELASTIC SCATTERING
  // Bertini Model for protons, pions and neutrons
  G4CascadeInterface * theBertiniModel = new G4CascadeInterface;
  // Set the min and max energy for the Bertini Model
  theBertiniModel -> SetMinEnergy(bertiniLowLimit);
  theBertiniModel -> SetMaxEnergy(bertiniHighLimit);

  // Binary Cascade for deuteron, triton, alpha particle, He3
  G4BinaryLightIonReaction* theBinaryCascade = new G4BinaryLightIonReaction();
  // Set the min and max energy for the Binary Cascade
  theBinaryCascade -> SetMinEnergy(binaryLightIonLowLimit);
  theBinaryCascade -> SetMaxEnergy(binaryLightIonHighLimit);

  // TRIPATHI CROSS SECTION
  // Implementation of formulas in analogy to NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* TripathiCrossSection = new G4TripathiCrossSection;
  
  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

//--------------------------------------------------------------------------------------
  // Proton BERTINI MODEL
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  // Model Registration 
  theIPProton.RegisterMe(theBertiniModel);
  // Active the Cross-sections for proton nuclear scattering up to 20 GeV
  theIPProton.AddDataSet(&thePXSec);

  // Active the proton inelastic scattering 
  pmanager -> AddDiscreteProcess(&theIPProton);
  // Active the Hadron Elastic Process 
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // deuteron 
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // Final state production model for deuteron inelastic scattering below 100 MeV
  G4LEDeuteronInelastic* theDIModel = new G4LEDeuteronInelastic;
  // Set the maximum energy for LEP model
  theDIModel -> SetMaxEnergy(LEPHighLimit);
  // Active the Tripathi and aShen Cross Section
  theIPdeuteron.AddDataSet(TripathiCrossSection);
  theIPdeuteron.AddDataSet(aShen);

  // Register the Parameterised Deuteron Inelastic Model and the Ion Binary Cascade Model
  theIPdeuteron.RegisterMe(theDIModel);
  theIPdeuteron.RegisterMe(theBinaryCascade);
  // Active the deuteron elastic and inelastic scattering 
  pmanager -> AddDiscreteProcess(&theIPdeuteron);
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering); 

 // triton
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  
  // Final state production model for Triton inelastic scattering below 100 MeV
  G4LETritonInelastic* theTIModel = new G4LETritonInelastic;
  // Set the maximum energy for LEP model
  theTIModel -> SetMaxEnergy(LEPHighLimit);
  // Active the Tripathi and aShen Cross Section
  theIPtriton.AddDataSet(TripathiCrossSection);
  theIPtriton.AddDataSet(aShen);
  // Register the Triton Inelastic and Binary Cascade Model
  theIPtriton.RegisterMe(theTIModel);
  theIPtriton.RegisterMe(theBinaryCascade);
  // Active the triton inelastic scattering using the triton inelastic and binary cascade model
  pmanager -> AddDiscreteProcess(&theIPtriton);
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering);

  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle->GetProcessManager();
  // Final state production model for Alpha inelastic scattering below 20 GeV
  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  // Set the maximum energy for LEP model
  theAIModel -> SetMaxEnergy(LEPHighLimit);
  // Register the Triton Inelastic and Binary Cascade Model
  theIPalpha.AddDataSet(TripathiCrossSection);
  theIPalpha.AddDataSet(aShen);
  // Register the Alpha Inelastic and Binary Cascade Model
  theIPalpha.RegisterMe(theAIModel);
  theIPalpha.RegisterMe(theBinaryCascade);
  // Active the alpha inelastic scattering using the alpha inelastic and binary cascade model
  pmanager -> AddDiscreteProcess(&theIPalpha);
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering); 

  // He3
  // particle = G4He3::He3();
//   pmanager = particle->GetProcessManager();
  
//   // Binary Cascade inelastic scattering for ions
//   G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
//   // Inelastic Scattering for ions
//   G4HadronInelasticProcess* theIPHe3 = new G4HadronInelasticProcess("He3Inelastic",particle);
//   // Active the Tripathi and aShen Cross Section
//   theIPHe3 -> AddDataSet(TripathiCrossSection);
//   theIPHe3 -> AddDataSet(aShen);
//   // Register the Alpha Binary Cascade Model
//   theIPHe3 -> RegisterMe(theGenIonBC);
//   // Active the Inelastic Process for He3
//   pmanager -> AddDiscreteProcess(theIPHe3);
//   // Active the Hadron Elastic Process
//   pmanager -> AddDiscreteProcess(elastic_scattering); 

// Neutron processes
  particle = G4Neutron::Neutron();
  pmanager = particle->GetProcessManager();
  // Register the Precompound model
  theIPNeutron.RegisterMe(theBertiniModel);
  // Active the Cross-sections for neutron nuclear scattering from 14 MeV up to 20 GeV
  theIPNeutron.AddDataSet(&theNXSec);
  // Active the neutron inelastic process
  pmanager -> AddDiscreteProcess(&theIPNeutron);
  // Active the Hadron Elastic Process
  pmanager -> AddDiscreteProcess(elastic_scattering);

// Pions plus processes
  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
  // Define the inelastic process for pions plus
  G4PionPlusInelasticProcess* thePionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");
  // Register the Low Energy Inelastic Model for pions plus
  thePionPlusInelasticProcess -> RegisterMe(theBertiniModel);
  // Active the inelastic process for pions plus
  pmanager->AddDiscreteProcess(thePionPlusInelasticProcess);
  pmanager -> AddDiscreteProcess(elastic_scattering);

  // Pion Minus processes
  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();
  // Define the inelastic process for pions minus
  G4PionMinusInelasticProcess* thePionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
  // Register the inelastic model for pion minus  
  thePionMinusInelasticProcess -> RegisterMe(theBertiniModel);
  // Active the inelastic process for pion minus
  pmanager -> AddDiscreteProcess(thePionMinusInelasticProcess); 
  // Active Absorption process for pion minus
  pmanager -> AddRestProcess(new G4PiMinusAbsorptionAtRest, ordDefault);
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



