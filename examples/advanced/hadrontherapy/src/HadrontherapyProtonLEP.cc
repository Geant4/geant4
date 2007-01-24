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
// $Id: HadrontherapyProtonLEP.cc; May 2005
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

#include "HadrontherapyProtonLEP.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4LEProtonInelastic.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4PiNuclearCrossSection.hh" 
#include "G4LFission.hh"
#include "G4LCapture.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4HadronInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4BinaryLightIonReaction.hh" 
#include "G4HadronInelasticProcess.hh"

//
// PRECOMPOUND PHYSICS LIST
//
// LEP MODEL

HadrontherapyProtonLEP::HadrontherapyProtonLEP(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<<"****** Proton LEP Physics List is active !!!!!! ******"
	<<G4endl;
}

HadrontherapyProtonLEP::~HadrontherapyProtonLEP()
{}

void HadrontherapyProtonLEP::ConstructProcess()
{
 G4ParticleDefinition* particle = 0;
 G4ProcessManager* pmanager = 0;

 // Physics for proton, neutron, pion+ and pion-

  G4double LEPLowLimit = 0. * MeV;
  G4double LEPHighLimit = 300. * MeV;

  //--------------------------------------------------------------------------------------
  // Proton processes
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();
  
  // Model Registration  
  G4LEProtonInelastic* theLEProtonModel = new G4LEProtonInelastic();
  theLEProtonModel->SetMinEnergy(LEPLowLimit);
  theLEProtonModel->SetMaxEnergy(LEPHighLimit);
  protonInelasticProcess.RegisterMe(theLEProtonModel); 
  //protonInelasticProcess.AddDataSet(&protonInelasticCrossSection);
   
  // Activate the proton inelastic scattering using the precompound model
   pmanager -> AddDiscreteProcess(&protonInelasticProcess);


  //--------------------------------------------------------------------------------------
  // Neutron processes

  particle = G4Neutron::Neutron();
  pmanager = particle -> GetProcessManager();
		  
  // Register the LEP model
  G4LENeutronInelastic* theLENeutronModel = new G4LENeutronInelastic();
  theLENeutronModel->SetMinEnergy(LEPLowLimit);
  theLENeutronModel->SetMaxEnergy(LEPHighLimit);

  neutronInelasticProcess.RegisterMe(theLENeutronModel);
  //neutronInelasticProcess.AddDataSet(&neutronInelasticCrossSection);

  // Activate the neutron inelastic process
  pmanager -> AddDiscreteProcess(&neutronInelasticProcess);
   
 // Neutron capture process
 
  // Energy limits 
  G4HadronCaptureProcess* neutronCapture = new G4HadronCaptureProcess();

  G4LCapture* captureModel = new G4LCapture();

  // Energy limit of the neutron fission and capture
  G4double neutronLowLimit = 0.*TeV;
  G4double neutronHighLimit = 100.*TeV;
  // Set the energy range for the capture model
  captureModel -> SetMinEnergy(neutronLowLimit);
  captureModel -> SetMaxEnergy(neutronHighLimit);

  // Register the capture model
  neutronCapture -> RegisterMe(captureModel);
  // Active the capture process
  pmanager -> AddDiscreteProcess(neutronCapture);

  // Process for induced fission
  G4HadronFissionProcess* fission = new G4HadronFissionProcess();
  //Final state production model for induced fission
  G4LFission* fissionModel = new G4LFission();
  // Set the energy range for the fission model
  fissionModel -> SetMinEnergy(neutronLowLimit);
  fissionModel -> SetMaxEnergy(neutronHighLimit);
  // Register the fission model
  fission -> RegisterMe(fissionModel); 
  // Active the fission process
  pmanager -> AddDiscreteProcess(fission);  

  //--------------------------
  // Pion Plus

  particle = G4PionPlus::PionPlus(); 
  pmanager = particle -> GetProcessManager();
 
  // Define the inelastic process for pions plus
  G4PionPlusInelasticProcess* pionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");
  
  // Register the LEP model
  G4LEPionPlusInelastic* theLEPionPlusModel = new G4LEPionPlusInelastic();
  theLEPionPlusModel->SetMinEnergy(LEPLowLimit);
  theLEPionPlusModel->SetMaxEnergy(LEPHighLimit);
  pionPlusInelasticProcess -> RegisterMe(theLEPionPlusModel);
 
  // Active the inelastic process for pions plus
  pmanager -> AddDiscreteProcess(pionPlusInelasticProcess);
 
  //-----------------------------
  // Pion Minus

  particle = G4PionMinus::PionMinus();
  pmanager = particle -> GetProcessManager();

  // Define the inelastic process for pions minus
  G4PionMinusInelasticProcess* pionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");
  // Register the inelastic model for pion minus 

  G4LEPionMinusInelastic* theLEPionMinusModel = new G4LEPionMinusInelastic();
  theLEPionMinusModel->SetMinEnergy(LEPLowLimit);
  theLEPionMinusModel->SetMaxEnergy(LEPHighLimit);
  pionMinusInelasticProcess -> RegisterMe(theLEPionMinusModel);

  // Active the inelastic process for pion minus
  pmanager -> AddDiscreteProcess(pionMinusInelasticProcess); 
 
  //--------------------------------------------------------------------------------------
  // Physics for ions
 
  // Energy limit of the LEP model for ions
  LEPHighLimit = 200.*MeV;
 
  // TRIPATHI CROSS SECTION
  // Implementation of formulas in analogy to NASA technical paper 3621 by 
  // Tripathi, et al. Cross-sections for ion ion scattering
  G4TripathiCrossSection* TripathiCrossSection = new G4TripathiCrossSection;
  
  // IONS SHEN CROSS SECTION
  // Implementation of formulas 
  // Shen et al. Nuc. Phys. A 491 130 (1989) 
  // Total Reaction Cross Section for Heavy-Ion Collisions
  G4IonsShenCrossSection* aShen = new G4IonsShenCrossSection;

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
  
  // Active the deuteron inelastic scattering using the deuteron inelastic and binary cascade model
  pmanager -> AddDiscreteProcess(&deuteronInelasticProcess);
 
  // Triton
  particle = G4Triton::Triton();
  pmanager = particle -> GetProcessManager();
  
  // Final state production model for Triton inelastic scattering below 100 MeV
  G4LETritonInelastic* tritonLEModel = new G4LETritonInelastic;
  // Set the maximum energy for LEP model
  tritonLEModel -> SetMaxEnergy(LEPHighLimit);

  // Active the Tripathi and aShen Cross Section
  tritonInelasticProcess.AddDataSet(TripathiCrossSection);
  tritonInelasticProcess.AddDataSet(aShen);

  // Register the triton inelastic scattering models
  tritonInelasticProcess.RegisterMe(tritonLEModel);

  // Active the triton inelastic scattering process
  pmanager -> AddDiscreteProcess(&tritonInelasticProcess);
 
  // Alpha particles
  particle = G4Alpha::Alpha();
  pmanager = particle -> GetProcessManager();
		  
  // Final state production model for Alpha inelastic scattering below 20 GeV
  G4LEAlphaInelastic* alphaLEModel = new G4LEAlphaInelastic;
  // Set the maximum energy for LEP model
  alphaLEModel -> SetMaxEnergy(LEPHighLimit);
 
 // Register the alpha inelastic scattering models
  alphaInelasticProcess.AddDataSet(TripathiCrossSection);
  alphaInelasticProcess.AddDataSet(aShen);
 
  // Register the alpha inelastic scattering models
  alphaInelasticProcess.RegisterMe(alphaLEModel);
 
  // Active the alpha inelastic scattering
  pmanager -> AddDiscreteProcess(&alphaInelasticProcess);
 
  // He3
  // particle = G4He3::He3();

  // G4HadronInelasticProcess* He3inelasticProcess = 
  //   new G4HadronInelasticProcess("He3Inelastic",particle);
  
  //   G4BinaryLightIonReaction * ionBinaryCascade= new G4BinaryLightIonReaction;
 
  // He3inelasticProcess -> AddDataSet(TripathiCrossSection);
  //He3inelasticProcess -> AddDataSet(aShen);
  //He3inelasticProcess -> RegisterMe(ionBinaryCascade);
 
  //pmanager = particle -> GetProcessManager();
  //pmanager -> AddDiscreteProcess(He3inelasticProcess);
  //pmanager -> AddDiscreteProcess(elasticScattering); 

}



