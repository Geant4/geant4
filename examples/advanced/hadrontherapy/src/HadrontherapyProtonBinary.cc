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
// $Id: HadrontherapyProtonBinary.cc; May 2005
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

#include "HadrontherapyProtonBinary.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleTable.hh"
#include "G4Material.hh"
#include "G4LElastic.hh"
#include "G4BinaryCascade.hh"
#include "G4CascadeInterface.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionPlusInelastic.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionPlusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4PionMinusInelasticProcess.hh"
#include "G4PionPlusInelasticProcess.hh"
#include "G4LEPionMinusInelastic.hh"
#include "G4HEPionMinusInelastic.hh"
#include "G4PiMinusAbsorptionAtRest.hh"
#include "G4AntiProtonInelasticProcess.hh"
#include "G4LEAntiProtonInelastic.hh"
#include "G4HEAntiProtonInelastic.hh"
#include "G4AntiProtonAnnihilationAtRest.hh"

// BINARY + PRECOMPOUND PHYSICS LIST
//
// BINARY + PRECOMPOUND + EVAPORATION(DEFAULT EVAPORATION) NO FERMI BREAK-UP FOR PROTONS, NEUTRONS AND PIONS
// 
// LEP MODEL UP TO 100 MEV AND BINARY ION MODEL BETWEEN 80 MEV AND 40. GEV 
// FOR  DEUTERON, TRITON, HE3, ALPHA
// 
//FISSION AND HADRON CAPTURE FOR NEUTRONS BETWEEN 0. MEV AND 100. TEV
//
HadrontherapyProtonBinary::HadrontherapyProtonBinary(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  // Energy limits of the model for ions
  binaryLightIonLowLimit = 80.*MeV;
  binaryLightIonHighLimit = 40.*GeV;
  LEPHighLimit = 100.*MeV;

  // Energy limits for protons, neutrons and pions
  precompoundLowLimit = 100.*MeV;
  precompoundHighLimit = 300.*MeV;
  neutronLowLimit = 0.*TeV;
  neutronHighLimit = 100.*TeV;
  
  binaryLowLimit = 0.*MeV;
  binaryHighLimit = 300.*MeV;
}

HadrontherapyProtonBinary::~HadrontherapyProtonBinary()
{}

void HadrontherapyProtonBinary::ConstructProcess()
{

  // ELASTIC SCATTERING FOR PROTON, NEUTRON, IONS
  G4LElastic* elastic_Model = new G4LElastic();
  G4HadronElasticProcess* elastic = new G4HadronElasticProcess();
  elastic -> RegisterMe(elastic_Model);

  // PRECOMPOUND MODEL
  G4PreCompoundModel* thePreEquilib = new G4PreCompoundModel(&theHandler);
  thePreEquilib -> SetMinEnergy(precompoundLowLimit);
  thePreEquilib -> SetMaxEnergy(precompoundHighLimit);

  /////////////////////////////////////////////////////////////////////////////
  // HADRONIC PHYSICS FOR PROTONS
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* pmanager = 0;
  particle = G4Proton::Proton();
  pmanager = particle -> GetProcessManager();

  // INELASTIC SCATTERING: 
  //BINARY - PRECOMPOUND + DEFAULT EVAPORATION NO FERMI BREAK-UP  
  G4BinaryCascade* thePBC = new G4BinaryCascade();
  thePBC -> SetMinEnergy(binaryLowLimit);
  thePBC -> SetMaxEnergy(binaryHighLimit); 
  theIPProton.RegisterMe(thePBC);
  theIPProton.RegisterMe(thePreEquilib);
  theIPProton.AddDataSet(&thePXSec);
  pmanager -> AddDiscreteProcess(&theIPProton);

  // ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic); //ELASIC SCATTERING

  /////////////////////////////////////////////////////////////////////////////
  //////////  VERIFICARE LA POSSIBILITA' DEL BINARY PER PIONI
  /////////////////////////////////////////////////////////////////////////////
  // HADRONIC PHYSICS FOR PION PLUS
  particle = G4PionPlus::PionPlus();

  // ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic);

  // INELASTIC SCATTERING: BINARY - PRECOMPOUND + DEFAULT EVAPORATION NO FERMI BREAK-UP  
  pmanager -> AddDiscreteProcess(elastic);
  G4PionPlusInelasticProcess* thePionPlusInelasticProcess = new G4PionPlusInelasticProcess("inelastic");

  G4PreCompoundModel* theLEPionPlusInelasticModel = new G4PreCompoundModel(&theHandler);

  thePionPlusInelasticProcess -> RegisterMe(theLEPionPlusInelasticModel);

  pmanager -> AddDiscreteProcess(thePionPlusInelasticProcess);

 
  /////////////////////////////////////////////////////////////////////////////
  //HADRONIC PHYSICS FOR PION MINUS
  particle = G4PionMinus::PionMinus();

  // ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic);

  // INELASTIC SCATTERING: BINARY - PRECOMPOUND + DEFAULT EVAPORATION NO FERMI BREAK-UP 
  G4PionMinusInelasticProcess* thePionMinusInelasticProcess = new G4PionMinusInelasticProcess("inelastic");

  G4PreCompoundModel* theLEPionMinusInelasticModel = new G4PreCompoundModel(&theHandler);

  thePionMinusInelasticProcess -> RegisterMe(theLEPionMinusInelasticModel);

  pmanager -> AddDiscreteProcess(thePionMinusInelasticProcess);

/////////////////////////////////////////////////////////////////////////////
  // Neutron
  particle = G4Neutron::Neutron();
  pmanager = particle->GetProcessManager();
  thePreEquilib -> SetMinEnergy(precompoundLowLimit);
  thePreEquilib -> SetMaxEnergy(precompoundHighLimit); 
  G4BinaryCascade* theNBC = new G4BinaryCascade();
  theNBC->SetMinEnergy(binaryLowLimit);
  theNBC->SetMaxEnergy(binaryHighLimit); 
  theIPNeutron.RegisterMe(theNBC);
  theIPNeutron.RegisterMe(thePreEquilib);
  theIPNeutron.AddDataSet(&theNXSec);
  pmanager -> AddDiscreteProcess(&theIPNeutron);
  pmanager -> AddDiscreteProcess(elastic); // ELASTIC SCATTERING

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

  /////////////////////////////////////////////////////////////////////////////
  // ION HADRONIC PHYSICS LIST
  G4TripathiCrossSection * TripathiCrossSection= new G4TripathiCrossSection;
  G4IonsShenCrossSection * aShen = new G4IonsShenCrossSection;

  // INELASTIC SCATTERING FOR IONS
  G4BinaryLightIonReaction* theBC = new G4BinaryLightIonReaction();
  theBC -> SetMinEnergy(binaryLightIonLowLimit);
  theBC -> SetMaxEnergy(binaryLightIonHighLimit);

  /////////////////////////////////////////////////////////////////////////////
  // deuteron
  particle = G4Deuteron::Deuteron();
  pmanager = particle -> GetProcessManager();

  // INELASTIC SCATTERING
  G4LEDeuteronInelastic* theDIModel = new G4LEDeuteronInelastic;
  theDIModel -> SetMaxEnergy(LEPHighLimit);
  theIPdeuteron.AddDataSet(TripathiCrossSection);
  theIPdeuteron.AddDataSet(aShen);
  theIPdeuteron.RegisterMe(theDIModel);
  theIPdeuteron.RegisterMe(theBC);
  pmanager -> AddDiscreteProcess(&theIPdeuteron);

  //ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic); 

  /////////////////////////////////////////////////////////////////////////////
  // triton
  particle = G4Triton::Triton();
  pmanager = particle->GetProcessManager();
  // INELASTIC SCATTERING
  G4LETritonInelastic* theTIModel = new G4LETritonInelastic;
  theTIModel -> SetMaxEnergy(LEPHighLimit);
  theIPtriton.AddDataSet(TripathiCrossSection);
  theIPtriton.AddDataSet(aShen);
  theIPtriton.RegisterMe(theTIModel);
  theIPtriton.RegisterMe(theBC);
  pmanager -> AddDiscreteProcess(&theIPtriton);

  //ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic); //ELASTIC SCATTERING

  /////////////////////////////////////////////////////////////////////////////
  // alpha
  particle = G4Alpha::Alpha();
  pmanager = particle->GetProcessManager();

  // INELASTIC SCATTERING
  G4LEAlphaInelastic* theAIModel = new G4LEAlphaInelastic;
  theAIModel -> SetMaxEnergy(LEPHighLimit);
  theIPalpha.AddDataSet(TripathiCrossSection);
  theIPalpha.AddDataSet(aShen);
  theIPalpha.RegisterMe(theAIModel);
  theIPalpha.RegisterMe(theBC);
  pmanager -> AddDiscreteProcess(&theIPalpha);
  //ELASTIC SCATTERING
  pmanager -> AddDiscreteProcess(elastic); //ELASTIC SCATTERING

  /////////////////////////////////////////////////////////////////////////////
  // He3
  particle = G4He3::He3();
  pmanager = particle->GetProcessManager();
  G4BinaryLightIonReaction * theGenIonBC= new G4BinaryLightIonReaction;
  G4HadronInelasticProcess* theIPHe3 =
    new G4HadronInelasticProcess("He3Inelastic",particle);
  theIPHe3 -> AddDataSet(TripathiCrossSection);
  theIPHe3 -> AddDataSet(aShen);
  theIPHe3 -> RegisterMe(theGenIonBC);
  pmanager -> AddDiscreteProcess(theIPHe3);
  pmanager -> AddDiscreteProcess(elastic); //ELASTIC SCATTERING
  
  
}



