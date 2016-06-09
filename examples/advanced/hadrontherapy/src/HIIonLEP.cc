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
// $Id: HadrontherapyProtonPrecompound.cc; November 2008
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

#include "HIIonLEP.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4TripathiCrossSection.hh"
#include "G4IonsShenCrossSection.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"
#include "G4LEDeuteronInelastic.hh"
#include "G4LETritonInelastic.hh"
#include "G4LEAlphaInelastic.hh"


HIIonLEP::HIIonLEP(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<< "HADRONIC INELASTIC PROCESS(ES): G4XXXInelasticProcess (all ions)" 
        << G4endl
        << "APPLIED MODEL(S): G4LEXXXInelastic" 
        << G4endl;
}

HIIonLEP::~HIIonLEP()
{}

void HIIonLEP::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

   
  // ***************************************************
  // *** Deuteron, Triton, Alpha: Common Definitions ***
  // ***************************************************

  G4TripathiCrossSection* ionTripathiCrossSection = new G4TripathiCrossSection;
  G4IonsShenCrossSection* ionShenCrossSection = new G4IonsShenCrossSection;

  G4double ionLEPMaxEnergy = 100. * MeV;


  // ****************
  // *** Deuteron ***
  // ****************

  G4DeuteronInelasticProcess* deuteronInelasticProcess = new G4DeuteronInelasticProcess;

  G4LEDeuteronInelastic* deuteronLEPModel = new G4LEDeuteronInelastic;
  deuteronLEPModel -> SetMaxEnergy(ionLEPMaxEnergy);

  deuteronInelasticProcess -> AddDataSet(ionTripathiCrossSection);
  deuteronInelasticProcess -> RegisterMe(deuteronLEPModel);

  particle = G4Deuteron::Deuteron();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(deuteronInelasticProcess);


  // **************
  // *** Triton ***
  // **************

  G4TritonInelasticProcess* tritonInelasticProcess = new G4TritonInelasticProcess;
  
  G4LETritonInelastic* tritonLEPModel = new G4LETritonInelastic;
  tritonLEPModel -> SetMaxEnergy(ionLEPMaxEnergy);

  tritonInelasticProcess -> AddDataSet(ionTripathiCrossSection);
  tritonInelasticProcess -> RegisterMe(tritonLEPModel);

  particle = G4Triton::Triton();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(tritonInelasticProcess);
 

  // *************
  // *** Alpha ***
  // *************

  G4AlphaInelasticProcess* alphaInelasticProcess = new G4AlphaInelasticProcess;

  G4LEAlphaInelastic* alphaLEPModel = new G4LEAlphaInelastic;
  alphaLEPModel -> SetMaxEnergy(ionLEPMaxEnergy);
 
  alphaInelasticProcess -> AddDataSet(ionTripathiCrossSection);
  alphaInelasticProcess -> RegisterMe(alphaLEPModel);
 
  particle = G4Alpha::Alpha();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(alphaInelasticProcess);
  
}



