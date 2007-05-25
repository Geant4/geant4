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

#include "HIProtonNeutronLEP.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4HadronCaptureProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4ProtonInelasticCrossSection.hh"
#include "G4NeutronInelasticCrossSection.hh"
#include "G4LEProtonInelastic.hh"
#include "G4LENeutronInelastic.hh"
#include "G4LFission.hh"
#include "G4LCapture.hh"


HIProtonNeutronLEP::HIProtonNeutronLEP(const G4String& name): 
  G4VPhysicsConstructor(name)
{
  G4cout<< "HADRONIC INELASTIC PROCESS(ES): G4XXXInelasticProcess (protons, neutrons)" 
        << G4endl
        << "APPLIED MODEL(S): G4LEXXXInelastic" 
        << G4endl;
}

HIProtonNeutronLEP::~HIProtonNeutronLEP()
{}

void HIProtonNeutronLEP::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

  // ******************************************
  // *** Proton,Neutron: Common Definitions ***
  // ******************************************

  G4double protonNeutronLEPMinEnergy = 0. * MeV;
  G4double protonNeutronLEPMaxEnergy = 100. * MeV;


  // **************
  // *** Proton ***
  // **************

  G4ProtonInelasticProcess* protonInelasticProcess = new G4ProtonInelasticProcess(); 
  G4ProtonInelasticCrossSection* protonInelasticCrossSection =  new G4ProtonInelasticCrossSection(); 

  G4LEProtonInelastic* protonLEPModel = new G4LEProtonInelastic();
  protonLEPModel -> SetMinEnergy(protonNeutronLEPMinEnergy);
  protonLEPModel -> SetMaxEnergy(protonNeutronLEPMaxEnergy);

  protonInelasticProcess -> RegisterMe(protonLEPModel);
  protonInelasticProcess -> AddDataSet(protonInelasticCrossSection);
   
  particle = G4Proton::Proton();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(protonInelasticProcess);


  // ***************
  // *** Neutron ***
  // ***************

  G4NeutronInelasticProcess* neutronInelasticProcess = new G4NeutronInelasticProcess; 		  
  G4NeutronInelasticCrossSection* neutronInelasticCrossSection = new G4NeutronInelasticCrossSection; 
  G4HadronCaptureProcess* neutronCaptureProcess = new G4HadronCaptureProcess();
  G4HadronFissionProcess* neutronFissionProcess = new G4HadronFissionProcess();

  G4LENeutronInelastic* neutronLEPModel = new G4LENeutronInelastic();
  neutronLEPModel -> SetMinEnergy(protonNeutronLEPMinEnergy);
  neutronLEPModel -> SetMaxEnergy(protonNeutronLEPMaxEnergy);

  G4double neutronCaptureFissionMinEnergy = 0. * TeV;
  G4double neutronCaptureFissionMaxEnergy = 100. * TeV;

  G4LCapture* neutronCaptureModel = new G4LCapture();
  neutronCaptureModel -> SetMinEnergy(neutronCaptureFissionMinEnergy);
  neutronCaptureModel -> SetMaxEnergy(neutronCaptureFissionMaxEnergy);

  G4LFission* neutronFissionModel = new G4LFission();
  neutronFissionModel -> SetMinEnergy(neutronCaptureFissionMinEnergy);
  neutronFissionModel -> SetMaxEnergy(neutronCaptureFissionMaxEnergy);

  neutronInelasticProcess -> RegisterMe(neutronLEPModel);
  neutronInelasticProcess -> AddDataSet(neutronInelasticCrossSection);
  neutronCaptureProcess -> RegisterMe(neutronCaptureModel);
  neutronFissionProcess -> RegisterMe(neutronFissionModel); 

  particle = G4Neutron::Neutron();
  processManager = particle -> GetProcessManager();
  processManager -> AddDiscreteProcess(neutronInelasticProcess);
  processManager -> AddDiscreteProcess(neutronCaptureProcess); 
  processManager -> AddDiscreteProcess(neutronFissionProcess);  

}



