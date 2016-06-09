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
// $Id: HadrontherapyHadronUElastic.cc; May 2005
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

#include "HEHadronIonUElastic.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4UHadronElasticProcess.hh"
#include "G4HadronElastic.hh"
#include "G4VQCrossSection.hh"


HEHadronIonUElastic::HEHadronIonUElastic(const G4String& name): 
   G4VPhysicsConstructor(name)
{
  G4cout<< "HADRONIC ELASTIC PROCESS(ES): G4UHadronElasticProcess (all considered hadrons and ions)" 
        << G4endl
        << "APPLIED MODEL(S): G4HadronElastic" 
        << G4endl;
 
  model = 0;
}

HEHadronIonUElastic::~HEHadronIonUElastic()
{
  delete model;
}

void HEHadronIonUElastic::ConstructProcess()
{
  G4ParticleDefinition* particle = 0;
  G4ProcessManager* processManager = 0;

  // **********************************************
  // *** Proton, Neutron, Pion plus, Pion minus ***
  // *** Deuteron, Triton, Alpha                ***
  // **********************************************

  G4UHadronElasticProcess* hadronIonElasticProcess = new G4UHadronElasticProcess("elasticProcess", false);     
  
  G4double hadronIonElasticMinEnergy = 0. * keV;
  G4double hadronIonElasticMaxEnergy = 100. * TeV;

  G4HadronElastic* hadronIonElasticModel = new G4HadronElastic();
  hadronIonElasticModel -> SetMinEnergy(hadronIonElasticMinEnergy);
  hadronIonElasticModel -> SetMaxEnergy(hadronIonElasticMaxEnergy);

  G4VQCrossSection* hadronIonElasticCrossSection = hadronIonElasticModel->GetCS();

  hadronIonElasticProcess -> SetQElasticCrossSection(hadronIonElasticCrossSection);
  hadronIonElasticProcess -> RegisterMe(hadronIonElasticModel);

  theParticleIterator -> reset();

  while( (*theParticleIterator)() )
  {
    particle = theParticleIterator -> value();
    G4String particleName = particle -> GetParticleName();

    if(particleName == "neutron"   || 
       particleName == "pi-"       || 
       particleName == "pi+"       || 
       particleName == "proton"    || 
       particleName == "alpha"     ||
       particleName == "deuteron"  ||
       particleName == "triton") 
      
      {
	processManager = particle -> GetProcessManager();
	processManager -> AddDiscreteProcess(hadronIonElasticProcess);
      } 
  }

}



