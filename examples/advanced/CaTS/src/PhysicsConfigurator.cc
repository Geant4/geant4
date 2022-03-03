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
// ********************************************************************
//
//  CaTS (Calorimetry and Tracking Simulation)
//
//  Authors : Hans Wenzel
//            Soon Yung Jun
//            (Fermi National Accelerator Laboratory)
//
// History
//   October 18th, 2021 : first implementation
//
// ********************************************************************
//
/// \file PhysicsConfigurator.cc
/// \brief Implementation of the CaTS::PhysicsConfigurator class

// Geant4 headers
#include "G4String.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhysListFactoryAlt.hh"
#include "G4PhysicsConstructorRegistry.hh"
#include "G4PhysListRegistry.hh"
#include "G4OpticalParameters.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4StepLimiter.hh"
#include "G4StepLimiterPhysics.hh"
#include "G4SystemOfUnits.hh"
// project Headers
#include "PhysicsConfigurator.hh"
#include "ConfigurationManager.hh"
// c++ headers
#include <stdlib.h> 
PhysicsConfigurator* PhysicsConfigurator::instance = 0;

G4VModularPhysicsList* PhysicsConfigurator::Construct(G4String physName)
{
  //
  // Access to registries and factories
  //
  G4PhysicsConstructorRegistry* g4pcr =
    G4PhysicsConstructorRegistry::Instance();
  G4PhysListRegistry* g4plr = G4PhysListRegistry::Instance();
  G4bool verbose = ConfigurationManager::getInstance()->isEnable_verbose();
  if(verbose)
  {
    G4cout << "Available Physics Constructors:  "
           << g4pcr->AvailablePhysicsConstructors().size() << G4endl;
    G4cout << "Available Physics Lists:         "
           << g4plr->AvailablePhysLists().size() << G4endl;
    G4cout << "Available Physics Extensions:    "
           << g4plr->AvailablePhysicsExtensions().size() << G4endl;
    G4cout << "Available Physics Lists Em:      "
           << g4plr->AvailablePhysListsEM().size() << G4endl;
    g4plr->SetVerbose(1);
  }
  else
  {
    g4plr->SetVerbose(0);
  }
  g4plr->AddPhysicsExtension("OPTICAL", "G4OpticalPhysics");
  g4plr->AddPhysicsExtension("STEPLIMIT", "G4StepLimiterPhysics");
  g4plr->AddPhysicsExtension("NEUTRONLIMIT", "G4NeutronTrackingCut");
  if(verbose)
  {
    g4pcr->PrintAvailablePhysicsConstructors();
    g4plr->PrintAvailablePhysLists();
  }
  g4alt::G4PhysListFactory factory;
  G4VModularPhysicsList* phys = nullptr;
  if(verbose)
    G4cout << "Physics configuration: " << physName << G4endl;
  //
  // currently using the Constructor names doesn't work otherwise it would be:
  // G4String physName = "FTFP_BERT+G4OpticalPhysics+G4StepLimiterPhysics";
  // using the name doesn't work either
  // G4String physName = "FTFP_BERT+Optical+stepLimiter";
  // reference PhysicsList via its name
  //
  if(factory.IsReferencePhysList(physName))
  {
    phys = factory.GetReferencePhysList(physName);
  }
  else
  {
    G4cout << "Not a reference physics list" << G4endl;
    g4plr->PrintAvailablePhysLists();
    exit(EXIT_FAILURE);
  }
  if(verbose)
  {
    G4cout << phys->GetPhysicsTableDirectory() << G4endl;
  }
  G4OpticalParameters::Instance()->SetProcessActivation("Cerenkov", true);
  G4OpticalParameters::Instance()->SetProcessActivation("Scintillation", true);
  G4OpticalParameters::Instance()->SetProcessActivation("OpAbsorption", true);
  G4OpticalParameters::Instance()->SetProcessActivation("OpRayleigh", true);
  G4OpticalParameters::Instance()->SetProcessActivation("OpMieHG", false);
  G4OpticalParameters::Instance()->SetProcessActivation("OpWLS", true);
  G4OpticalParameters::Instance()->SetProcessActivation("OpWLS2", false);

  G4OpticalParameters::Instance()->SetCerenkovStackPhotons(false);
  G4OpticalParameters::Instance()->SetScintStackPhotons(false);
  G4OpticalParameters::Instance()->SetScintTrackSecondariesFirst(
    true);  // only relevant if we actually stack and trace the optical photons
  G4OpticalParameters::Instance()->SetCerenkovTrackSecondariesFirst(
    true);  // only relevant if we actually stack and trace the optical photons
  G4OpticalParameters::Instance()->SetCerenkovMaxPhotonsPerStep(100);
  G4OpticalParameters::Instance()->SetCerenkovMaxBetaChange(10.0);
  if(verbose)
  {
    phys->DumpList();
  }
  return phys;
}

PhysicsConfigurator* PhysicsConfigurator::getInstance()
{
  if(instance == 0)
    instance = new PhysicsConfigurator();
  return instance;
}
