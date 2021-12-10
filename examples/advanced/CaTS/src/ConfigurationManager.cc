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
/// \file ConfigurationManager.cc
/// \brief Implementation of the CaTS::ConfigurationManager class

// Geant4 headers
#include <G4ios.hh>
#include "G4GenericMessenger.hh"
// project headers
#include "ConfigurationManager.hh"
// c++ headers
#include <ostream>
#include <string>
ConfigurationManager* ConfigurationManager::fginstance = nullptr;
std::once_flag ConfigurationManager::fginitInstanceFlag;

ConfigurationManager::ConfigurationManager() { DefineCommands(); }

ConfigurationManager::~ConfigurationManager() { delete fMessenger; }

void ConfigurationManager::Print()
{
  G4cout << "--------------------------------------------------" << G4endl;
  G4cout << "CaTS configuration: " << G4endl;
  G4cout << "====================" << G4endl;
  G4cout << G4endl;
  G4cout << "fenable_verbose:   " << fenable_verbose << G4endl;
  G4cout << "fdumpgdml:         " << fdumpgdml << G4endl;
  G4cout << "fGDMLFileName:     " << fGDMLFileName << G4endl;
#ifdef WITH_ROOT
  G4cout << "fdoAnalysis:       " << fdoAnalysis << G4endl;
  G4cout << "fHistoFileName:    " << fHistoFileName << G4endl;
  G4cout << "fwriteHits:        " << fwriteHits << G4endl;
  G4cout << "fname:            " << fname << G4endl;
#endif
  G4cout << "fenable_opticks:   " << fenable_opticks << G4endl;
  G4cout << "fMaxPhotons:       " << fMaxPhotons << G4endl;
  G4cout << "--------------------------------------------------" << G4endl;
}

void ConfigurationManager::DefineCommands()
{
  fMessenger = new G4GenericMessenger(this, "/CaTS/", "Configuring CaTS");
  //
  // Commands defining RootIO
#ifdef WITH_ROOT
  auto& HistoFileNameCmd =
    fMessenger->DeclareProperty("HistoFileName", fHistoFileName);
  HistoFileNameCmd.SetGuidance("Filename for Analysis Histograms");
  HistoFileNameCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  HistoFileNameCmd.SetDefaultValue("histograms.root");
  auto& writeHitsCmd = fMessenger->DeclareProperty("writeHits", fwriteHits);
  writeHitsCmd.SetGuidance("Write out Hits collection");
  writeHitsCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  writeHitsCmd.SetDefaultValue("false");
  auto& doAnalysisCmd = fMessenger->DeclareProperty("doAnalysis", fdoAnalysis);
  doAnalysisCmd.SetGuidance("Do Analysis");
  doAnalysisCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  doAnalysisCmd.SetDefaultValue("false");
  auto& fnameCmd = fMessenger->DeclareProperty("fname", fname);
  fnameCmd.SetGuidance("set fname");
  fnameCmd.SetStates(G4State_PreInit);
  fnameCmd.SetDefaultValue("Hits");
#endif
  //
  // Commands enabling G4Opticks and frequency calling it:
  //
  auto& enable_opticksCmd =
    fMessenger->DeclareProperty("enable_opticks", fenable_opticks);
  enable_opticksCmd.SetGuidance(
    "use opticks for generating and tracing of optical photons");
  enable_opticksCmd.SetDefaultValue("false");
  enable_opticksCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  auto& MaxPhotonsCmd = fMessenger->DeclareProperty("MaxPhotons", fMaxPhotons);
  MaxPhotonsCmd.SetGuidance(
    "set number of photons to be collecetd befores invoking G4Opticks");
  MaxPhotonsCmd.SetDefaultValue("1000000");
  MaxPhotonsCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  //
  // general command to control verbosity and writing out geometry to a gdml
  // file:
  //
  fMessenger->DeclareMethod("list", &ConfigurationManager::Print)
    .SetGuidance("Print all configuration parameters")
    .SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  auto& verboseCmd = fMessenger->DeclareProperty("verbose", fenable_verbose);
  verboseCmd.SetGuidance("Set flag for enabling verbose diagnostic printout");
  verboseCmd.SetDefaultValue("false");
  verboseCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  auto& GDMLFileNameCmd =
    fMessenger->DeclareProperty("GDMLFileName", fGDMLFileName);
  GDMLFileNameCmd.SetGuidance(
    "Set Filename to dump GDML representation of Geometry");
  GDMLFileNameCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  GDMLFileNameCmd.SetDefaultValue("dump.gdml_G4");
  auto& dumpgdmlCmd = fMessenger->DeclareProperty("dumpgdml", fdumpgdml);
  dumpgdmlCmd.SetGuidance(
    "Set flag for enabling dumping the Geometry to a gdml file");
  dumpgdmlCmd.SetStates(G4State_PreInit, G4State_Init, G4State_Idle);
  dumpgdmlCmd.SetDefaultValue("false");
}
