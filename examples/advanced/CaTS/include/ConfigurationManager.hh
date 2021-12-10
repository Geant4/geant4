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
/// \file ConfigurationManager.hh
/// \brief Definition of the CaTS::ConfigurationManager class

#pragma once

#include <mutex>
#include "G4String.hh"
class G4GenericMessenger;

class ConfigurationManager
{
 private:
  static ConfigurationManager* fginstance;
  static std::once_flag fginitInstanceFlag;
#ifdef WITH_ROOT
  G4bool fdoAnalysis{ false };  // variable determines if we are doing analysis
  G4String fHistoFileName{
    "histograms.root"
  };  // File name for histos and  ntuples
  G4bool fwriteHits{
    false
  };  // variable determines if hits are written out into Root File
  G4String fname{ "Hits" };  // full File name for root io
#endif
  G4bool fenable_opticks{ true };  // use opticks if available
  unsigned int fMaxPhotons{ 1000000 };
  G4bool fenable_verbose{ false };  // switch on/off diagnostic printouts
  G4bool fdumpgdml{ false };        // write out Detector to gdml file
  G4String fGDMLFileName{ "dump.gdml_G4" };
  ConfigurationManager();

 public:
  ~ConfigurationManager();
  /// Define UI commands: choice of primary generators
  void DefineCommands();
  /// Pointer to the messenger for UI commands
  G4GenericMessenger* fMessenger = nullptr;
  static ConfigurationManager* getInstance()
  {
    std::call_once(fginitInstanceFlag,
                   ConfigurationManager::initConfigurationManager);
    return fginstance;
  }
  static void initConfigurationManager()
  {
    fginstance = new ConfigurationManager();
  }
#ifdef WITH_ROOT
  inline G4String getHistoFileName() const { return fHistoFileName; }
  inline G4bool isWriteHits() const { return fwriteHits; }
  inline G4bool isdoAnalysis() const { return fdoAnalysis; }
  inline void setfname(G4String name) { fname = name; }
  inline G4String getfname() const { return fname; }
#endif
  inline G4bool isEnable_opticks() const { return fenable_opticks; };
  inline unsigned int getMaxPhotons() const { return fMaxPhotons; }
  inline G4bool isEnable_verbose() const { return fenable_verbose; };
  inline G4String getGDMLFileName() const { return fGDMLFileName; }
  inline void setGDMLFileName(G4String GDMLFileName)
  {
    fGDMLFileName = GDMLFileName;
  }
  inline G4bool isDumpgdml() const { return fdumpgdml; }
  void Print();
};
