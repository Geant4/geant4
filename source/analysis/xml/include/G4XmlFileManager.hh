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

// The manager for Xml file output operations.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4XmlFileManager_h
#define G4XmlFileManager_h 1

#include "G4VTFileManager.hh"
#include "G4TNtupleDescription.hh"
#include "G4TFileManager.hh"
#include "globals.hh"

#include "tools/waxml/ntuple"

#include <fstream>
#include <memory>
#include <string_view>

// Type aliases
using XmlNtupleDescription = G4TNtupleDescription<tools::waxml::ntuple, std::ofstream>;

class G4AnalysisManagerState;

class G4XmlFileManager : public G4VTFileManager<std::ofstream>
{
  public:
    explicit G4XmlFileManager(const G4AnalysisManagerState& state);
    G4XmlFileManager() = delete;
    ~G4XmlFileManager() override = default;

    using G4BaseFileManager::GetNtupleFileName;
    using G4VTFileManager<std::ofstream>::WriteFile;
    using G4VTFileManager<std::ofstream>::CloseFile;

    // Methods to manipulate output files
    G4bool OpenFile(const G4String& fileName) final;

    G4String GetFileType() const final { return "xml"; }

    // Specific methods for files per objects
    G4bool CreateNtupleFile(XmlNtupleDescription* ntupleDescription);
    G4bool CloseNtupleFile(XmlNtupleDescription* ntupleDescription);

  protected:
    // Methods derived from templated base class
    std::shared_ptr<std::ofstream> CreateFileImpl(const G4String& fileName) final;
    G4bool WriteFileImpl(std::shared_ptr<std::ofstream> file) final;
    G4bool CloseFileImpl(std::shared_ptr<std::ofstream> file) final;

  private:
    // Utility method
    G4String GetNtupleFileName(XmlNtupleDescription* ntupleDescription);

    // Static data members
    static constexpr std::string_view fkClass { "G4XmlFileManager" };
};

#endif
