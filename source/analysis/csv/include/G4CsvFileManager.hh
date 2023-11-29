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

// The manager for Csv output file operations.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4CsvFileManager_h
#define G4CsvFileManager_h 1

#include "G4VTFileManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include "tools/wcsv_ntuple"

#include <string_view>

// Type aliases
using CsvNtupleDescription = G4TNtupleDescription<tools::wcsv::ntuple, std::ofstream>;

class G4CsvFileManager : public G4VTFileManager<std::ofstream>
{
  public:
    explicit G4CsvFileManager(const G4AnalysisManagerState& state);
    G4CsvFileManager() = delete;
    ~G4CsvFileManager() override = default;

    using G4BaseFileManager::GetNtupleFileName;
    using G4VTFileManager<std::ofstream>::WriteFile;
    using G4VTFileManager<std::ofstream>::CloseFile;

    // Methods to manipulate output files
    G4bool OpenFile(const G4String& fileName) final;

    G4bool SetHistoDirectoryName(const G4String& dirName) final;
    G4bool SetNtupleDirectoryName(const G4String& dirName) final;

    G4String GetFileType() const final { return "csv"; }

    // Specific methods for files per objects
    G4bool NotifyNtupleFile(CsvNtupleDescription* ntupleDescription);
    G4bool CreateNtupleFile(CsvNtupleDescription* ntupleDescription);
    G4bool CloseNtupleFile(CsvNtupleDescription* ntupleDescription);

    G4bool IsHistoDirectory() const;
    G4bool IsNtupleDirectory() const;

  protected:
    // Methods derived from templated base class
    std::shared_ptr<std::ofstream> CreateFileImpl(const G4String& fileName) final;
    G4bool WriteFileImpl(std::shared_ptr<std::ofstream> file) final;
    G4bool CloseFileImpl(std::shared_ptr<std::ofstream> file) final;

  private:
    // Utility method
    G4String GetNtupleFileName(CsvNtupleDescription* ntupleDescription);

    // Static data members
    static constexpr std::string_view fkClass { "G4CsvFileManager" };

    // Data members
    G4bool fIsHistoDirectory { false };
    G4bool fIsNtupleDirectory { false };
};

// inline functions

inline G4bool G4CsvFileManager::IsHistoDirectory() const
{ return fIsHistoDirectory; }

inline G4bool G4CsvFileManager::IsNtupleDirectory() const
{ return fIsNtupleDirectory; }

#endif
