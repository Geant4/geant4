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

// The manager for Root output file operations.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4GenericFileManager_h
#define G4GenericFileManager_h 1

#include "G4VFileManager.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include <vector>
#include <memory>
#include <string_view>

class G4HnInformation;
class G4VNtupleFileManager;

class G4CsvFileManager;
#ifdef TOOLS_USE_HDF5
class G4Hdf5FileManager;
#endif
class G4RootFileManager;
class G4XmlFileManager;

class G4GenericFileManager : public G4VFileManager
{
  public:
    explicit G4GenericFileManager(const G4AnalysisManagerState& state);
    ~G4GenericFileManager() override = default;

    // Method to manipulate a default file
    G4bool OpenFile(const G4String& fileName) final;

    // Methods applied to all registered files
    virtual G4bool OpenFiles() final;
    G4bool WriteFiles() final;
    G4bool CloseFiles() final;
    G4bool DeleteEmptyFiles() final;

    // Clear all data
    void Clear() final;

    // Methods applied to file per name
    G4bool CreateFile(const G4String& fileName) final;
    G4bool WriteFile(const G4String& fileName) final;
    G4bool CloseFile(const G4String& fileName) final;
    G4bool SetIsEmpty(const G4String& fileName, G4bool isEmpty) final;

    G4bool SetHistoDirectoryName(const G4String& dirName) override;
    G4bool SetNtupleDirectoryName(const G4String& dirName) override;

    G4String GetFileType() const final { return ""; }

    // Set default output type (backward compatibility)
    // this type will be used for file names without extension
    void SetDefaultFileType(const G4String& value);
    G4String GetDefaultFileType() const;

    // Get file manager for a specific output type
    std::shared_ptr<G4VFileManager> GetFileManager(const G4String& fileName);

    template <typename HT>
    G4bool WriteTExtra(const G4String& fileName, HT* ht, const G4String& htName);

    // Methods
    std::shared_ptr<G4VNtupleFileManager> CreateNtupleFileManager(G4AnalysisOutput output);

  private:
    // Methods
    void CreateFileManager(G4AnalysisOutput output);
    std::shared_ptr<G4VFileManager> GetFileManager(G4AnalysisOutput output) const;

    // Static data members
    static constexpr std::string_view fkClass { "G4GenericFileManager" };
    // inline static const G4String fgkDefaultFileType { "root" };
    // inline static const G4String fgkDefaultFileType { "undefined" };

    // Data members
    G4String fDefaultFileType;
    std::shared_ptr<G4VFileManager> fDefaultFileManager { nullptr };
    std::vector<std::shared_ptr<G4VFileManager>> fFileManagers {
       // Csv,  Hdf5,    Root,    Xml
       nullptr, nullptr, nullptr, nullptr
     };
    std::shared_ptr<G4CsvFileManager>  fCsvFileManager { nullptr };
#ifdef TOOLS_USE_HDF5
    std::shared_ptr<G4Hdf5FileManager> fHdf5FileManager { nullptr };
#endif
    std::shared_ptr<G4RootFileManager> fRootFileManager { nullptr };
    std::shared_ptr<G4XmlFileManager>  fXmlFileManager { nullptr };
    G4bool fHdf5Warn { true };
};

#include "G4GenericFileManager.icc"

#endif

