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

// The manager for Hdf5 file input operations.

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5RFileManager_h
#define G4Hdf5RFileManager_h 1

#include "G4VRFileManager.hh"
#include "globals.hh"

#include "toolx/hdf5/ntuple"

#include <map>
#include <tuple>
#include <string_view>

using G4Hdf5File = std::tuple<hid_t, hid_t, hid_t>;

class G4Hdf5RFileManager : public G4VRFileManager
{
  public:
    explicit G4Hdf5RFileManager(const G4AnalysisManagerState& state);
    G4Hdf5RFileManager() = delete;
    ~G4Hdf5RFileManager() override = default;

    G4String GetFileType() const final { return "hdf5"; }

    // Methods from base class
    void CloseFiles() final;

    // Get methods
    G4Hdf5File* GetRFile(const G4String& fileName, G4bool isPerThread);
    hid_t GetHistoRDirectory(const G4String& fileName, const G4String& dirName,
                   G4bool isPerThread);
    hid_t GetNtupleRDirectory(const G4String& fileName, const G4String& dirName,
                    G4bool isPerThread);

  private:
    // Methods
    hid_t OpenRFile(const G4String& fileName, G4bool isPerThread);
    hid_t OpenDirectory(hid_t file, const G4String& directoryName);
    hid_t GetRDirectory(const G4String& directoryType,
                   const G4String& fileName, const G4String& dirName,
                   G4bool isPerThread);

    // Static data members
    static constexpr std::string_view fkClass { "G4Hdf5RFileManager" };
    inline static const G4String fgkDefaultDirectoryName { "default" };

    // Data members
    std::map<G4String, G4Hdf5File> fRFiles;
};

#endif
