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

// The manager for Hdf5 file output operations.

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5FileManager_h
#define G4Hdf5FileManager_h 1

#include "G4VTFileManager.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include "toolx/hdf5/ntuple" // for hid_t

#include <memory>
#include <tuple>
#include <string_view>

using G4Hdf5File = std::tuple<hid_t, hid_t, hid_t>;
using Hdf5NtupleDescription = G4TNtupleDescription<toolx::hdf5::ntuple, G4Hdf5File>;

class G4Hdf5FileManager : public G4VTFileManager<G4Hdf5File>
{
  public:
    explicit G4Hdf5FileManager(const G4AnalysisManagerState& state);
    G4Hdf5FileManager() = delete;
    ~G4Hdf5FileManager() override = default;

    using G4BaseFileManager::GetNtupleFileName;
    using G4VTFileManager<G4Hdf5File>::WriteFile;
    using G4VTFileManager<G4Hdf5File>::CloseFile;

    // Methods to manipulate output files
    G4bool OpenFile(const G4String& fileName) final;

    G4String GetFileType() const final { return "hdf5"; }

    // Specific methods for files per objects
    G4bool CreateNtupleFile(Hdf5NtupleDescription* ntupleDescription);
    G4bool CloseNtupleFile(Hdf5NtupleDescription* ntupleDescription);

    // Set methods
    void  SetBasketSize(unsigned int basketSize);

    // Get methods
    hid_t GetHistoDirectory() const;
    hid_t GetNtupleDirectory() const;
    unsigned int GetBasketSize() const;

  protected:
    // // Methods derived from base class
    std::shared_ptr<G4Hdf5File> CreateFileImpl(const G4String& fileName) final;
    G4bool WriteFileImpl(std::shared_ptr<G4Hdf5File> file) final;
    G4bool CloseFileImpl(std::shared_ptr<G4Hdf5File> file) final;

  private:
    hid_t CreateDirectory(hid_t& file, const G4String& directoryName,
             const G4String& objectType);
    G4String GetNtupleFileName(Hdf5NtupleDescription* ntupleDescription);

    // Static data members
    static constexpr std::string_view fkClass { "G4Hdf5FileManager" };
    inline static const G4String fgkDefaultDirectoryName { "default" };

    // data members
    unsigned int fBasketSize { G4Analysis::kDefaultBasketSize };
};

// inline functions

//_____________________________________________________________________________
inline void G4Hdf5FileManager::SetBasketSize(unsigned int basketSize)
{ fBasketSize = basketSize; }

//_____________________________________________________________________________
inline unsigned int G4Hdf5FileManager::GetBasketSize() const
{ return fBasketSize; }

#endif
