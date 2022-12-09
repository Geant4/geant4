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

#ifndef G4RootFileManager_h
#define G4RootFileManager_h 1

#include "G4VTFileManager.hh"
#include "G4RootFileDef.hh"
#include "G4AnalysisUtilities.hh"
#include "globals.hh"

#include <vector>
#include <memory>
#include <tuple>
#include <string_view>

namespace tools {
namespace wroot{
class ntuple;
}
}

// Types alias
using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;

class G4RootFileManager : public G4VTFileManager<G4RootFile>
{
  public:
    explicit G4RootFileManager(const G4AnalysisManagerState& state);
    G4RootFileManager() = delete;
    ~G4RootFileManager() override = default;

    using G4BaseFileManager::GetNtupleFileName;
    using G4VTFileManager<G4RootFile>::WriteFile;
    using G4VTFileManager<G4RootFile>::CloseFile;

    // Methods to manipulate file from base classes
    G4bool OpenFile(const G4String& fileName) final;

    G4String GetFileType() const final { return "root"; }
    G4bool HasCycles() const final { return true; }

    // Specific methods for files per objects
    std::shared_ptr<G4RootFile> CreateNtupleFile(RootNtupleDescription* ntupleDescription,
                                  G4int mainNumber = -1);
    std::shared_ptr<G4RootFile> GetNtupleFile(RootNtupleDescription* ntupleDescription,
                                  G4bool perThread = true,
                                  G4int mainNumber = -1) const;
    G4bool CloseNtupleFile(RootNtupleDescription* ntupleDescription,
                                  G4int mainNumber = -1);

    // Set methods
    void  SetBasketSize(unsigned int basketSize);
    void  SetBasketEntries(unsigned int basketEntries);

    // Get methods
    unsigned int GetBasketSize() const;
    unsigned int GetBasketEntries() const;

  protected:
    // Methods derived from templated base class
    std::shared_ptr<G4RootFile> CreateFileImpl(const G4String& fileName) final;
    G4bool WriteFileImpl(std::shared_ptr<G4RootFile> file) final;
    G4bool CloseFileImpl(std::shared_ptr<G4RootFile> file) final;

  private:
    // Methods
    tools::wroot::directory* CreateDirectory(
                               tools::wroot::file* rfile,
                               const G4String& directoryName,
                               const G4String& objectType) const;
    G4String GetNtupleFileName(
                RootNtupleDescription* ntupleDescription,
                G4bool perThread = true,
                G4int mainNumber = -1) const;

    // Static data members
    static constexpr std::string_view fkClass { "G4RootFileManager" };

    // Data members
    unsigned int fBasketSize { G4Analysis::kDefaultBasketSize };
    unsigned int fBasketEntries { G4Analysis::kDefaultBasketEntries };
};

// inline functions

//_____________________________________________________________________________
inline void  G4RootFileManager::SetBasketSize(unsigned int basketSize)
{ fBasketSize = basketSize; }

//_____________________________________________________________________________
inline void  G4RootFileManager::SetBasketEntries(unsigned int basketEntries)
{ fBasketEntries = basketEntries; }

//_____________________________________________________________________________
inline unsigned int G4RootFileManager::GetBasketSize() const
{ return fBasketSize; }

//_____________________________________________________________________________
inline unsigned int G4RootFileManager::GetBasketEntries() const
{ return fBasketEntries; }

#endif

