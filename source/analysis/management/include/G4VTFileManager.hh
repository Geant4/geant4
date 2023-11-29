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

// The common implementaion of G4VFileManager functions via G4TFileManager.

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VTFileManager_h
#define G4VTFileManager_h 1

#include "G4VFileManager.hh"
#include "G4TNtupleDescription.hh"
#include "G4TFileManager.hh"
#include "globals.hh"

#include "tools/wcsv_ntuple"

template <typename FT>
class G4VTFileManager : public G4VFileManager,
                        public G4TFileManager<FT>
{
  public:
    explicit G4VTFileManager(const G4AnalysisManagerState& state)
      : G4VFileManager(state), G4TFileManager<FT>(state) {}
    ~G4VTFileManager() override = default;

    using G4VFileManager::WriteFile;
    using G4VFileManager::CloseFile;

    // Methods applied to file per name
    G4bool CreateFile(const G4String& fileName) final;
    G4bool WriteFile(const G4String& fileName) final;
    G4bool CloseFile(const G4String& fileName) final;
    G4bool SetIsEmpty(const G4String& fileName, G4bool isEmpty) final;

    // Methods applied to all registered files
    G4bool OpenFiles() final;
    G4bool WriteFiles() final;
    G4bool CloseFiles() final;
    G4bool DeleteEmptyFiles() final;

    // Clear all data
    void Clear() final;

    // Get method
    std::shared_ptr<FT> GetFile() const;

  protected:
    // Data members
    // Default file - created at OpenFile call
    std::shared_ptr<FT> fFile { nullptr };
};

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::CreateFile(const G4String& fileName)
{
  return (G4TFileManager<FT>::CreateTFile(fileName) != nullptr);
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::WriteFile(const G4String& fileName)
{
  return G4TFileManager<FT>::WriteTFile(fileName);
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::CloseFile(const G4String& fileName)
{
  return G4TFileManager<FT>::CloseTFile(fileName);
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::SetIsEmpty(const G4String& fileName, G4bool isEmpty)
{
  return G4TFileManager<FT>::SetIsEmpty(fileName, isEmpty);
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::OpenFiles()
{
  return G4TFileManager<FT>::OpenFiles();
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::WriteFiles()
{
  return G4TFileManager<FT>::WriteFiles();
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::CloseFiles()
{
  auto result = G4TFileManager<FT>::CloseFiles();

  fIsOpenFile = false;
  fFile.reset();

  return result;
}

//_____________________________________________________________________________
template <typename FT>
inline
G4bool G4VTFileManager<FT>::DeleteEmptyFiles()
{
  return G4TFileManager<FT>::DeleteEmptyFiles();
}

//_____________________________________________________________________________
template <typename FT>
inline
void G4VTFileManager<FT>::Clear()
{
  G4TFileManager<FT>::ClearData();
  UnlockDirectoryNames();
}

//_____________________________________________________________________________
template <typename FT>
inline std::shared_ptr<FT> G4VTFileManager<FT>::GetFile() const
{
  return fFile;
}

#endif
