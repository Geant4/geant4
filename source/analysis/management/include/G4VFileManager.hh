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

// Base class for File manager.
//
// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VFileManager_h
#define G4VFileManager_h 1

#include "G4BaseFileManager.hh"
#include "G4VTHnFileManager.hh"
#include "globals.hh"

// check whether directories exist
#include <sys/types.h>
#include <sys/stat.h>

namespace tools {
namespace histo { 
class h1d; 
class h2d; 
class h3d; 
class p1d; 
class p2d; 
}
}

class G4VFileManager : public G4BaseFileManager
{
  public:
    explicit G4VFileManager(const G4AnalysisManagerState& state);
    virtual ~G4VFileManager();
   
    // Method to open default file
    virtual G4bool OpenFile(const G4String& fileName) = 0;

    // Methods applied to file per name
    virtual G4bool CreateFile(const G4String& fileName) = 0;
    virtual G4bool WriteFile(const G4String& fileName) = 0;
    virtual G4bool CloseFile(const G4String& fileName) = 0;
    virtual G4bool SetIsEmpty(const G4String& fileName, G4bool isEmpty) = 0;

    // Methods applied to all registered files
    virtual G4bool WriteFiles() = 0;
    virtual G4bool CloseFiles() = 0;
    virtual G4bool DeleteEmptyFiles() = 0;

    // Methods for handling files and directories names
    virtual G4bool SetFileName(const G4String& fileName) final;
    G4bool SetHistoDirectoryName(const G4String& dirName);
    G4bool SetNtupleDirectoryName(const G4String& dirName);

    void LockDirectoryNames();

    G4bool IsOpenFile() const;
    G4String GetHistoDirectoryName() const;
    G4String GetNtupleDirectoryName() const;
    G4String GetHistoDirectoryNameIfExists() const;
    G4String GetNtupleDirectoryNameIfExists() const;

    // Access to helpers
    template <typename HT>
    std::shared_ptr<G4VTHnFileManager<HT>> GetHnFileManager() const;

  protected:
    // data members
    G4String fHistoDirectoryName;
    G4String fNtupleDirectoryName; 
    G4bool   fIsOpenFile;
    G4bool   fLockDirectoryNames;     

    // FileManagers per object type
    std::shared_ptr<G4VTHnFileManager<tools::histo::h1d>> fH1FileManager;
    std::shared_ptr<G4VTHnFileManager<tools::histo::h2d>> fH2FileManager;
    std::shared_ptr<G4VTHnFileManager<tools::histo::h3d>> fH3FileManager;
    std::shared_ptr<G4VTHnFileManager<tools::histo::p1d>> fP1FileManager;
    std::shared_ptr<G4VTHnFileManager<tools::histo::p2d>> fP2FileManager;
};

// inline functions

inline void G4VFileManager::LockDirectoryNames()
{ fLockDirectoryNames = true; }

inline G4bool G4VFileManager::IsOpenFile() const
{ return fIsOpenFile; }

inline G4String G4VFileManager::GetHistoDirectoryName() const 
{ return fHistoDirectoryName; }  

inline G4String G4VFileManager::GetNtupleDirectoryName() const 
{ return fNtupleDirectoryName; }

inline G4String G4VFileManager::GetHistoDirectoryNameIfExists() const 
{  
  G4String dirName = fHistoDirectoryName;

  struct stat info;

  int statRC = stat( dirName, &info );

  if (statRC==0 and info.st_mode and S_IFDIR) {
    return dirName;
  } else {
    return "";
  }
}

inline G4String G4VFileManager::GetNtupleDirectoryNameIfExists() const
{
  G4String dirName = fNtupleDirectoryName;

  struct stat info;

  int statRC = stat( dirName, &info );

  if (statRC==0 and info.st_mode and S_IFDIR) {
    return dirName;
  } else {
    return "";
  }
}

template <>
inline
std::shared_ptr<G4VTHnFileManager<tools::histo::h1d>> 
G4VFileManager::GetHnFileManager<tools::histo::h1d>() const
{   return fH1FileManager; }

template <>
inline
std::shared_ptr<G4VTHnFileManager<tools::histo::h2d>> 
G4VFileManager::GetHnFileManager<tools::histo::h2d>() const
{ return fH2FileManager; }

template <>
inline
std::shared_ptr<G4VTHnFileManager<tools::histo::h3d>> 
G4VFileManager::GetHnFileManager<tools::histo::h3d>() const
{ return fH3FileManager; }

template <>
inline
std::shared_ptr<G4VTHnFileManager<tools::histo::p1d>> 
G4VFileManager::GetHnFileManager<tools::histo::p1d>() const
{ return fP1FileManager; }

template <>
inline
std::shared_ptr<G4VTHnFileManager<tools::histo::p2d>> 
G4VFileManager::GetHnFileManager<tools::histo::p2d>() const
{ return fP2FileManager; }

#endif
