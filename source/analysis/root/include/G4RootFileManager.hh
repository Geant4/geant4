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
#include "globals.hh"

#include <vector>
#include <memory>
#include <tuple>

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
    virtual ~G4RootFileManager();

    using G4BaseFileManager::GetNtupleFileName;
    using G4VTFileManager<G4RootFile>::WriteFile;
    using G4VTFileManager<G4RootFile>::CloseFile;

    // Methods to manipulate file from base classes
    virtual G4bool OpenFile(const G4String& fileName) final;

    virtual G4String GetFileType() const final { return "root"; }

    // Specific methods for files per objects
    std::shared_ptr<G4RootFile> CreateNtupleFile(RootNtupleDescription* ntupleDescription, 
                                  G4int mainNumber = -1);
    std::shared_ptr<G4RootFile> GetNtupleFile(RootNtupleDescription* ntupleDescription, 
                                  G4bool perThread = true,
                                  G4int mainNumber = -1) const;
    G4bool WriteNtupleFile(RootNtupleDescription* ntupleDescription);
    G4bool CloseNtupleFile(RootNtupleDescription* ntupleDescription); 

    // Set methods
    void  SetBasketSize(unsigned int basketSize);
    void  SetBasketEntries(unsigned int basketEntries);

    // Get methods
    unsigned int GetBasketSize() const;
    unsigned int GetBasketEntries() const;

  protected:
    // Methods derived from templated base class
    virtual std::shared_ptr<G4RootFile> CreateFileImpl(const G4String& fileName) final;
    virtual G4bool WriteFileImpl(std::shared_ptr<G4RootFile> file) final;
    virtual G4bool CloseFileImpl(std::shared_ptr<G4RootFile> file) final;    

  private:
    // methods
    tools::wroot::directory* CreateDirectory(
                               std::shared_ptr<tools::wroot::file> rfile,
                               const G4String& directoryName, 
                               const G4String& objectType) const;
    G4String GetNtupleFileName(
                RootNtupleDescription* ntupleDescription, 
                G4bool perThread = true,
                G4int mainNumber = -1) const;

    // data members
    unsigned int fBasketSize;
    unsigned int fBasketEntries;
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

