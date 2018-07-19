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
// $Id$

// The manager for Hdf5 file output operations.

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5FileManager_h
#define G4Hdf5FileManager_h 1

#include "G4VFileManager.hh"
#include "G4TNtupleDescription.hh"
#include "globals.hh"

#include "tools/hdf5/ntuple"

#include <fstream>
#include <memory>

class G4AnalysisManagerState;

class G4Hdf5FileManager : public G4VFileManager
{
  public:
    explicit G4Hdf5FileManager(const G4AnalysisManagerState& state);
    ~G4Hdf5FileManager();

    // Type aliases
    using NtupleType = tools::hdf5::ntuple;
    using NtupleDescriptionType = G4TNtupleDescription<NtupleType>;

    // Methods to manipulate output files
    virtual G4bool OpenFile(const G4String& fileName) final;
    virtual G4bool WriteFile() final;
    virtual G4bool CloseFile() final; 

    G4bool WriteHistoDirectory();
    G4bool WriteNtupleDirectory();
    void   CloseAfterHnWrite();

    // Set methods
    void  SetBasketSize(unsigned int basketSize);

    // Get methods
    hid_t GetFile() const;
    hid_t GetHistoDirectory() const;
    hid_t GetNtupleDirectory() const;
    unsigned int GetBasketSize() const; 
    
   private:
    G4bool CreateDirectory(const G4String& directoryType, 
                           const G4String& directoryName, hid_t& directory);
    G4bool WriteDirectory(const G4String& directoryType, 
                          const G4String& directoryName, hid_t& directory);

    // constants
    static const G4String fgkDefaultDirectoryName;
    
    // data members
    hid_t  fFile;
    hid_t  fHistoDirectory;
    hid_t  fNtupleDirectory;
    unsigned int fBasketSize;
};

// inline functions

inline void  G4Hdf5FileManager::SetBasketSize(unsigned int basketSize)  
{ fBasketSize = basketSize; }

inline hid_t G4Hdf5FileManager::GetFile() const
{ return fFile; }

inline hid_t G4Hdf5FileManager::GetHistoDirectory() const
{ return fHistoDirectory; }

inline hid_t G4Hdf5FileManager::GetNtupleDirectory() const
{ return fNtupleDirectory; }

inline unsigned int G4Hdf5FileManager::GetBasketSize() const
{ return fBasketSize; }

#endif
