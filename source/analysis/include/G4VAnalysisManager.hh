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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifndef G4VAnalysisManager_h
#define G4VAnalysisManager_h 1

#include "G4AnalysisVerbose.hh"
#include "globals.hh"

class G4VAnalysisManager
{
  public:
    G4VAnalysisManager(const G4String& type = "");
    virtual ~G4VAnalysisManager();
   
    // Methods to manipulate files
    virtual G4bool OpenFile(const G4String& fileName) = 0;
    virtual G4bool Write() = 0; 
    virtual G4bool CloseFile() = 0; 
    
    // Methods for handling directories
    virtual void SetHistoDirectoryName(const G4String& dirName);
    virtual void SetNtupleDirectoryName(const G4String& dirName);
    
    // Methods for handling histogrammes, ntuples
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax) = 0;
    virtual G4int CreateH2(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax) = 0;
                           
    virtual void  CreateNtuple(const G4String& name, const G4String& title) = 0;
    virtual G4int CreateNtupleIColumn(const G4String& name) = 0;
    virtual G4int CreateNtupleFColumn(const G4String& name) = 0;
    virtual G4int CreateNtupleDColumn(const G4String& name) = 0;
    virtual void  FinishNtuple() = 0;   
    
    // The ids of histograms and ntuples are generated automatically
    // starting from 0; with following functions it is possible to
    // change the first Id to start from other value
    virtual void SetFirstHistoId(G4int firstId);
    virtual void SetFirstNtupleId(G4int firstId);
  
    // Methods to fill histogrammes, ntuples
    virtual G4bool FillH1(G4int id, G4double value, G4double weight) = 0;
    virtual G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                          G4double weight) = 0;
    virtual G4bool FillNtupleIColumn(G4int id, G4int value) = 0;
    virtual G4bool FillNtupleFColumn(G4int id, G4float value) = 0;
    virtual G4bool FillNtupleDColumn(G4int id, G4double value) = 0;
    virtual G4bool AddNtupleRow() = 0;
    
    // Verbosity
    virtual G4int GetVerboseLevel() const;
    virtual void  SetVerboseLevel(G4int verboseLevel);
    
  protected:
    G4int    fVerboseLevel;
    G4int    fFirstHistoId;
    G4int    fFirstNtupleId;
    G4String fHistoDirectoryName;
    G4String fNtupleDirectoryName;      
    G4AnalysisVerbose  fVerboseL1;
    G4AnalysisVerbose  fVerboseL2;
    G4AnalysisVerbose* fpVerboseL1;
    G4AnalysisVerbose* fpVerboseL2;
};

// inline functions

inline void G4VAnalysisManager::SetHistoDirectoryName(const G4String& dirName) {
  fHistoDirectoryName = dirName;
}  

inline void G4VAnalysisManager::SetNtupleDirectoryName(const G4String& dirName) {
  fNtupleDirectoryName = dirName;
}  

inline void G4VAnalysisManager::SetFirstHistoId(G4int firstId) {
  fFirstHistoId = firstId;
}  

inline void G4VAnalysisManager::SetFirstNtupleId(G4int firstId) {
  fFirstNtupleId = firstId;
}

inline G4int G4VAnalysisManager::GetVerboseLevel() const {
  return fVerboseLevel;
}  

#endif

