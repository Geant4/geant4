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
#include "G4HnInformation.hh"
#include "globals.hh"

#include <vector>
#include <fstream>

class G4AnalysisMessenger;

class G4VAnalysisManager
{
  public:
    // Supported object types
    enum ObjectType { kH1, kH2, kNtuple };

  public:
    G4VAnalysisManager(const G4String& type = "");
    virtual ~G4VAnalysisManager();
   
    // Methods to manipulate files
    virtual G4bool OpenFile();
    virtual G4bool OpenFile(const G4String& fileName) = 0;
    virtual G4bool Write() = 0; 
    virtual G4bool CloseFile() = 0; 
    
    // Methods for handling files and directories names
    virtual G4bool SetFileName(const G4String& fileName);
    virtual G4bool SetHistoDirectoryName(const G4String& dirName);
    virtual G4bool SetNtupleDirectoryName(const G4String& dirName);
    virtual G4String GetFileName() const;
    virtual G4String GetFullFileName() const;
    virtual G4String GetHistoDirectoryName() const;
    virtual G4String GetNtupleDirectoryName() const;
    
    // Methods for handling histogrammes, ntuples
    virtual G4int CreateH1(const G4String& name, const G4String& title,
                           G4int nbins, G4double xmin, G4double xmax,
                           G4double unit = 1.0) = 0;
    virtual G4int CreateH2(const G4String& name, const G4String& title,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4double unit = 1.0, G4double yunit = 1.0) = 0;
                           
    virtual G4bool SetH1(G4int id,
                           G4int nbins, G4double xmin, G4double xmax,
                           G4double unit = 1.0) = 0;
    virtual G4bool SetH2(G4int id,
                           G4int nxbins, G4double xmin, G4double xmax, 
                           G4int nybins, G4double ymin, G4double ymax,
                           G4double xunit = 1.0, G4double yunit = 1.0) = 0;
                           
    virtual void  CreateNtuple(const G4String& name, const G4String& title) = 0;
    virtual G4int CreateNtupleIColumn(const G4String& name) = 0;
    virtual G4int CreateNtupleFColumn(const G4String& name) = 0;
    virtual G4int CreateNtupleDColumn(const G4String& name) = 0;
    virtual void  FinishNtuple() = 0;   
    
    // The ids of histograms and ntuples are generated automatically
    // starting from 0; with following functions it is possible to
    // change the first Id to start from other value
    virtual G4bool SetFirstHistoId(G4int firstId);
    virtual G4bool SetFirstNtupleColumnId(G4int firstId);
  
    // Methods to fill histogrammes, ntuples
    virtual G4bool FillH1(G4int id, G4double value, G4double weight) = 0;
    virtual G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                          G4double weight) = 0;
    virtual G4bool FillNtupleIColumn(G4int id, G4int value) = 0;
    virtual G4bool FillNtupleFColumn(G4int id, G4float value) = 0;
    virtual G4bool FillNtupleDColumn(G4int id, G4double value) = 0;
    virtual G4bool AddNtupleRow() = 0;
    
    // Activation option
    
    // When this option is enabled, only the histograms marked as activated
    // are returned, filled or saved on file.
    // No warning is issued when Get or Fill is called on inactive histogram.
    void   SetActivation(G4bool activation);
    G4bool GetActivation() const;

    // Return false if activation is enabled and there is no object activated,
    // return true otherwise
    G4bool IsActive() const;

    // ASCII option

    // Return false if there is no object selected for ASCII output,
    // return true otherwise
    G4bool IsAscii() const;

    // Access methods
    virtual G4int GetNofH1s() const;
    virtual G4int GetNofH2s() const;
    
    // Methods to manipulate additional information

    // Setters for additional information by fields
    // (other fields are set by SetH1, SetH2 methods)
    void  SetActivation(ObjectType type, G4int id, G4bool activation);
    void  SetAscii(ObjectType type, G4int id, G4bool ascii);

    // Access to additional information by fields
    G4String GetName(ObjectType type, G4int id) const;
    G4double GetXUnit(ObjectType type, G4int id) const;
    G4double GetYUnit(ObjectType type, G4int id) const;
    G4bool   GetActivation(ObjectType type, G4int id) const;
    G4bool   GetAscii(ObjectType type, G4int id) const;
    
    // Verbosity
    virtual G4int GetVerboseLevel() const;
    virtual void  SetVerboseLevel(G4int verboseLevel);

    // The manager type (starts with an uppercase letter)
    G4String GetType() const;                 
    // The manager file type (starts with a lowercase letter)
    G4String GetFileType() const;                 
    
  protected:
    // methods to manipulate additional information
    void   AddH1Information(const G4String& name, 
                            G4double unit);
    void   AddH2Information(const G4String& name,
                            G4double xunit, G4double yunit); 
    
    // Methods to manipulate additional information
    G4HnInformation* GetH1Information(G4int id) const;
    G4HnInformation* GetH2Information(G4int id) const;
    G4HnInformation* GetInformation(ObjectType type, G4int id) const;
    
    G4bool  WriteAscii(); 
    virtual G4bool WriteOnAscii(std::ofstream& output) = 0;
    
    // data members
    G4int    fVerboseLevel;
    G4bool   fActivation;
    G4int    fFirstHistoId;
    G4int    fFirstNtupleColumnId;
    G4String fFileName;
    G4String fHistoDirectoryName;
    G4String fNtupleDirectoryName; 
    G4bool   fLockFirstHistoId;     
    G4bool   fLockFirstNtupleColumnId;     
    G4bool   fLockFileName;     
    G4bool   fLockHistoDirectoryName;     
    G4bool   fLockNtupleDirectoryName;

    // Additional histograms properties not included in tools
    G4AnalysisVerbose  fVerboseL1;
    G4AnalysisVerbose  fVerboseL2;
    G4AnalysisVerbose  fVerboseL3;
    G4AnalysisVerbose* fpVerboseL1;
    G4AnalysisVerbose* fpVerboseL2;
    G4AnalysisVerbose* fpVerboseL3;
 
  private:
    G4AnalysisMessenger* fMessenger;    
    G4int   fNofActiveObjects;
    G4int   fNofAsciiObjects;

    // Additional histograms/ntuple properties not included in tools
    std::vector<G4HnInformation*> fH1Informations;
    std::vector<G4HnInformation*> fH2Informations;
};

// inline functions

inline G4String G4VAnalysisManager::GetFileName() const {
  return fFileName;
}  

inline G4String G4VAnalysisManager::GetHistoDirectoryName() const {
  return fHistoDirectoryName;
}  

inline G4String G4VAnalysisManager::GetNtupleDirectoryName() const {
  return fNtupleDirectoryName;
}  

inline G4int G4VAnalysisManager::GetVerboseLevel() const {
  return fVerboseLevel;
}  

inline  G4String G4VAnalysisManager::GetType() const {
  return fVerboseL1.GetType();
}                 

inline  void  G4VAnalysisManager::SetActivation(G4bool activation) {
  fActivation = activation;
}

inline G4bool G4VAnalysisManager::GetActivation() const {
  return fActivation;
}  

inline  G4int  G4VAnalysisManager::GetNofH1s() const { 
  return fH1Informations.size();
}
  
inline  G4int  G4VAnalysisManager::GetNofH2s() const {
  return fH2Informations.size();
}
  
#endif

