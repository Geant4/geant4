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
// $Id: G4VAnalysisManager.hh 71635 2013-06-19 13:48:28Z ihrivnac $

// The nonvirtual public interface class to g4tools based analysis.
// It is defined as a composite of object manager base classes.
// Individual use of the component managers is disabled
// (except for file manager and Hn manager which are also used from
//  other object managers).
// The functions which has to be implemented in concrete managers 
// are declared as virtual protected.

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#ifndef G4VAnalysisManager_h
#define G4VAnalysisManager_h 1

#include "G4AnalysisManagerState.hh"
#include "globals.hh"

#include <vector>
#include <fstream>


class G4AnalysisMessenger;
class G4HnManager;
class G4VH1Manager;
class G4VH2Manager;
class G4VNtupleManager;
class G4VFileManager;

class G4VAnalysisManager
{
  public:
    G4VAnalysisManager(const G4String& type, G4bool isMaster);
    virtual ~G4VAnalysisManager();
   
    // Methods for handling files 
    G4bool OpenFile();
    G4bool OpenFile(const G4String& fileName);
    G4bool Write(); 
    G4bool CloseFile(); 

    // Methods for handling files and directories names    
    G4bool SetFileName(const G4String& fileName);
    G4bool SetHistoDirectoryName(const G4String& dirName);
    G4bool SetNtupleDirectoryName(const G4String& dirName);
    G4String GetFileName() const;
    G4String GetHistoDirectoryName() const;
    G4String GetNtupleDirectoryName() const;
    
    // Methods for handling histogrammes
    G4int CreateH1(const G4String& name, const G4String& title,
                   G4int nbins, G4double xmin, G4double xmax,
                   const G4String& unitName = "none",
                   const G4String& fcnName = "none",
                   const G4String& binSchemeName = "linear");

    G4int CreateH1(const G4String& name, const G4String& title,
                   const std::vector<G4double>& edges,
                   const G4String& unitName = "none",
                   const G4String& fcnName = "none");

    G4int CreateH2(const G4String& name, const G4String& title,
                   G4int nxbins, G4double xmin, G4double xmax, 
                   G4int nybins, G4double ymin, G4double ymax,
                   const G4String& xunitName = "none", 
                   const G4String& yunitName = "none",
                   const G4String& xfcnName = "none", 
                   const G4String& yfcnName = "none",
                   const G4String& xbinScheme = "linear",
                   const G4String& ybinScheme = "linear");
   
    G4int CreateH2(const G4String& name, const G4String& title,
                   const std::vector<G4double>& xedges,
                   const std::vector<G4double>& yedges,
                   const G4String& xunitName = "none", 
                   const G4String& yunitName = "none",
                   const G4String& xfcnName = "none", 
                   const G4String& yfcnName = "none");
                         
    G4bool SetH1(G4int id,
                   G4int nbins, G4double xmin, G4double xmax,
                   const G4String& unitName = "none",
                   const G4String& fcnName = "none",
                   const G4String& binSchemeName = "linear");

    G4bool SetH1(G4int id,
                   const std::vector<G4double>& edges,
                   const G4String& unitName = "none",
                   const G4String& fcnName = "none");

    G4bool SetH2(G4int id,
                   G4int nxbins, G4double xmin, G4double xmax, 
                   G4int nybins, G4double ymin, G4double ymax,
                   const G4String& xunitName = "none", 
                   const G4String& yunitName = "none",
                   const G4String& xfcnName = "none", 
                   const G4String& yfcnName = "none",
                   const G4String& xbinSchemeName = "linear",
                   const G4String& ybinSchemeName = "linear");

    G4bool SetH2(G4int id,
                   const std::vector<G4double>& xedges,
                   const std::vector<G4double>& yedges,
                   const G4String& xunitName = "none", 
                   const G4String& yunitName = "none",
                   const G4String& xfcnName = "none", 
                   const G4String& yfcnName = "none");

    G4bool ScaleH1(G4int id, G4double factor);
    G4bool ScaleH2(G4int id, G4double factor);
                           
    // Methods for handling ntuples
    G4int CreateNtuple(const G4String& name, const G4String& title);
    // Create columns in the last created ntuple
    G4int CreateNtupleIColumn(const G4String& name);
    G4int CreateNtupleFColumn(const G4String& name);
    G4int CreateNtupleDColumn(const G4String& name);
    void  FinishNtuple();   
    // Create columns in the ntuple with given id
    G4int CreateNtupleIColumn(G4int ntupleId, const G4String& name);
    G4int CreateNtupleFColumn(G4int ntupleId, const G4String& name);
    G4int CreateNtupleDColumn(G4int ntupleId, const G4String& name);
    void  FinishNtuple(G4int ntupleId); 
    
    // The ids of histograms and ntuples are generated automatically
    // starting from 0; with following functions it is possible to
    // change the first Id to start from other value
    G4bool SetFirstHistoId(G4int firstId);
    G4bool SetFirstH1Id(G4int firstId);
    G4bool SetFirstH2Id(G4int firstId);
    G4bool SetFirstNtupleId(G4int firstId);
    G4bool SetFirstNtupleColumnId(G4int firstId);
  
    // Methods to fill histogrammes, ntuples
    G4bool FillH1(G4int id, G4double value, G4double weight = 1.0);
    G4bool FillH2(G4int id, G4double xvalue, G4double yvalue,
                  G4double weight = 1.0);
    // Methods for ntuple with id = FirstNtupleId                     
    G4bool FillNtupleIColumn(G4int id, G4int value);
    G4bool FillNtupleFColumn(G4int id, G4float value);
    G4bool FillNtupleDColumn(G4int id, G4double value);
    G4bool AddNtupleRow();
    // Methods for ntuple with id > FirstNtupleId (when more ntuples exist)                      
    G4bool FillNtupleIColumn(G4int ntupleId, G4int columnId, G4int value);
    G4bool FillNtupleFColumn(G4int ntupleId, G4int columnId, G4float value);
    G4bool FillNtupleDColumn(G4int ntupleId, G4int columnId, G4double value);
    G4bool AddNtupleRow(G4int ntupleId);
    
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
    G4int GetNofH1s() const;
    G4int GetNofH2s() const;
    G4int GetNofNtuples() const;

    // Access methods via names
    G4int GetH1Id(const G4String& name, G4bool warn = true) const;
    G4int GetH2Id(const G4String& name, G4bool warn = true) const;

   
    // Methods to manipulate additional information
    void  SetH1Activation(G4bool activation);
    void  SetH1Activation(G4int id, G4bool activation);
    void  SetH1Ascii(G4int id, G4bool ascii);
    void  SetH2Activation(G4bool activation);
    void  SetH2Activation(G4int id, G4bool activation);
    void  SetH2Ascii(G4int id, G4bool ascii);
    
    // Access to H1 parameters
    G4int    GetH1Nbins(G4int id) const;
    G4double GetH1Xmin(G4int id) const;
    G4double GetH1Xmax(G4int id) const;
    G4double GetH1Width(G4int id) const;
    
    // Access to H2 parameters
    G4int    GetH2Nxbins(G4int id) const;
    G4double GetH2Xmin(G4int id) const;
    G4double GetH2Xmax(G4int id) const;
    G4double GetH2XWidth(G4int id) const;
    G4int    GetH2Nybins(G4int id) const;
    G4double GetH2Ymin(G4int id) const;
    G4double GetH2Ymax(G4int id) const;
    G4double GetH2YWidth(G4int id) const;

    // Access to H1 additional information
    G4String GetH1Name(G4int id) const;
    G4double GetH1Unit(G4int id) const;
    G4bool   GetH1Activation(G4int id) const;
    G4bool   GetH1Ascii(G4int id) const;

    // Access to H2 additional information
    G4String GetH2Name(G4int id) const;
    G4double GetH2XUnit(G4int id) const;
    G4double GetH2YUnit(G4int id) const;
    G4bool   GetH2Activation(G4int id) const;
    G4bool   GetH2Ascii(G4int id) const;
    
    // Setters for attributes for plotting
    G4bool SetH1Title(G4int id, const G4String& title);
    G4bool SetH1XAxisTitle(G4int id, const G4String& title);
    G4bool SetH1YAxisTitle(G4int id, const G4String& title);
    G4bool SetH2Title(G4int id, const G4String& title);
    G4bool SetH2XAxisTitle(G4int id, const G4String& title);
    G4bool SetH2YAxisTitle(G4int id, const G4String& title);
    G4bool SetH2ZAxisTitle(G4int id, const G4String& title);

    // Access attributes for plotting
    G4String GetH1Title(G4int id) const;
    G4String GetH1XAxisTitle(G4int id) const;
    G4String GetH1YAxisTitle(G4int id) const;
    G4String GetH2Title(G4int id) const;
    G4String GetH2XAxisTitle(G4int id) const;
    G4String GetH2YAxisTitle(G4int id) const;
    G4String GetH2ZAxisTitle(G4int id) const;

    // Verbosity
    void  SetVerboseLevel(G4int verboseLevel);
    G4int GetVerboseLevel() const;

    // The manager type (starts with an uppercase letter)
    G4String GetType() const;                 
    // The manager file type (starts with a lowercase letter)
    G4String GetFileType() const;                 
   
  protected:
    // virtual methods
    virtual G4bool OpenFileImpl(const G4String& fileName) = 0;
    virtual G4bool WriteImpl() = 0; 
    virtual G4bool CloseFileImpl() = 0; 
  
    // methods
    void SetH1Manager(G4VH1Manager* h1Manager);
    void SetH2Manager(G4VH2Manager* h2Manager);
    void SetNtupleManager(G4VNtupleManager* ntupleManager);
    void SetFileManager(G4VFileManager* fileManager);

    // Methods to manipulate additional information
    G4bool  WriteAscii(const G4String& fileName); 

    // data members
    G4AnalysisManagerState fState;

  private:
    // data members
    G4AnalysisMessenger* fMessenger;
    G4HnManager*      fH1HnManager;
    G4HnManager*      fH2HnManager;
    G4VH1Manager*     fVH1Manager;
    G4VH2Manager*     fVH2Manager;
    G4VNtupleManager* fVNtupleManager;
    G4VFileManager*   fVFileManager;
};

// inline functions

#include "G4VAnalysisManager.icc"
 
#endif

