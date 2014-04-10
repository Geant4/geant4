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
// $Id: G4VAnalysisReader.hh 71635 2013-06-19 13:48:28Z ihrivnac $

// The nonvirtual public interface reader class for g4tools based analysis.
// It is defined as a composite of object manager base classes.
// Individual use of the component managers is disabled
// (except for file manager and Hn manager which are also used from
//  other object managers).
// The functions which has to be implemented in concrete managers 
// are declared as virtual protected.

// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4VAnalysisReader_h
#define G4VAnalysisReader_h 1

#include "G4AnalysisManagerState.hh"
#include "globals.hh"

#include <vector>
#include <fstream>


class G4HnManager;
class G4VH1Manager;
class G4VH2Manager;
class G4VRNtupleManager;
class G4VFileManager;

class G4VAnalysisReader
{
  public:
    G4VAnalysisReader(const G4String& type, G4bool isMaster);
    virtual ~G4VAnalysisReader();
   
    // Methods for handling files 
    G4bool OpenFile(const G4String& fileName = "");

    // Methods for handling files and directories names  
    G4bool SetFileName(const G4String& fileName);
    G4String GetFileName() const;
    
    // Methods to read histograms from a file
    G4int GetH1(const G4String& h1Name, const G4String& fileName = "");
    G4int GetH2(const G4String& h2Name, const G4String& fileName = "");
                           
    // The ids of histograms and ntuples are generated automatically
    // starting from 0; with following functions it is possible to
    // change the first Id to start from other value
    G4bool SetFirstHistoId(G4int firstId);
    G4bool SetFirstH1Id(G4int firstId);
    G4bool SetFirstH2Id(G4int firstId);
    G4bool SetFirstNtupleId(G4int firstId);

    // Methods to read ntuple from a file
    G4int GetNtuple(const G4String& ntupleName, const G4String& fileName = "");
    
    // Methods for ntuple with id = FirstNtupleId
    G4bool SetNtupleIColumn(const G4String& columnName, G4int& value);
    G4bool SetNtupleFColumn(const G4String& columnName, G4float& value);
    G4bool SetNtupleDColumn(const G4String& columnName, G4double& value);
    // Bind the ntuple colums of vector type
    G4bool SetNtupleIColumn(const G4String& columnName, 
                            std::vector<int>& vector);
    G4bool SetNtupleFColumn(const G4String& columnName, 
                            std::vector<float>& vector);
    G4bool SetNtupleDColumn(const G4String& columnName, 
                            std::vector<double>& vector);
    // Methods for ntuple with id > FirstNtupleId   
    G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            G4int& value);
    G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            G4float& value);
    G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            G4double& value);
    // Bind the ntuple colums of vector type
    G4bool SetNtupleIColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<int>& vector);
    G4bool SetNtupleFColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<float>& vector);
    G4bool SetNtupleDColumn(G4int ntupleId, const G4String& columnName, 
                            std::vector<double>& vector);

    G4bool GetNtupleRow();
    G4bool GetNtupleRow(G4int ntupleId);

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

    virtual G4int  GetH1Impl(const G4String& /*h1Name*/, 
                             const G4String& /*fileName*/) { return false; }
    virtual G4int  GetH2Impl(const G4String& /*h2Name*/, 
                             const G4String& /*fileName*/) { return false; }
    virtual G4int  GetNtupleImpl(const G4String& /*ntupleName*/, 
                             const G4String& /*fileName*/) { return false; }
 
    // methods
    void SetH1Manager(G4VH1Manager* h1Manager);
    void SetH2Manager(G4VH2Manager* h2Manager);
    void SetNtupleManager(G4VRNtupleManager* ntupleManager);
    void SetFileManager(G4VFileManager* fileManager);

    // data members
    G4AnalysisManagerState fState;

  private:
    // data members
    G4HnManager*      fH1HnManager;
    G4HnManager*      fH2HnManager;
    G4VH1Manager*     fVH1Manager;
    G4VH2Manager*     fVH2Manager;
    G4VRNtupleManager* fVNtupleManager;
    G4VFileManager*    fVFileManager;
};

// inline functions

#include "G4VAnalysisReader.icc"
 
#endif

