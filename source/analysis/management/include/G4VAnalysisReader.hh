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
#include <memory>

class G4HnManager;
class G4BaseFileManager;
class G4VH1Manager;
class G4VH2Manager;
class G4VH3Manager;
class G4VP1Manager;
class G4VP2Manager;
class G4VRNtupleManager;

class G4VAnalysisReader
{
  public:
    G4VAnalysisReader(const G4String& type, G4bool isMaster);
    virtual ~G4VAnalysisReader();
   
    // Methods for handling files 
    // G4bool OpenFile(const G4String& fileName = "");

    // Methods for handling files and directories names  
    void SetFileName(const G4String& fileName);
    G4String GetFileName() const;
    
    // Methods to read histograms from a file
    G4int ReadH1(const G4String& h1Name, const G4String& fileName = "", const G4String& dirName = "");
    G4int ReadH2(const G4String& h2Name, const G4String& fileName = "", const G4String& dirName = "");
    G4int ReadH3(const G4String& h3Name, const G4String& fileName = "", const G4String& dirName = "");
    G4int ReadP1(const G4String& h1Name, const G4String& fileName = "", const G4String& dirName = "");
    G4int ReadP2(const G4String& h2Name, const G4String& fileName = "", const G4String& dirName = "");
                           
    // The ids of histograms and ntuples are generated automatically
    // starting from 0; with following functions it is possible to
    // change the first Id to start from other value
    G4bool SetFirstHistoId(G4int firstId);
    G4bool SetFirstH1Id(G4int firstId);
    G4bool SetFirstH2Id(G4int firstId);
    G4bool SetFirstH3Id(G4int firstId);
    G4bool SetFirstProfileId(G4int firstId);
    G4bool SetFirstP1Id(G4int firstId);
    G4bool SetFirstP2Id(G4int firstId);
    G4bool SetFirstNtupleId(G4int firstId);

    // Methods to read ntuple from a file
    G4int GetNtuple(const G4String& ntupleName, const G4String& fileName = "", const G4String& dirName = "");
    
    // Methods for ntuple with id = FirstNtupleId
    G4bool SetNtupleIColumn(const G4String& columnName, G4int& value);
    G4bool SetNtupleFColumn(const G4String& columnName, G4float& value);
    G4bool SetNtupleDColumn(const G4String& columnName, G4double& value);
    G4bool SetNtupleSColumn(const G4String& columnName, G4String& value);
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
    G4bool SetNtupleSColumn(G4int ntupleId, const G4String& columnName, 
                            G4String& value);
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
    G4int GetNofH3s() const;
    G4int GetNofP1s() const;
    G4int GetNofP2s() const;
    G4int GetNofNtuples() const;

    // Access methods via names
    G4int GetH1Id(const G4String& name, G4bool warn = true) const;
    G4int GetH2Id(const G4String& name, G4bool warn = true) const;
    G4int GetH3Id(const G4String& name, G4bool warn = true) const;
    G4int GetP1Id(const G4String& name, G4bool warn = true) const;
    G4int GetP2Id(const G4String& name, G4bool warn = true) const;
    
    // Access to histogram & profiles parameters
    //
    G4int    GetH1Nbins(G4int id) const;
    G4double GetH1Xmin(G4int id) const;
    G4double GetH1Xmax(G4int id) const;
    G4double GetH1Width(G4int id) const;
    //
    G4int    GetH2Nxbins(G4int id) const;
    G4double GetH2Xmin(G4int id) const;
    G4double GetH2Xmax(G4int id) const;
    G4double GetH2XWidth(G4int id) const;
    G4int    GetH2Nybins(G4int id) const;
    G4double GetH2Ymin(G4int id) const;
    G4double GetH2Ymax(G4int id) const;
    G4double GetH2YWidth(G4int id) const;
    //
    G4int    GetH3Nxbins(G4int id) const;
    G4double GetH3Xmin(G4int id) const;
    G4double GetH3Xmax(G4int id) const;
    G4double GetH3XWidth(G4int id) const;
    G4int    GetH3Nybins(G4int id) const;
    G4double GetH3Ymin(G4int id) const;
    G4double GetH3Ymax(G4int id) const;
    G4double GetH3YWidth(G4int id) const;
    G4int    GetH3Nzbins(G4int id) const;
    G4double GetH3Zmin(G4int id) const;
    G4double GetH3Zmax(G4int id) const;
    G4double GetH3ZWidth(G4int id) const;
    //
    G4int    GetP1Nbins(G4int id) const;
    G4double GetP1Xmin(G4int id) const;
    G4double GetP1Xmax(G4int id) const;
    G4double GetP1XWidth(G4int id) const;
    G4double GetP1Ymin(G4int id) const;
    G4double GetP1Ymax(G4int id) const;
    //
    G4int    GetP2Nxbins(G4int id) const;
    G4double GetP2Xmin(G4int id) const;
    G4double GetP2Xmax(G4int id) const;
    G4double GetP2XWidth(G4int id) const;
    G4int    GetP2Nybins(G4int id) const;
    G4double GetP2Ymin(G4int id) const;
    G4double GetP2Ymax(G4int id) const;
    G4double GetP2YWidth(G4int id) const;
    G4double GetP2Zmin(G4int id) const;
    G4double GetP2Zmax(G4int id) const;

    // Access histogram & profiles attributes for plotting
    //
    G4String GetH1Title(G4int id) const;
    G4String GetH1XAxisTitle(G4int id) const;
    G4String GetH1YAxisTitle(G4int id) const;
    //
    G4String GetH2Title(G4int id) const;
    G4String GetH2XAxisTitle(G4int id) const;
    G4String GetH2YAxisTitle(G4int id) const;
    G4String GetH2ZAxisTitle(G4int id) const;
    //
    G4String GetH3Title(G4int id) const;
    G4String GetH3XAxisTitle(G4int id) const;
    G4String GetH3YAxisTitle(G4int id) const;
    G4String GetH3ZAxisTitle(G4int id) const;
    //
    G4String GetP1Title(G4int id) const;
    G4String GetP1XAxisTitle(G4int id) const;
    G4String GetP1YAxisTitle(G4int id) const;
    G4String GetP1ZAxisTitle(G4int id) const;
    //
    G4String GetP2Title(G4int id) const;
    G4String GetP2XAxisTitle(G4int id) const;
    G4String GetP2YAxisTitle(G4int id) const;
    G4String GetP2ZAxisTitle(G4int id) const;

    // Verbosity
    void  SetVerboseLevel(G4int verboseLevel);
    G4int GetVerboseLevel() const;

    // The manager type (starts with an uppercase letter)
    G4String GetType() const;                 
    // The manager file type (starts with a lowercase letter)
    G4String GetFileType() const;                 
   
  protected:
    // virtual methods
    virtual G4int  ReadH1Impl(const G4String& h1Name, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
    virtual G4int  ReadH2Impl(const G4String& h2Name, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
    virtual G4int  ReadH3Impl(const G4String& h3Name, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
    virtual G4int  ReadP1Impl(const G4String& p1Name, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
    virtual G4int  ReadP2Impl(const G4String& p2Name, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
    virtual G4int  ReadNtupleImpl(const G4String& ntupleName, const G4String& fileName,
                              const G4String& dirName, G4bool isUserFileName) = 0;
 
    // methods
    void SetH1Manager(G4VH1Manager* h1Manager);
    void SetH2Manager(G4VH2Manager* h2Manager);
    void SetH3Manager(G4VH3Manager* h3Manager);
    void SetP1Manager(G4VP1Manager* p1Manager);
    void SetP2Manager(G4VP2Manager* p2Manager);
    void SetNtupleManager(G4VRNtupleManager* ntupleManager);
    void SetFileManager(G4BaseFileManager* fileManager);

    // data members
    G4AnalysisManagerState fState;

  private:
    // data members
    std::unique_ptr<G4VH1Manager>     fVH1Manager;
    std::unique_ptr<G4VH2Manager>     fVH2Manager;
    std::unique_ptr<G4VH3Manager>     fVH3Manager;
    std::unique_ptr<G4VP1Manager>     fVP1Manager;
    std::unique_ptr<G4VP2Manager>     fVP2Manager;
    std::unique_ptr<G4VRNtupleManager> fVNtupleManager;
    std::unique_ptr<G4BaseFileManager> fFileManager;
};

// inline functions

#include "G4VAnalysisReader.icc"
 
#endif

