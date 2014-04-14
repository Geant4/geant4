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
// $Id: G4VAnalysisReader.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#include "G4VAnalysisReader.hh"
#include "G4AnalysisUtilities.hh"
#include "G4HnManager.hh"
#include "G4VH1Manager.hh"
#include "G4VH2Manager.hh"
#include "G4VRNtupleManager.hh"
#include "G4VFileManager.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4VAnalysisReader::G4VAnalysisReader(const G4String& type, G4bool isMaster)
 : fState(type, isMaster),
   fH1HnManager(0),
   fH2HnManager(0),
   fVH1Manager(0),
   fVH2Manager(0),
   fVNtupleManager(0),
   fVFileManager(0)
{}

//_____________________________________________________________________________
G4VAnalysisReader::~G4VAnalysisReader()
{
  delete fVH1Manager;
  delete fVH2Manager;
  delete fVNtupleManager;
  delete fVFileManager;
}

// 
// protected methods
//
    
//_____________________________________________________________________________
void G4VAnalysisReader::SetH1Manager(G4VH1Manager* h1Manager)
{
  fVH1Manager = h1Manager;
  fH1HnManager = h1Manager->fHnManager;
} 

//_____________________________________________________________________________
void G4VAnalysisReader::SetH2Manager(G4VH2Manager* h2Manager)
{
  fVH2Manager = h2Manager;
  fH2HnManager = h2Manager->fHnManager;
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetNtupleManager(G4VRNtupleManager* ntupleManager)
{
  fVNtupleManager = ntupleManager;
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetFileManager(G4VFileManager* fileManager)
{
  fVFileManager = fileManager;
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4VAnalysisReader::OpenFile(const G4String& fileName)
{
  if ( fileName != "" ) {
    return OpenFileImpl(fileName);
  }
  else {  
    if ( fVFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot open file. File name is not defined.";
      G4Exception("G4VFileManager::OpenFile()",
                  "Analysis_W009", JustWarning, description);
      return false;
    }           
    return OpenFileImpl(fVFileManager->GetFileName());
  }  
} 

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFileName(const G4String& fileName)
{ 
  return fVFileManager->SetFileName(fileName); 
}

//_____________________________________________________________________________
G4String G4VAnalysisReader::GetFileName() const 
{  
  return fVFileManager->GetFileName(); 
}

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetH1(const G4String& h1Name, 
                               const G4String& fileName)
{
  if ( fileName != "" ) {
    return GetH1Impl(h1Name, fileName);
  }
  else {  
    if ( fVFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get H1. File name has to be set first.";
      G4Exception("G4VAnalysisReader::GetH1()",
                  "Analysis_WR001", JustWarning, description);
      return -1;
    }           
    return GetH1Impl(h1Name, fVFileManager->GetFileName());
  }  
}                                      

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetH2(const G4String& h2Name, 
                               const G4String& fileName)
{
  if ( fileName != "" ) {
    return GetH2Impl(h2Name, fileName);
  }
  else {  
    if ( fVFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get H2. File name has to be set first.";
      G4Exception("G4VAnalysisReader::GetH2()",
                  "Analysis_WR001", JustWarning, description);
      return -1;
    }           
    return GetH2Impl(h2Name, fVFileManager->GetFileName());
  }  
}                                      

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstHistoId(G4int firstId) 
{
  G4bool finalResult = true;
  G4bool result = fVH1Manager->SetFirstId(firstId);
  finalResult = finalResult && result;
  
  result = fH1HnManager->SetFirstId(firstId);
  finalResult = finalResult && result;
  
  result = fVH2Manager->SetFirstId(firstId);
  finalResult = finalResult && result;

  result = fH2HnManager->SetFirstId(firstId);
  finalResult = finalResult && result;
   
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstH1Id(G4int firstId) 
{
  G4bool finalResult = true;
  G4bool result = fVH1Manager->SetFirstId(firstId);
  finalResult = finalResult && result;
  
  result = fH1HnManager->SetFirstId(firstId);
  finalResult = finalResult && result;
  
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstH2Id(G4int firstId) 
{
  G4bool finalResult = true;
  G4bool result = fVH2Manager->SetFirstId(firstId);
  finalResult = finalResult && result;

  result = fH2HnManager->SetFirstId(firstId);
  finalResult = finalResult && result;
   
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstNtupleId(G4int firstId) 
{
  return fVNtupleManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNtuple(const G4String& ntupleName, 
                                   const G4String& fileName)
{
  if ( fileName != "" ) {
    return GetNtupleImpl(ntupleName, fileName);
  }
  else {  
    if ( fVFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get Ntuple. File name has to be set first.";
      G4Exception("G4VAnalysisReader::GetNtuple()",
                  "Analysis_WR001", JustWarning, description);
      return -1;
    }           
    return GetNtupleImpl(ntupleName, fVFileManager->GetFileName());
  }  
}                                      

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleIColumn(const G4String& columnName, 
                                            G4int& value)
{
  return fVNtupleManager->SetNtupleIColumn(columnName, value);
}                                              

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleFColumn(const G4String& columnName, 
                                            G4float& value)
{
  return fVNtupleManager->SetNtupleFColumn(columnName, value);
}                                              
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleDColumn(const G4String& columnName, 
                                            G4double& value)
{
  return fVNtupleManager->SetNtupleDColumn(columnName, value);
}                                              

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleIColumn(const G4String& columnName, 
                                            std::vector<int>& vector)
{
  return fVNtupleManager->SetNtupleIColumn(columnName, vector);
}                                              
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleFColumn(const G4String& columnName, 
                                            std::vector<float>& vector)
{
  return fVNtupleManager->SetNtupleFColumn(columnName, vector);
}                                              
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleDColumn(const G4String& columnName, 
                                            std::vector<double>& vector)
{
  return fVNtupleManager->SetNtupleDColumn(columnName, vector);
}                                              
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleIColumn(G4int ntupleId,
                                            const G4String& columnName, 
                                            G4int& value)
{
  return fVNtupleManager->SetNtupleIColumn(ntupleId, columnName, value);
}                                              

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleFColumn(G4int ntupleId,
                                            const G4String& columnName, 
                                            G4float& value)
{
  return fVNtupleManager->SetNtupleFColumn(ntupleId, columnName, value);
}                                              
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleDColumn(G4int ntupleId,
                                            const G4String& columnName, 
                                            G4double& value)
{
  return fVNtupleManager->SetNtupleDColumn(ntupleId, columnName, value);
}                                              

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleIColumn(G4int ntupleId, 
                                            const G4String& columnName, 
                                            std::vector<int>& vector)
{
  return fVNtupleManager->SetNtupleIColumn(ntupleId, columnName, vector);
}
                                            
//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleFColumn(G4int ntupleId, 
                                            const G4String& columnName, 
                                            std::vector<float>& vector)
{
  return fVNtupleManager->SetNtupleFColumn(ntupleId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetNtupleDColumn(G4int ntupleId, 
                                            const G4String& columnName, 
                                            std::vector<double>& vector)
{
  return fVNtupleManager->SetNtupleDColumn(ntupleId, columnName, vector);
}

//_____________________________________________________________________________
G4bool G4VAnalysisReader::GetNtupleRow()
{
  return fVNtupleManager->GetNtupleRow();
}                                              


//_____________________________________________________________________________
G4bool G4VAnalysisReader::GetNtupleRow(G4int ntupleId)
{
  return fVNtupleManager->GetNtupleRow(ntupleId);
}                                              

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofH1s() const
{
  return fH1HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofH2s() const
{
  return fH2HnManager->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofNtuples() const
{
  return fVNtupleManager->GetNofNtuples();
}  

// GetH1Id(), GetH2Id in .icc 

// Access methods in .icc

//_____________________________________________________________________________
void G4VAnalysisReader::SetVerboseLevel(G4int verboseLevel) 
{
  fState.SetVerboseLevel(verboseLevel);
} 

// GetVerboseLevel() in .icc 
