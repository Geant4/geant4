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
#include "G4BaseFileManager.hh"
#include "G4HnManager.hh"
#include "G4VH1Manager.hh"
#include "G4VH2Manager.hh"
#include "G4VH3Manager.hh"
#include "G4VP1Manager.hh"
#include "G4VP2Manager.hh"
#include "G4VRNtupleManager.hh"

#include <iostream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4VAnalysisReader::G4VAnalysisReader(const G4String& type, G4bool isMaster)
 : fState(type, isMaster),
   fVH1Manager(nullptr),
   fVH2Manager(nullptr),
   fVH3Manager(nullptr),
   fVP1Manager(nullptr),
   fVP2Manager(nullptr),
   fVNtupleManager(nullptr),
   fFileManager(nullptr)
{}

//_____________________________________________________________________________
G4VAnalysisReader::~G4VAnalysisReader()
{}

// 
// protected methods
//
    
//_____________________________________________________________________________
void G4VAnalysisReader::SetH1Manager(G4VH1Manager* h1Manager)
{
  fVH1Manager.reset(h1Manager);
} 

//_____________________________________________________________________________
void G4VAnalysisReader::SetH2Manager(G4VH2Manager* h2Manager)
{
  fVH2Manager.reset(h2Manager);
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetH3Manager(G4VH3Manager* h3Manager)
{
  fVH3Manager.reset(h3Manager);
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetP1Manager(G4VP1Manager* p1Manager)
{
  fVP1Manager.reset(p1Manager);
} 

//_____________________________________________________________________________
void G4VAnalysisReader::SetP2Manager(G4VP2Manager* p2Manager)
{
  fVP2Manager.reset(p2Manager);
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetNtupleManager(G4VRNtupleManager* ntupleManager)
{
  fVNtupleManager.reset(ntupleManager);
}  

//_____________________________________________________________________________
void G4VAnalysisReader::SetFileManager(G4BaseFileManager* fileManager)
{
  fFileManager.reset(fileManager);
}  

// 
// public methods
//

//_____________________________________________________________________________
void G4VAnalysisReader::SetFileName(const G4String& fileName)
{ 
  fFileManager->SetFileName(fileName); 
}

//_____________________________________________________________________________
G4String G4VAnalysisReader::GetFileName() const 
{  
  return fFileManager->GetFileName(); 
}

//_____________________________________________________________________________
G4int G4VAnalysisReader::ReadH1(const G4String& h1Name, 
                                const G4String& fileName,
                                const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadH1Impl(h1Name, fileName, dirName, true);
  }
  else {  
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get H1. File name has to be set first.";
      G4Exception("G4VAnalysisReader::ReadH1()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadH1Impl(h1Name, fFileManager->GetFileName(), dirName, false);
  }  
}                                      

//_____________________________________________________________________________
G4int G4VAnalysisReader::ReadH2(const G4String& h2Name, 
                                const G4String& fileName,
                                const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadH2Impl(h2Name, fileName, dirName, true);
  }
  else {  
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get H2. File name has to be set first.";
      G4Exception("G4VAnalysisReader::ReadH2()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadH2Impl(h2Name, fFileManager->GetFileName(), dirName, false);
  }  
}                                      

//_____________________________________________________________________________
G4int G4VAnalysisReader::ReadH3(const G4String& h3Name, 
                                const G4String& fileName,
                                const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadH3Impl(h3Name, fileName, dirName, true);
  }
  else {  
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get H3. File name has to be set first.";
      G4Exception("G4VAnalysisReader::ReadH2()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadH3Impl(h3Name, fFileManager->GetFileName(), dirName, false);
  }  
}                                      

//_____________________________________________________________________________
G4int G4VAnalysisReader::ReadP1(const G4String& p1Name, 
                                const G4String& fileName,
                                const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadP1Impl(p1Name, fileName, dirName, true);
  }
  else {  
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get P1. File name has to be set first.";
      G4Exception("G4VAnalysisReader::ReadP1()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadP1Impl(p1Name, fFileManager->GetFileName(), dirName, false);
  }  
}                                      

//_____________________________________________________________________________
G4int G4VAnalysisReader::ReadP2(const G4String& p2Name, 
                                const G4String& fileName,
                                const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadP2Impl(p2Name, fileName, dirName, true);
  }
  else {  
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get P2. File name has to be set first.";
      G4Exception("G4VAnalysisReader::ReadP2()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadP2Impl(p2Name, fFileManager->GetFileName(), dirName, false);
  }  
}                                      

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstHistoId(G4int firstId) 
{
  G4bool finalResult = true;

  G4bool result = SetFirstH1Id(firstId);
  finalResult = finalResult && result;
  
  result = SetFirstH2Id(firstId);
  finalResult = finalResult && result;

  result = SetFirstH3Id(firstId);
  finalResult = finalResult && result;
   
  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstH1Id(G4int firstId) 
{
  return fVH1Manager->GetHnManager()->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstH2Id(G4int firstId) 
{
  return fVH2Manager->GetHnManager()->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstH3Id(G4int firstId) 
{
  return fVH3Manager->GetHnManager()->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstProfileId(G4int firstId) 
{
  G4bool finalResult = true;

  G4bool result = SetFirstP1Id(firstId);
  finalResult = finalResult && result;
  
  result = SetFirstP2Id(firstId);
  finalResult = finalResult && result;

  return finalResult; 
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstP1Id(G4int firstId) 
{
  return fVP1Manager->GetHnManager()->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstP2Id(G4int firstId) 
{
  return fVP2Manager->GetHnManager()->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4bool G4VAnalysisReader::SetFirstNtupleId(G4int firstId) 
{
  return fVNtupleManager->SetFirstId(firstId);
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNtuple(const G4String& ntupleName, 
                                   const G4String& fileName,
                                   const G4String& dirName)
{
  if ( fileName != "" ) {
    return ReadNtupleImpl(ntupleName, fileName, dirName, true);
  }
  else { 
    // Check if fileName was set
    if ( fFileManager->GetFileName() == "" ) {
      G4ExceptionDescription description;
      description 
        << "Cannot get Ntuple. File name has to be set first.";
      G4Exception("G4VAnalysisReader::GetNtuple()",
                  "Analysis_WR011", JustWarning, description);
      return kInvalidId;
    }           
    return ReadNtupleImpl(ntupleName, fFileManager->GetFileName(), dirName, false);
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
G4bool G4VAnalysisReader::SetNtupleSColumn(const G4String& columnName, 
                                            G4String& value)
{
  return fVNtupleManager->SetNtupleSColumn(columnName, value);
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
G4bool G4VAnalysisReader::SetNtupleSColumn(G4int ntupleId,
                                            const G4String& columnName, 
                                            G4String& value)
{
  return fVNtupleManager->SetNtupleSColumn(ntupleId, columnName, value);
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
  return fVH1Manager->GetHnManager()->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofH2s() const
{
  return fVH2Manager->GetHnManager()->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofH3s() const
{
  return fVH3Manager->GetHnManager()->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofP1s() const
{
  return fVP1Manager->GetHnManager()->GetNofHns();
}  

//_____________________________________________________________________________
G4int G4VAnalysisReader::GetNofP2s() const
{
  return fVP2Manager->GetHnManager()->GetNofHns();
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
