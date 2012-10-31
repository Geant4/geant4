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

#include "G4VAnalysisManager.hh"
#include "G4AnalysisMessenger.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4VAnalysisManager::G4VAnalysisManager(const G4String& type)
  : fVerboseLevel(0),
    fActivation(false),
    fFirstHistoId(0),
    fFirstNtupleColumnId(0),
    fFileName(""), 
    fHistoDirectoryName(""), 
    fNtupleDirectoryName(""),
    fLockFirstHistoId(false),
    fLockFirstNtupleColumnId(false),
    fLockFileName(false),
    fLockHistoDirectoryName(false), 
    fLockNtupleDirectoryName(false),
    fVerboseL1(type,1),
    fVerboseL2(type,2),
    fVerboseL3(type,3),
    fVerboseL4(type,4),
    fpVerboseL1(0),
    fpVerboseL2(0),
    fpVerboseL3(0),
    fpVerboseL4(0),
    fMessenger(0),
    fNofActiveObjects(0),
    fNofAsciiObjects(0),
    fH1Informations(),
    fH2Informations()
{
  fMessenger = new G4AnalysisMessenger(this);
}

//_____________________________________________________________________________
G4VAnalysisManager::~G4VAnalysisManager()
{
  delete fMessenger;
  // add delete G4HnInformation objects
}

// 
// protected methods
//

//_____________________________________________________________________________
void G4VAnalysisManager::AddH1Information(const G4String& name,
                                          const G4String& unitName,
                                          const G4String& fcnName,
                                          G4double unit, G4Fcn fcn)
{
  fH1Informations.push_back(
    new G4HnInformation(name, unitName, unitName, fcnName, fcnName,
                        unit, unit, fcn, fcn));
  ++fNofActiveObjects;
}  

//_____________________________________________________________________________
void  G4VAnalysisManager::AddH2Information(const G4String& name,
                                           const G4String& xunitName, 
                                           const G4String& yunitName,
                                           const G4String& xfcnName,
                                           const G4String& yfcnName,
                                           G4double xunit, G4double yunit,
                                           G4Fcn xfcn, G4Fcn yfcn) 
{
  fH2Informations.push_back(
    new G4HnInformation(name, xunitName, yunitName, xfcnName, yfcnName,
                        xunit, yunit, xfcn, yfcn));
  ++fNofActiveObjects;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::WriteAscii()
{
  // Replace or add file extension .ascii
  G4String name(fFileName);
  if ( name.find(".") != std::string::npos ) { 
    name.erase(name.find("."), name.length()); 
  }
  name.append(".ascii");

#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("write ASCII", "file", name);
#endif
     
  std::ofstream output(name, std::ios::out);
  if ( ! output ) {
    G4ExceptionDescription description;
    description 
      << "Cannot open file. File name is not defined.";
    G4Exception("G4VAnalysisManager::WriteAscii()",
                "Analysis_W009", JustWarning, description);
    return false;
  }
  output.setf( std::ios::scientific, std::ios::floatfield );

  G4bool result = WriteOnAscii(output);

#ifdef G4VERBOSE
    if ( fpVerboseL1 ) 
      fpVerboseL1->Message("write ASCII", "file",  name, result);
#endif
  
  return result;
}     

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetName(ObjectType type, G4int id) const
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return "";
    
  return info->fName;
}    

//_____________________________________________________________________________
G4double G4VAnalysisManager::GetXUnit(ObjectType type, G4int id) const
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return 1.0;
  
  return info->fXUnit;
}    

//_____________________________________________________________________________
G4double G4VAnalysisManager::GetYUnit(ObjectType type, G4int id) const
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return 1.0;
  
  return info->fYUnit;
}    

//_____________________________________________________________________________
G4bool G4VAnalysisManager::GetActivation(ObjectType type, G4int id) const
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return true;
  
  return info->fActivation;
}    

//_____________________________________________________________________________
G4bool G4VAnalysisManager::GetAscii(ObjectType type, G4int id) const
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return false;
  
  return info->fAscii;
}    

//_____________________________________________________________________________
G4double G4VAnalysisManager::GetUnitValue(const G4String& unit) const
{
   G4double value = 1.;
   if ( unit != "none" ) {
     value = G4UnitDefinition::GetValueOf(unit);
     if ( value == 0. ) value = 1.; 
   }  
   return value;
}   

//_____________________________________________________________________________
G4Fcn G4VAnalysisManager::GetFunction(const G4String& fcnName) const
{
  G4Fcn fcn = G4FcnIdentity;
   if ( fcnName != "none" ) {
    if      ( fcnName == "log" )  fcn = std::log;
    else if ( fcnName == "log10") fcn = std::log10;
    else if ( fcnName == "exp" )  fcn = std::exp;
    else {
      G4ExceptionDescription description;
      description 
        << "    \"" << fcnName << "\" function is not supported." << G4endl
        << "    " << "No function will be applied to h1 values.";
      G4Exception("G4AnalysisMessenger::GetFunction",
                "Analysis_W013", JustWarning, description);
    }              
  }
  return fcn;            
}
    
// 
// public methods
//

//_____________________________________________________________________________
void G4VAnalysisManager::SetVerboseLevel(G4int verboseLevel) 
{
  if ( verboseLevel == fVerboseLevel || verboseLevel < 0 ) return;
  
  fVerboseLevel = verboseLevel;
  
  if ( verboseLevel == 0 ) {
    fpVerboseL1 = 0;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 1 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 2 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 3 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = &fVerboseL3;
    fpVerboseL4 = 0;
  }
  else {
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = &fVerboseL3;
    fpVerboseL4 = &fVerboseL4;
  }
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::OpenFile()
{
  if ( fFileName == "" ) {
    G4ExceptionDescription description;
    description 
      << "Cannot open file. File name is not defined.";
    G4Exception("G4VAnalysisManager::OpenFile()",
                "Analysis_W009", JustWarning, description);
    return false;
  }           
  
  return OpenFile(fFileName);
}     

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFileName(const G4String& fileName) 
{
  if ( fLockFileName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set File name as its value was already used.";
    G4Exception("G4VAnalysisManager::SetFileName()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fFileName = fileName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetHistoDirectoryName(const G4String& dirName) 
{
  if ( fLockHistoDirectoryName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Histo directory name as its value was already used.";
    G4Exception("G4VAnalysisManager::SetHistoDirectoryName()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fHistoDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetNtupleDirectoryName(const G4String& dirName) 
{
  if ( fLockNtupleDirectoryName ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set Ntuple directory name as its value was already used.";
    G4Exception("G4VAnalysisManager::SetNtupleDirectoryName()",
                "Analysis_W010", JustWarning, description);
    return false;
  }              

  fNtupleDirectoryName = dirName;
  return true;
}  

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetFullFileName() const 
{  
  G4String name(fFileName);
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }  

  return name;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstHistoId(G4int firstId) 
{
  if ( fLockFirstHistoId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstHistoId as its value was already used.";
    G4Exception("G4VAnalysisManager::SetFirstHistoId()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fFirstHistoId = firstId;
  return true;
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::SetFirstNtupleColumnId(G4int firstId) 
{
  if ( fLockFirstNtupleColumnId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstNtupleColumnId as its value was already used.";
    G4Exception("G4VAnalysisManager::SetFirstHistoId()",
                "Analysis_W010", JustWarning, description);
    return false;
  }              

  fFirstNtupleColumnId = firstId;
  return true;
}

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsActive() const
{
  if ( ! fActivation ) return true;
  
  return ( fNofActiveObjects > 0 );
}  

//_____________________________________________________________________________
G4bool G4VAnalysisManager::IsAscii() const
{
  return ( fNofAsciiObjects > 0 );
}  

//_____________________________________________________________________________
G4HnInformation* G4VAnalysisManager::GetH1Information(G4int id) const
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= GetNofH1s() ) {
    G4ExceptionDescription description;
    description << "      " << "histo " << id << " does not exist.";
    G4Exception("G4VAnalysisManager::GetH1Information()",
                "Analysis_W007", JustWarning, description);
    return 0;         
  }
  return fH1Informations[index];
}    

//_____________________________________________________________________________
G4HnInformation* G4VAnalysisManager::GetH2Information(G4int id) const
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= GetNofH2s() ) {
    G4ExceptionDescription description;
    description << "      " << "histo " << id << " does not exist.";
    G4Exception("G4VAnalysisManager::GetH2Information()",
                "Analysis_W007", JustWarning, description);
    return 0;         
  }
  return fH2Informations[index];
}    
    
//_____________________________________________________________________________
G4HnInformation* G4VAnalysisManager::GetInformation(ObjectType objType, G4int id) const
{
  switch ( objType ) {
    case kH1: 
      return GetH1Information(id);
      break;
      
    case kH2: 
      return GetH2Information(id);
      break;
      
    case kNtuple:
    default:
      return 0;
      break;
  }    

  // Cannot reach this line
  G4ExceptionDescription description;
  description << "Wrong object type.";
  G4Exception("G4VAnalysisManager::SetFirstHistoId()",
              "Analysis_W010", FatalException, description);
  return 0;
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetActivation(ObjectType type, G4int id, 
                                        G4bool activation)
{
// Set activation to a given object

  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return;

  // Do nothing if activation does not change
  if ( info->fActivation == activation ) return;
  
  // Change activation and account it in fNofActiveObjects
  info->fActivation = activation;
  if ( activation ) 
    fNofActiveObjects++;
  else
    fNofActiveObjects--;   
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetActivation(ObjectType type, G4bool activation)
{
// Set activation to all objects of the given type

  std::vector<G4HnInformation*>* informations;
  if ( type == kH1 ) 
    informations = &fH1Informations;
  else if ( type == kH2 )
    informations = &fH2Informations;
  else  if ( type == kNtuple ) {
    return;
  }
  else {  
    G4ExceptionDescription description;
    description << "Wrong object type.";
    G4Exception("G4VAnalysisManager::SetActivation()",
              "Analysis_W010", FatalException, description);
    return;
  }   
  
  std::vector<G4HnInformation*>::iterator it;
  for ( it = informations->begin(); it != informations->end(); it++ ) {
    G4HnInformation* info = *it;

    // Do nothing if activation does not change
    if ( info->fActivation == activation ) continue;
  
    // Change activation and account it in fNofActiveObjects
    info->fActivation = activation;
    if ( activation ) 
      fNofActiveObjects++;
    else
      fNofActiveObjects--; 
  }     
}    

//_____________________________________________________________________________
void  G4VAnalysisManager::SetAscii(ObjectType type, G4int id, G4bool ascii)
{
  G4HnInformation* info = GetInformation(type, id);

  if ( ! info ) return;

  // Do nothing if ascii does not change
  if ( info->fAscii == ascii ) return;
  
  // Change ascii and account it in fNofAsciiObjects
  info->fAscii = ascii;
  if ( ascii ) 
    fNofAsciiObjects++;
  else
    fNofAsciiObjects--;   
}    

//_____________________________________________________________________________
G4String G4VAnalysisManager::GetFileType() const 
{
  G4String fileType = fVerboseL1.GetType();
  fileType.toLower();
  return fileType;
}                 
