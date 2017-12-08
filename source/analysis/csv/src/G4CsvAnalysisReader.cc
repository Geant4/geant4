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
// $Id: G4CsvAnalysisReader.cc 74257 2013-10-02 14:24:55Z gcosmo $

// Author: Ivana Hrivnacova, 05/09/2014 (ivana@ipno.in2p3.fr)

#include "G4CsvAnalysisReader.hh"
#include "G4CsvRFileManager.hh"
#include "G4CsvRNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Threading.hh"

#include <tools/aida_ntuple>
#include <tools/rcsv_histo>

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

G4CsvAnalysisReader* G4CsvAnalysisReader::fgMasterInstance = nullptr;
G4ThreadLocal G4CsvAnalysisReader* G4CsvAnalysisReader::fgInstance = nullptr;

//
// utility functions
//

namespace {

//_____________________________________________________________________________
void*  ReadObject(std::istream& hnFile,
                  const G4String& objectType,
                  const G4String& fileName,
                  const G4String& inFunction)
{
  tools::rcsv::histo handler(hnFile);
  std::string objectTypeInFile;
  void* object;
  auto verbose = false;
  if ( ! handler.read(G4cout, objectTypeInFile, object, verbose) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Cannot get "<< objectType << " in file " << fileName; 
    G4String inFunctionFull = "G4CsvAnalysisReader::";
    inFunctionFull.append(inFunction);
    G4Exception(inFunctionFull, "Analysis_WRnullptr11", JustWarning, description);
    return nullptr;
  }
  if ( objectTypeInFile != objectType ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Object type read in "<< fileName
      << " does not match" << G4endl; 
    G4String inFunctionFull = "G4CsvAnalysisReader::";
    inFunctionFull.append(inFunction);
    G4Exception(inFunctionFull, "Analysis_WR011", JustWarning, description);
    return nullptr;
  }
  
  return object;
}

}

//_____________________________________________________________________________
G4CsvAnalysisReader* G4CsvAnalysisReader::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4CsvAnalysisReader(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4CsvAnalysisReader::G4CsvAnalysisReader(G4bool isMaster)
 : G4ToolsAnalysisReader("Csv", isMaster),
   fNtupleManager(nullptr),
   fFileManager(nullptr)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4CsvAnalysisReader already exists." 
      << "Cannot create another instance.";
    G4Exception("G4CsvAnalysisReader::G4CsvAnalysisReader()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // Create managers
  fNtupleManager = new G4CsvRNtupleManager(fState);
  fFileManager = new G4CsvRFileManager(fState);
      // The managers will be deleted by the base class
  
  // Set managers to base class
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4CsvAnalysisReader::~G4CsvAnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
G4String G4CsvAnalysisReader::GetHnFileName(
                                 const G4String& hnType, 
                                 const G4String& hnName,
                                 const G4String& fileName,
                                 G4bool isUserFileName) const
{
  if ( isUserFileName ) {
    return fFileManager->GetFullFileName(fileName);
  }
  else {  
    return fFileManager->GetHnFileName(hnType, hnName);
  }
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisReader::Reset()
{
// Reset histograms and ntuple

  auto finalResult = true;
  
  auto result = G4ToolsAnalysisReader::Reset();
  finalResult = finalResult && result;
  
  result = fNtupleManager->Reset();
  finalResult = finalResult && result;
  
  return finalResult;
}  
 
// 
// protected methods
//

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadH1Impl(const G4String& h1Name, 
                                      const G4String& fileName,
                                      const G4String& /*dirName*/, 
                                      G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("get", "h1", h1Name);
#endif

  // open file
  auto h1FileName = GetHnFileName("h1", h1Name, fileName, isUserFileName);
  std::ifstream hnFile(h1FileName);
  if ( ! hnFile.is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << h1FileName;
    G4Exception("G4CsvAnalysisReader::ReadH1Impl()",
                "Analysis_WR001", JustWarning, description);
    return kInvalidId;
  }
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read file", h1FileName);
#endif

  void* object 
    = ReadObject(hnFile, tools::histo::h1d::s_class(), h1FileName, "ReadH1Impl");
  if ( ! object ) return kInvalidId;
        
  auto h1 = static_cast<tools::histo::h1d*>(object);
  auto id = fH1Manager->AddH1(h1Name, h1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h1", h1Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadH2Impl(const G4String& h2Name, 
                                      const G4String& fileName,
                                      const G4String& /*dirName*/, 
                                      G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h2", h2Name);
#endif

  // open file
  auto h2FileName = GetHnFileName("h2", h2Name, fileName, isUserFileName);
  std::ifstream hnFile(h2FileName);
  if ( ! hnFile.is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << h2FileName;
    G4Exception("G4CsvAnalysisReader::ReadH2Impl()",
                "Analysis_WR001", JustWarning, description);
    return kInvalidId;
  }
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read file", h2FileName);
#endif

  void* object 
    = ReadObject(hnFile, tools::histo::h2d::s_class(), h2FileName, "ReadH2Impl");
  if ( ! object ) return kInvalidId;
        
  auto h2 = static_cast<tools::histo::h2d*>(object);
  auto id = fH2Manager->AddH2(h2Name, h2);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h2", h2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadH3Impl(const G4String& h3Name, 
                                      const G4String& fileName,
                                      const G4String& /*dirName*/, 
                                      G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h3", h3Name);
#endif

  // open file
  auto h3FileName = GetHnFileName("h3", h3Name, fileName, isUserFileName);
  std::ifstream hnFile(h3FileName);
  if ( ! hnFile.is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << h3FileName;
    G4Exception("G4CsvAnalysisReader::ReadH3Impl()",
                "Analysis_WR001", JustWarning, description);
    return kInvalidId;
  }
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read file", h3FileName);
#endif

  void* object 
    = ReadObject(hnFile, tools::histo::h3d::s_class(), h3FileName, "ReadH3Impl");
  if ( ! object ) return kInvalidId;
        
  auto h3 = static_cast<tools::histo::h3d*>(object);
  auto id = fH3Manager->AddH3(h3Name, h3);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h3", h3Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadP1Impl(const G4String& p1Name, 
                                      const G4String& fileName,
                                      const G4String& /*dirName*/, 
                                      G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p1", p1Name);
#endif

  // open file
  G4String p1FileName = GetHnFileName("p1", p1Name, fileName, isUserFileName);
  std::ifstream hnFile(p1FileName);
  if ( ! hnFile.is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << p1FileName;
    G4Exception("G4CsvAnalysisReader::ReadP1Impl()",
                "Analysis_WR001", JustWarning, description);
    return kInvalidId;
  }
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read file", p1FileName);
#endif

  void* object 
    = ReadObject(hnFile, tools::histo::p1d::s_class(), fileName, "ReadP1Impl");
  if ( ! object ) return kInvalidId;
        
  auto p1 = static_cast<tools::histo::p1d*>(object);
  auto id = fP1Manager->AddP1(p1Name, p1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p1", p1Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadP2Impl(const G4String& p2Name, 
                                      const G4String& fileName,
                                      const G4String& /*dirName*/, 
                                      G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p2", p2Name);
#endif

  // open file
  G4String p2FileName = GetHnFileName("p2", p2Name, fileName, isUserFileName);
  std::ifstream hnFile(p2FileName);
  if ( ! hnFile.is_open() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << p2FileName;
    G4Exception("G4CsvAnalysisReader::ReadP2Impl()",
                "Analysis_WR001", JustWarning, description);
    return kInvalidId;
  }
#ifdef G4VERBOSE
  if ( fState.GetVerboseL1() ) 
    fState.GetVerboseL1()
      ->Message("open", "read file", p2FileName);
#endif

  void* object 
    = ReadObject(hnFile, tools::histo::p2d::s_class(), p2FileName, "ReadP2Impl");
  if ( ! object ) return kInvalidId;
        
  auto p2 = static_cast<tools::histo::p2d*>(object);
  auto id = fP2Manager->AddP2(p2Name, p2);


#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p2", p2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4CsvAnalysisReader::ReadNtupleImpl(const G4String& ntupleName, 
                                          const G4String& fileName,
                                          const G4String& /*dirName*/, 
                                          G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "ntuple", ntupleName);
#endif

  // Ntuples are saved per object and per thread
  // but apply the ntuple name and the thread suffixes
  // only if fileName is not provided explicitly
  G4String fullFileName = fileName;
  if ( ! isUserFileName ) {
    fullFileName = fFileManager->GetNtupleFileName(ntupleName);
  }  

  // Open file
  if ( ! fFileManager->OpenRFile(fullFileName) ) return kInvalidId;
  auto ntupleFile = fFileManager->GetRFile(fullFileName);
  
  // Create ntuple 
  auto rntuple = new tools::rcsv::ntuple(*ntupleFile);
  auto id = fNtupleManager->SetNtuple(new G4TRNtupleDescription<tools::rcsv::ntuple>(rntuple));
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "ntuple", ntupleName, id > kInvalidId);
#endif
  
  return id;
}  
