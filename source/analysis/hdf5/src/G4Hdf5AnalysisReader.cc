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

// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#include "G4Hdf5AnalysisReader.hh"
#include "G4Hdf5RNtupleManager.hh"
#include "G4AnalysisVerbose.hh"
#include "G4AnalysisUtilities.hh"
#include "G4Threading.hh"

#include <iostream>
#include <cstdio>

using namespace G4Analysis;

G4Hdf5AnalysisReader* G4Hdf5AnalysisReader::fgMasterInstance = nullptr;
G4ThreadLocal G4Hdf5AnalysisReader* G4Hdf5AnalysisReader::fgInstance = nullptr;

//_____________________________________________________________________________
G4Hdf5AnalysisReader* G4Hdf5AnalysisReader::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4Hdf5AnalysisReader(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4Hdf5AnalysisReader::G4Hdf5AnalysisReader(G4bool isMaster)
 : G4ToolsAnalysisReader("Hdf5", isMaster),
   fNtupleManager(nullptr),
   fFileManager(nullptr)
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4Hdf5AnalysisReader already exists." 
      << "Cannot create another instance.";
    G4Exception("G4Hdf5AnalysisReader::G4Hdf5AnalysisReader()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;

  // Create managers
  fNtupleManager = new G4Hdf5RNtupleManager(fState);
  fFileManager = new G4Hdf5RFileManager(fState);
  
  // Set managers to base class
  SetNtupleManager(fNtupleManager);
  SetFileManager(fFileManager);
}

//_____________________________________________________________________________
G4Hdf5AnalysisReader::~G4Hdf5AnalysisReader()
{
  if ( fState.GetIsMaster() ) fgMasterInstance = nullptr;
  fgInstance = nullptr;
}

// 
// private methods
//

//_____________________________________________________________________________
G4bool G4Hdf5AnalysisReader::Reset()
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
G4int G4Hdf5AnalysisReader::ReadH1Impl(const G4String& h1Name, 
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h1", h1Name);
#endif

  auto h1 = ReadHnImpl<tools::histo::h1d>(h1Name, fileName, dirName);
 
  if ( ! h1 ) return kInvalidId;
  
  auto id = fH1Manager->AddH1(h1Name, h1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h1", h1Name, id > kInvalidId);
#endif
  
  return id;
}  

//_____________________________________________________________________________
G4int G4Hdf5AnalysisReader::ReadH2Impl(const G4String& h2Name, 
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h2", h2Name);
#endif

  auto h2 = ReadHnImpl<tools::histo::h2d>(h2Name, fileName, dirName);
 
  if ( ! h2 ) return kInvalidId;
  
  auto id = fH2Manager->AddH2(h2Name, h2);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h2", h2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4Hdf5AnalysisReader::ReadH3Impl(const G4String& h3Name, 
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool /*isUserFileName*/)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "h3", h3Name);
#endif

  auto h3 = ReadHnImpl<tools::histo::h3d>(h3Name, fileName, dirName);
 
  if ( ! h3 ) return kInvalidId;
    
  auto id = fH3Manager->AddH3(h3Name, h3);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "h3", h3Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4Hdf5AnalysisReader::ReadP1Impl(const G4String& p1Name, 
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool /*isUserFileName*/)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p1", p1Name);
#endif

  auto p1 = ReadPnImpl<tools::histo::p1d>(p1Name, fileName, dirName);
 
  if ( ! p1 ) return kInvalidId;
  
  auto id = fP1Manager->AddP1(p1Name, p1);

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p1", p1Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4Hdf5AnalysisReader::ReadP2Impl(const G4String& p2Name, 
                                       const G4String& fileName,
                                       const G4String& dirName,
                                       G4bool /*isUserFileName*/)
{

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "p2", p2Name);
#endif

  auto p2 = ReadPnImpl<tools::histo::p2d>(p2Name, fileName, dirName);
 
  if ( ! p2 ) return kInvalidId;
  
  auto id = fP2Manager->AddP2(p2Name, p2);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "p2", p2Name, id > kInvalidId);
#endif
  
  return id;  
}  

//_____________________________________________________________________________
G4int G4Hdf5AnalysisReader::ReadNtupleImpl(const G4String& ntupleName, 
                                           const G4String& fileName,
                                           const G4String& dirName,
                                           G4bool isUserFileName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("read", "ntuple", ntupleName);
#endif

  // Ntuples are saved in files per thread
  // but apply thethe thread suffix only if fileName is not provided explicitly
  G4String fullFileName = fileName;
  if ( ! isUserFileName ) {
    fullFileName = fFileManager->GetFullFileName();
  }  

  // Get directory
  auto directory = fFileManager->GetNtupleRDirectory(fullFileName, dirName, false);
  if ( directory < 0 ) return kInvalidId;

  // Create ntuple 
  auto rntuple = new tools::hdf5::ntuple(G4cout, directory, ntupleName);
  auto rntupleDescription = new G4TRNtupleDescription<tools::hdf5::ntuple>(rntuple);
  auto id = fNtupleManager->SetNtuple(rntupleDescription);
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("read", "ntuple", ntupleName, id > kInvalidId);
#endif

  return id;
}  
