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
//
// Author: Ivana Hrivnacova, 30/10/2018  (ivana@ipno.in2p3.fr)
// ---------------------------------------------------------------------

#include "G4ScoreNtupleWriter.hh"
#include "G4ScoreNtupleWriterMessenger.hh"
#include "G4THitsMap.hh"
#include "g4csv.hh"
#include "g4root.hh"
#include "g4xml.hh"
#include "G4Threading.hh"

//_____________________________________________________________________________
//
// ctor, dtor
//

//_____________________________________________________________________________
G4ScoreNtupleWriter::G4ScoreNtupleWriter(const G4String& outputType)
 : G4VScoreNtupleWriter(),
   fMessenger(nullptr),
   fHCIds(),
   fAnalysisManager(nullptr),
   fOutputTypeName(outputType),
   fFileName("scoring"),
   fVerboseLevel(1),
   fHasAnalysisManager(false),
   fHasAnalysisFile(false),
   fIsBooked(false),
   fIsInitialized(false),
   fFirstNtupleId(0)
{
  // if ( G4Threading::IsMasterThread() ) {
    fMessenger = new G4ScoreNtupleWriterMessenger(this);
  // }
}

//_____________________________________________________________________________
G4ScoreNtupleWriter::~G4ScoreNtupleWriter()
{
  if ( fHasAnalysisManager ) {
    delete fAnalysisManager;
  }

  delete fMessenger;
}

//
// private methods
//

//_____________________________________________________________________________
void G4ScoreNtupleWriter::CreateAnalysisManager()
{
  if ( fAnalysisManager ) return;

  switch ( GetOutputType() ) {
    case G4AnalysisOutput::kCsv:
      // create Csv analysis manager
      if ( ! G4CsvAnalysisManager::IsInstance() ) fHasAnalysisManager = true;
      fAnalysisManager = G4CsvAnalysisManager::Instance();
      break;
    case G4AnalysisOutput::kRoot:
      // create Root analysis manager 
      if ( ! G4RootAnalysisManager::IsInstance() ) fHasAnalysisManager = true;
      fAnalysisManager = G4RootAnalysisManager::Instance();
      break;
    case G4AnalysisOutput::kXml:
      // create Xml analysis manager
      if ( ! G4XmlAnalysisManager::IsInstance() ) fHasAnalysisManager = true;
      fAnalysisManager = G4XmlAnalysisManager::Instance();
      break;
    case G4AnalysisOutput::kNone:
      // do nothing
      break;  
  }   

  fAnalysisManager->SetVerboseLevel(fVerboseLevel);
}

//
// protected methods
//

//_____________________________________________________________________________
G4VScoreNtupleWriter* G4ScoreNtupleWriter::CreateInstance() const
{
  auto instance = new G4ScoreNtupleWriter();
  instance->SetOutputType(fOutputTypeName);
  instance->SetFileName(fFileName);
  instance->SetVerboseLevel(fVerboseLevel);

  return instance;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4ScoreNtupleWriter::Book(G4HCofThisEvent* hce)
{
#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- G4ScoreNtupleWriter::Book" << G4endl;
  }
#endif
  if ( ! hce) return false;

  // book collection ID if DoubleHitsMap
  if (fHCIds.size() == 0) {
#ifdef G4VERBOSE
    if ( fVerboseLevel > 1 ) {
        G4cout << "--- going to fill fHCIds " << G4endl;
    }
#endif
    for (G4int i=0; i<hce->GetNumberOfCollections(); ++i) {
      auto hitsMap 
        = dynamic_cast<G4THitsMap<G4double>*>(hce->GetHC(i));

      if ( hitsMap ) {
#ifdef G4VERBOSE
    if ( fVerboseLevel > 1 ) {
       G4cout << "--- adding hcID = " << i << G4endl; 
    }
#endif
       fHCIds.push_back(i);
      }
    }
  }

  // Create analysis manager
  if ( fHCIds.size() ) {
    CreateAnalysisManager();
  }

  // Create ntuples (only once)
  if ( ! fIsInitialized ) {
#ifdef G4VERBOSE
    if ( fVerboseLevel > 1 ) {
      G4cout << "-- going to create ntuples " << G4endl;
    }
#endif
    auto first = true; 
    for (auto id : fHCIds) {
      auto hitsMap 
        = static_cast<G4THitsMap<G4double>*>(hce->GetHC(id));

      // Create ntuple for this primitive
      G4String ntupleName(hitsMap->GetSDname());
      ntupleName.append("_");
      ntupleName.append(hitsMap->GetName());

      G4int ntupleId 
        = fAnalysisManager->CreateNtuple(ntupleName, ntupleName);
      if ( first ) {
        fFirstNtupleId = ntupleId;
#ifdef G4VERBOSE
        if ( fVerboseLevel > 1 ) {
          G4cout << "-- set first ntuple Id " <<  fFirstNtupleId << G4endl;
        }
#endif
        first = false;
      }

      G4String colName(ntupleName);
      colName.append("_eventId");
      fAnalysisManager->CreateNtupleIColumn(colName);
      colName = ntupleName;
      colName.append("_cell");
      fAnalysisManager->CreateNtupleIColumn(colName);
      colName = ntupleName;
      colName.append("_score");
      fAnalysisManager->CreateNtupleDColumn(colName);
      fAnalysisManager->FinishNtuple();
      fIsBooked = true;    
    }
    fIsInitialized = true;
  }

  return fIsBooked;
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::OpenFile()
{
#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
      G4cout << "--- G4ScoreNtupleWriter::OpenFile" << G4endl;
  }
#endif

  if ( fAnalysisManager->IsOpenFile() )  return;

#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
      G4cout << "--- G4ScoreNtupleWriter::OpenFile executing" << G4endl;
  }
#endif

  if ( fAnalysisManager->GetFileName() == "" ) {
    fAnalysisManager->SetFileName(fFileName);
  }
  fAnalysisManager->OpenFile();

#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- G4ScoreNtupleWriter::OpenFile isOpenFile: "
           << fAnalysisManager->IsOpenFile() << G4endl;
  }
#endif

  fHasAnalysisFile = fAnalysisManager->IsOpenFile();
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::Fill(G4HCofThisEvent* hce, G4int eventNumber)
{

#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- start G4ScoreNtupleWriter::Fill" << G4endl;
  }
#endif

  auto counter = 0;
  for (auto id : fHCIds) {
#ifdef G4VERBOSE
    if ( fVerboseLevel > 1 ) {
      G4cout << "in loop over fHCIds, counter " << counter << G4endl;
  }
#endif
    auto hitsMap 
      = static_cast<G4THitsMap<G4double>*>(hce->GetHC(id));

    //G4cout << eventNumber << ".. go to fill ntuple " << counter + fFirstNtupleId << G4endl;

    // fill hits in ntuple
    std::map<G4int, G4double*>::iterator it;
    for ( it = hitsMap->GetMap()->begin(); it != hitsMap->GetMap()->end(); it++ ) {
      fAnalysisManager->FillNtupleIColumn(counter + fFirstNtupleId, 0, eventNumber);    
      fAnalysisManager->FillNtupleIColumn(counter + fFirstNtupleId, 1, it->first);    
      fAnalysisManager->FillNtupleDColumn(counter + fFirstNtupleId, 2, *(it->second)); 
      fAnalysisManager->AddNtupleRow(counter + fFirstNtupleId);
    }
    counter++;
  }

#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- done G4ScoreNtupleWriter::Fill" << G4endl;
  }
#endif
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::Write()
{
#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- start G4ScoreNtupleWriter::Write" << G4endl;
  }
#endif

  if ( fHasAnalysisFile ) {
#ifdef G4VERBOSE
    if ( fVerboseLevel > 1 ) {
      G4cout << "--- sG4ScoreNtupleWriter::Write - has file" << G4endl;
    }
#endif
    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
    fAnalysisManager->SetFileName("");
  }

#ifdef G4VERBOSE
  if ( fVerboseLevel > 1 ) {
    G4cout << "--- done G4ScoreNtupleWriter::Write" << G4endl;
  }
#endif
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::SetOutputType(const G4String& outputTypeName)
{ 
  if ( fAnalysisManager ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set output type. The analysis manager has not been created.";
    G4Exception("G4HitsNtupleWriter::SetOutputType", "", JustWarning, description);
    return;
  }

  fOutputTypeName = outputTypeName; 
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::SetVerboseLevel(G4int value)
{ 
  fVerboseLevel = value; 
  if ( fAnalysisManager ) {
    fAnalysisManager->SetVerboseLevel(value);
  }
}

//_____________________________________________________________________________
void G4ScoreNtupleWriter::SetFileName(const G4String& fileName)
{ 
  fFileName = fileName; 
}

//_____________________________________________________________________________
G4AnalysisOutput G4ScoreNtupleWriter::GetOutputType() const
{ 
  return G4Analysis::GetOutput(fOutputTypeName); 
}

