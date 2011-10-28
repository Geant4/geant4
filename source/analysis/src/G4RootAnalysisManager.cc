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

#include "G4RootAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

G4RootAnalysisManager* G4RootAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4RootAnalysisManager* G4RootAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new G4RootAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4RootAnalysisManager::G4RootAnalysisManager()
 : G4VAnalysisManager("Root"),
   fFile(0),
   fHistoDirectory(0),
   fNtupleDirectory(0),
   fH1Vector(),   
   fH1MapByName(),
   fH2Vector(),   
   fH2MapByName(),
   fNtupleName(),
   fNtupleTitle(),
   fNtuple(0),
   fNtupleIColumnMap(),
   fNtupleFColumnMap(),
   fNtupleDColumnMap()
{
  if ( fgInstance ) {
    G4ExceptionDescription description;
    description << "      " 
                << "G4RootAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("G4RootAnalysisManager::G4RootAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  fgInstance = this;
}

//_____________________________________________________________________________
G4RootAnalysisManager::~G4RootAnalysisManager()
{  
/*
  // CHECK
  std::vector<tools::histo::h1d*>::iterator it;
  for ( it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  }  
  std::vector<tools::histo::h2d*>::iterator it;
  for ( it = fH2Vector.begin(); it != fH2Vector.end(); it++ ) {
    delete *it;
  }  
  delete fNtuple;
*/
  delete fFile;  

  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
tools::wroot::ntuple::column<int>*    
G4RootAnalysisManager::GetNtupleIColumn(G4int id) const
{
  std::map<G4int, tools::wroot::ntuple::column<int>* >::const_iterator it
    = fNtupleIColumnMap.find(id);
  if ( it == fNtupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::GetNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wroot::ntuple::column<float>*  
G4RootAnalysisManager::GetNtupleFColumn(G4int id) const
{
  std::map<G4int, tools::wroot::ntuple::column<float>* >::const_iterator it
    = fNtupleFColumnMap.find(id);
  if ( it == fNtupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wroot::ntuple::column<double>* 
G4RootAnalysisManager::GetNtupleDColumn(G4int id) const
{
  std::map<G4int, tools::wroot::ntuple::column<double>* >::const_iterator it
    = fNtupleDColumnMap.find(id);
  if ( it == fNtupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::GetNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::OpenFile(const G4String& fileName)
{
  // Add file extension .root if no extension is given
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) name.append(".root");
  
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("open", "analysis file", name);
#endif
  
  fFile = new tools::wroot::file(std::cout, name);

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("open", "analysis file", name);
#endif
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::Write() 
{
  // h1 histograms
  std::map<G4String, tools::histo::h1d*>::iterator it;
  for ( it = fH1MapByName.begin(); it != fH1MapByName.end(); it++ ) {
#ifdef G4VERBOSE
    if ( fpVerboseL2 ) 
      fpVerboseL2->Message("write", "h1d", it->first);
#endif
    G4bool result
      = to(*fHistoDirectory,*(it->second),it->first);
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "saving histo " << it->first << " failed";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W003", JustWarning, description);
      return false;       
    } 
  }

  // h2 histograms
  std::map<G4String, tools::histo::h2d*>::iterator it2;
  for ( it2 = fH2MapByName.begin(); it2 != fH2MapByName.end(); it2++ ) {
#ifdef G4VERBOSE
    if ( fpVerboseL2 ) 
      fpVerboseL2->Message("write", "h2d", it2->first);
#endif
    G4bool result
      = to(*fHistoDirectory,*(it2->second),it2->first);
    if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "saving histo " << it2->first << " failed";
      G4Exception("G4RootAnalysisManager::Write()",
                "Analysis_W003", JustWarning, description);
      return false;       
    } 
  }

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("write", "file", "");
#endif

  unsigned int n;
  G4bool result = fFile->write(n);

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("write", "file", "", result);
#endif

  return result;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::CloseFile()
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("close", "file", "");
#endif

  fFile->close();  

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", "");
#endif

  return true;
} 
   
//_____________________________________________________________________________
G4int G4RootAnalysisManager::CreateH1(const G4String& name,  const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax)
{

  if ( ! fHistoDirectory ) {
#ifdef G4VERBOSE
    if ( fpVerboseL2 ) 
      fpVerboseL2->Message("create", "directory for histograms", fHistoDirectoryName);
#endif
  
    fHistoDirectory = fFile->dir().mkdir(fHistoDirectoryName);
    if ( ! fHistoDirectory ) {
      G4ExceptionDescription description;
      description << "      " 
                  << "cannot create directory " << fHistoDirectoryName;
      G4Exception("G4RootAnalysisManager::CreateH1()",
                "Analysis_W002", JustWarning, description);
      return -1;       
    }       
#ifdef G4VERBOSE
    else {
      if ( fpVerboseL1 ) 
        fpVerboseL1->Message("create", "directory for histograms", fHistoDirectoryName);
    }    
#endif
 }

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "H1", name);
#endif
  G4int index = fH1Vector.size();
  tools::histo::h1d* h1 = new tools::histo::h1d(title, nbins, xmin, xmax);
  fH1Vector.push_back(h1);
  fH1MapByName[name] = h1;
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "H1", name);
#endif
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4int G4RootAnalysisManager::CreateH2(const G4String& name,  const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax)
                               
{

  if ( ! fHistoDirectory ) {
#ifdef G4VERBOSE
    if ( fpVerboseL2 ) 
      fpVerboseL2->Message("create", "directory for histograms", fHistoDirectoryName);
#endif
    fHistoDirectory = fFile->dir().mkdir(fHistoDirectoryName);
    if ( ! fHistoDirectory ) {
      G4ExceptionDescription description;
      description << "      " 
                  << "cannot create directory " << fHistoDirectoryName;
      G4Exception("G4RootAnalysisManager::CreateH2()",
                "Analysis_W002", JustWarning, description);
      return -1;       
    }       
#ifdef G4VERBOSE
    else {
      if ( fpVerboseL1 ) 
        fpVerboseL1->Message("create", "directory for histograms", fHistoDirectoryName);
    }    
#endif
  }

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "H2", name);
#endif
  G4int index = fH2Vector.size();
  tools::histo::h2d* h2 
    = new tools::histo::h2d(title, nxbins, xmin, xmax, nybins, ymin, ymax);
  fH2Vector.push_back(h2);
  fH2MapByName[name] = h2;
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "H2", name);
#endif
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
void G4RootAnalysisManager::CreateNtuple(const G4String& name, 
                                         const G4String& title)
{
  if ( fNtuple ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple already exists. "
                << "(Only one ntuple is currently supported.)";
    G4Exception("G4RootAnalysisManager::CreateNtuple()",
                "Analysis_W006", JustWarning, description);
    return;       
  }

  if ( ! fNtupleDirectory ) {
#ifdef G4VERBOSE
    if ( fpVerboseL2 ) 
      fpVerboseL2
        ->Message("create", "directory for ntuples", fNtupleDirectoryName);
#endif
    fNtupleDirectory = fFile->dir().mkdir(fNtupleDirectoryName);
    if ( ! fNtupleDirectory ) {
      G4ExceptionDescription description;
      description << "      " 
                  << "cannot create directory " << fNtupleDirectoryName;
      G4Exception("G4RootAnalysisManager::CreateNtuple()",
                  "Analysis_W002", JustWarning, description);
      return;       
    }       
#ifdef G4VERBOSE
    else {
      if ( fpVerboseL1 ) 
        fpVerboseL1
          ->Message("create", "directory for ntuples", fNtupleDirectoryName);
    }    
#endif
  }

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple", name);
#endif

  fNtuple = new tools::wroot::ntuple(*fNtupleDirectory, name, title);
  fNtupleName = name;
  fNtupleTitle = title;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple", name);
#endif
}                                         

//_____________________________________________________________________________
G4int G4RootAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple I column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wroot::ntuple::column<int>* column = fNtuple->create_column<int>(name);  
  fNtupleIColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple I column", name);
#endif

  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4RootAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple F column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wroot::ntuple::column<float>* column = fNtuple->create_column<float>(name);  
  fNtupleFColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple F column", name);
#endif

  return index + fFirstNtupleId;       
}                                         


//_____________________________________________________________________________
G4int G4RootAnalysisManager::CreateNtupleDColumn(const G4String& name)   
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple D column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wroot::ntuple::column<double>* column = fNtuple->create_column<double>(name);  
  fNtupleDColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple D column", name);
#endif

  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
void G4RootAnalysisManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::FillH1(G4int id, G4double value, G4double weight)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "H1", description);
  }  
#endif

  tools::histo::h1d* h1d = GetH1(id);
  if ( ! h1d ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::FillH1()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  h1d->fill(value, weight);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL1->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::FillH2(G4int id, 
                                     G4double xvalue, G4double yvalue, 
                                     G4double weight)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
    fpVerboseL2->Message("fill", "H2", description);
  }  
#endif

  tools::histo::h2d* h2d = GetH2(id);
  if ( ! h2d ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::FillH2()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  h2d->fill(xvalue, yvalue, weight);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
    fpVerboseL1->Message("fill", "H2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple I column", description);
  }  
#endif

  tools::wroot::ntuple::column<int>* column = GetNtupleIColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::FillNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL1->Message("fill", "ntuple I column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple F column", description);
  }  
#endif

  tools::wroot::ntuple::column<float>* column = GetNtupleFColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::FillNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL1->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4RootAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple D column", description);
  }  
#endif

  tools::wroot::ntuple::column<double>* column = GetNtupleDColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4RootAnalysisManager::FillNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL1->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4RootAnalysisManager::AddNtupleRow()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL2 )
    fpVerboseL2->Message("add", "ntuple row", "");
#endif

  if ( ! fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4RootAnalysisManager::AddNtupleRow()",
                "Analysis_W008", JustWarning, description);
    return false;
  }  
  
  G4bool result =fNtuple->add_row();
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "      " << "adding row has failed.";
    G4Exception("G4RootAnalysisManager::AddNtupleRow()",
                "Analysis_W004", JustWarning, description);
  }         
#ifdef G4VERBOSE
  if ( fpVerboseL1 )
    fpVerboseL1->Message("add", "ntuple row", "", result);
#endif

  return result;
}
 
//_____________________________________________________________________________
tools::histo::h1d*  G4RootAnalysisManager::GetH1(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4RootAnalysisManager::GetH1()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  return fH1Vector[index];
}

//_____________________________________________________________________________
tools::histo::h2d*  G4RootAnalysisManager::GetH2(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4RootAnalysisManager::GetH2()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  return fH2Vector[index];
}


