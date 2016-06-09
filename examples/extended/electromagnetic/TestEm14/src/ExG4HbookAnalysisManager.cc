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
/// \file electromagnetic/TestEm14/src/ExG4HbookAnalysisManager.cc
/// \brief Implementation of the ExG4HbookAnalysisManager class
//
// $Id$
//
/// \file ExG4HbookAnalysisManager.cc
/// \brief Implementation of the ExG4HbookAnalysisManager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

extern "C" int setpawc();
extern "C" int setntuc();

ExG4HbookAnalysisManager* ExG4HbookAnalysisManager::fgInstance = 0;
const G4int ExG4HbookAnalysisManager::fgkDefaultH2HbookIdOffset = 100;
const G4int ExG4HbookAnalysisManager::fgkDefaultNtupleHbookId = 1;
const G4String ExG4HbookAnalysisManager::fgkDefaultNtupleDirectoryName = "ntuple";

//_____________________________________________________________________________
ExG4HbookAnalysisManager* ExG4HbookAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new ExG4HbookAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
ExG4HbookAnalysisManager::ExG4HbookAnalysisManager()
 : G4VAnalysisManager("Hbook"),
   fH1HbookIdOffset(-1),
   fH2HbookIdOffset(-1),
   fNtupleHbookId(-1),
   fFile(0),
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
                << "G4HbookAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("G4HbookAnalysisManager::G4HbookAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  fgInstance = this;
  
  // Initialize HBOOK :
  tools::hbook::CHLIMIT(setpawc());
  setntuc(); //for ntuple.
}

//_____________________________________________________________________________
ExG4HbookAnalysisManager::~ExG4HbookAnalysisManager()
{  
  std::vector<tools::hbook::h1*>::iterator it;
  for ( it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  }  
  std::vector<tools::hbook::h2*>::iterator it2;
  for ( it2 = fH2Vector.begin(); it2 != fH2Vector.end(); it2++ ) {
    delete *it2;
  }  
  delete fNtuple;
  delete fFile;  

  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
tools::hbook::wntuple::column<int>*    
ExG4HbookAnalysisManager::GetNtupleIColumn(G4int id) const
{
  std::map<G4int, tools::hbook::wntuple::column<int>* >::const_iterator it
    = fNtupleIColumnMap.find(id);
  if ( it == fNtupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::hbook::wntuple::column<float>*  
ExG4HbookAnalysisManager::GetNtupleFColumn(G4int id) const
{
  std::map<G4int, tools::hbook::wntuple::column<float>* >::const_iterator it
    = fNtupleFColumnMap.find(id);
  if ( it == fNtupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::hbook::wntuple::column<double>* 
ExG4HbookAnalysisManager::GetNtupleDColumn(G4int id) const
{
  std::map<G4int, tools::hbook::wntuple::column<double>* >::const_iterator it
    = fNtupleDColumnMap.find(id);
  if ( it == fNtupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::OpenFile(const G4String& fileName)
{
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }  
  
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("open", "analysis file", name);
#endif
  
  tools::hbook::CHCDIR("//PAWC"," ");
  
  unsigned int unit = 1;
  fFile = new tools::hbook::wfile(std::cout, name, unit);
  if ( ! fFile->is_valid() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4HbookAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;       
  }

  // At this point, in HBOOK, we should have :
  //   - created a //LUN1 directory attached to the file
  //   - created a //PAWC/LUN1 in memory
  //   - be in the directory //PAWC/LUN1.

  // create an "histo" HBOOK directory both in memory and in the file :
  if ( fHistoDirectoryName != "" ) {
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    tools::hbook::CHMDIR(fHistoDirectoryName.data()," ");
    tools::hbook::CHCDIR("//LUN1"," ");
    tools::hbook::CHMDIR(fHistoDirectoryName.data()," ");
  }
  
  fLockHistoDirectoryName = true;

  // the five upper lines could have been done with :
  //fFile->cd_home();
  //fFile->mkcd("histo");

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("open", "analysis file", name);
#endif
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::Write() 
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("write", "file", "");
#endif

  // ntuple 
  //if ( fNtuple ) fNtuple->add_row_end();

  G4bool result = fFile->write();  

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("write", "file", "", result);
#endif

  return result;  
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::CloseFile()
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("close", "file", "");
#endif

  //WARNING : have to delete the ntuple before closing the file.
  delete fNtuple;
  fNtuple = 0;

  G4bool result = fFile->close();  

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", "", result);
#endif

  return result;
} 
   
//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateH1(const G4String& name, const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "H1", name);
#endif

  // Go to histograms directory
  if ( fHistoDirectoryName != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fHistoDirectoryName.data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  

  // Set  fH1HbookIdOffset if needed
  if (  fH1Vector.size() == 0 ) {
    if ( fH1HbookIdOffset == -1 ) {
      if ( fFirstHistoId > 0 ) 
        fH1HbookIdOffset = 0;
      else
        fH1HbookIdOffset = 1;
          
      if ( fH1HbookIdOffset > 0 ) {
        G4ExceptionDescription description;
        description << "H1 will be defined in HBOOK with ID = G4_firstHistoId + 1";
        G4Exception("ExG4HbookAnalysisManager::CreateH1()",
                    "Analysis_W011", JustWarning, description);
      }              
    }
  }  
  
  // Create histogram    
  G4int index = fH1Vector.size();
  G4int hbookIndex = fH1HbookIdOffset + fH1Vector.size() + fFirstHistoId;
  tools::hbook::h1* h1 = new tools::hbook::h1(hbookIndex, title, nbins, xmin, xmax);
  fH1Vector.push_back(h1);
  fH1MapByName[name] = h1;
 
  if ( fHistoDirectoryName != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
  
  fLockFirstHistoId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) { 
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << hbookIndex; 
    fpVerboseL1->Message("create", "H1", description);
  }  
#endif

  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateH2(const G4String& name, const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "H2", name);
#endif

  // Go to histograms directory
  if ( fHistoDirectoryName != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fHistoDirectoryName.data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  

  // Set fH2HbookIdOffset if needed
  if (  fH2Vector.size() == 0 ) {
    if ( fH2HbookIdOffset == -1 ) { 
      if ( fFirstHistoId > 0 ) 
        fH2HbookIdOffset = fgkDefaultH2HbookIdOffset;
      else
        fH2HbookIdOffset = fgkDefaultH2HbookIdOffset + 1;
           
      if ( fH2HbookIdOffset != fgkDefaultH2HbookIdOffset ) {
        G4ExceptionDescription description;
        description 
          << "H2 will be defined in HBOOK with ID = " 
          << fgkDefaultH2HbookIdOffset << " + G4_firstHistoId + 1";
        G4Exception("ExG4HbookAnalysisManager::CreateH2()",
                    "Analysis_W011", JustWarning, description);
      }
    }                
  }    

  G4int index = fH2Vector.size();
  G4int hbookIndex = fH2HbookIdOffset + fH2Vector.size() + fFirstHistoId;
  tools::hbook::h2* h2 
    = new tools::hbook::h2(hbookIndex, title, nxbins, xmin, xmax, nybins, ymin, ymax);
  fH2Vector.push_back(h2);
  fH2MapByName[name] = h2;
 
  // Return to //PAWC/LUN1 :
  if ( fHistoDirectoryName != "" ) {
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  

  fLockFirstHistoId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << hbookIndex; 
    fpVerboseL1->Message("create", "H2", description);
  }  
#endif

  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::CreateNtuple(const G4String& name, 
                                          const G4String& title)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "ntuple", name);
#endif

  if ( fNtuple ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple already exists. "
                << "(Only one ntuple is currently supported.)";
    G4Exception("G4HbookAnalysisManager::CreateNtuple()",
                "Analysis_W006", JustWarning, description);
    return;       
  }

  // Create an "ntuple" directory both in memory and in the file
  fFile->cd_home();      //go under //PAWC/LUN1
  if ( fNtupleDirectoryName == "" )
    fFile->mkcd(fgkDefaultNtupleDirectoryName.data());
  else  
    fFile->mkcd(fNtupleDirectoryName.data());

  // Define ntuple ID in HBOOK
  if ( fNtupleHbookId == -1 ) fNtupleHbookId = fgkDefaultNtupleHbookId;
  
  // We should be under //PAWC/LUN1/ntuple
  fNtuple = new tools::hbook::wntuple(fNtupleHbookId, name);
  fNtupleName = name;
  fNtupleTitle = title;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) {
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << fNtupleHbookId; 
    fpVerboseL1->Message("create", "ntuple", description);
  }  
#endif
}                                         

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "ntuple I column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::hbook::wntuple::column<int>* column = fNtuple->create_column<int>(name);  
  fNtupleIColumnMap[index] = column;
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple I column", name);
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "ntuple F column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::hbook::wntuple::column<float>* column = fNtuple->create_column<float>(name);  
  fNtupleFColumnMap[index] = column;
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple F column", name);
#endif

  return index + fFirstNtupleColumnId;       
}                                         


//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateNtupleDColumn(const G4String& name) 
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create", "ntuple D column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::hbook::wntuple::column<double>* column = fNtuple->create_column<double>(name);  
  fNtupleDColumnMap[index] = column;
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple D column", name);
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::FinishNtuple()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("finish", "ntuple", fNtupleName);
#endif

  // Return to //PAWC/LUN1 :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");

  //fNtuple->add_row_beg();
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("finish", "ntuple", fNtupleName);
#endif
}
  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillH1(G4int id, G4double value, G4double weight)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL3->Message("fill", "H1", description);
  }  
#endif

  tools::hbook::h1* h1 = GetH1(id);
  if ( ! h1 ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillH1()",
                "Analysis_W007", JustWarning, description);
    return false;
  }  

  h1->fill(value, weight);
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillH2(G4int id, 
                                       G4double xvalue, G4double yvalue,
                                       G4double weight)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
    fpVerboseL3->Message("fill", "H2", description);
  }  
#endif

  tools::hbook::h2* h2 = GetH2(id);
  if ( ! h2 ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillH2()",
                "Analysis_W007", JustWarning, description);
    return false;
  }  

  h2->fill(xvalue, yvalue, weight);
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
    fpVerboseL2->Message("fill", "H2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL3->Message("fill", "ntuple I column", description);
  }  
#endif

  tools::hbook::wntuple::column<int>* column = GetNtupleIColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
 #ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple I column", description);
  }  
#endif
 return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL3->Message("fill", "ntuple F column", description);
  }  
#endif

  tools::hbook::wntuple::column<float>* column = GetNtupleFColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL3 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL3->Message("fill", "ntuple D column", description);
  }  
#endif

  tools::hbook::wntuple::column<double>* column = GetNtupleDColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::FillNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return false;
  }  
  
  column->fill(value);
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::AddNtupleRow()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL3 )
    fpVerboseL3->Message("add", "ntuple row", "");
#endif

  //G4cout << "Hbook: Going to add Ntuple row ..." << G4endl;
  if ( ! fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4HbookAnalysisManager::AddNtupleRow()",
                "Analysis_W008", JustWarning, description);
    return false;
  }  
  
  //fNtuple->add_row_fast();
  fNtuple->add_row();
#ifdef G4VERBOSE
  if ( fpVerboseL2 )
    fpVerboseL2->Message("add", "ntuple row", "");
#endif
  return true;
}
 
//_____________________________________________________________________________
tools::hbook::h1*  ExG4HbookAnalysisManager::GetH1(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH1()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  return fH1Vector[index];
}

//_____________________________________________________________________________
tools::hbook::h2*  ExG4HbookAnalysisManager::GetH2(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH2()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  return fH2Vector[index];
}

//_____________________________________________________________________________
tools::hbook::wntuple* ExG4HbookAnalysisManager::GetNtuple() const
{
  return fNtuple;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH1HbookIdOffset(G4int offset) 
{
  if ( fH1Vector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set H1HbookIdOffset as some H1 histogramms already exist.";
    G4Exception("G4HbookAnalysisManager::SetH1HbookIdOffset()",
                 "Analysis_W009", JustWarning, description);
    return false;             
  }
  
  if ( fFirstHistoId + offset < 1 ) {
    G4ExceptionDescription description;
    description << "The first histogram HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetH1HbookIdOffset()",
                 "Analysis_W009", JustWarning, description);
    return false;             
  }
  
  fH1HbookIdOffset = offset;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2HbookIdOffset(G4int offset) 
{
  if ( fH2Vector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set H2HbookIdOffset as some H2 histogramms already exist.";
    G4Exception("G4HbookAnalysisManager::SetH2HbookIdOffset()",
                 "Analysis_W009", JustWarning, description);
    return false;             
  }

  if ( fFirstHistoId + offset < 1 ) {
    G4ExceptionDescription description;
    description << "The first histogram HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetH1HbookIdOffset()",
                 "Analysis_W009", JustWarning, description);
    return false;             
  }
  
  fH2HbookIdOffset = offset;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetNtupleHbookId(G4int ntupleId)
{
  if ( fNtuple ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set NtupleHbookId as an ntuple already exists.";
    G4Exception("G4HbookAnalysisManager::SetNtupleHbookId()",
                 "Analysis_W010", JustWarning, description);
    return false;             
  }

  if ( ntupleId < 1 ) {
    G4ExceptionDescription description;
    description << "The ntuple HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetNtupleHbookId()",
                 "Analysis_W010", JustWarning, description);
    return false;             
  }
  
  fNtupleHbookId = ntupleId;
  return true;
}

#endif
