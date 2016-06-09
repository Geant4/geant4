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

#include "G4CsvAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include "tools/waxml/begend"
#include "tools/waxml/histos"

#include <iostream>

G4CsvAnalysisManager* G4CsvAnalysisManager::fgInstance = 0;

//_____________________________________________________________________________
G4CsvAnalysisManager* G4CsvAnalysisManager::Instance()
{
  if ( fgInstance == 0 ) {
    fgInstance = new G4CsvAnalysisManager();
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4CsvAnalysisManager::G4CsvAnalysisManager()
 : G4VAnalysisManager("Csv"),
   fFile(0),
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
                << "G4CsvAnalysisManager already exists." 
                << "Cannot create another instance.";
    G4Exception("G4CsvAnalysisManager::G4CsvAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }              
   
  fgInstance = this;
}

//_____________________________________________________________________________
G4CsvAnalysisManager::~G4CsvAnalysisManager()
{  
  delete fNtuple;
  delete fFile;

  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
tools::wcsv::ntuple::column<int>*    
G4CsvAnalysisManager::GetNtupleIColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<int>* >::const_iterator it
    = fNtupleIColumnMap.find(id);
  if ( it == fNtupleIColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleIColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
    
//_____________________________________________________________________________
tools::wcsv::ntuple::column<float>*  
G4CsvAnalysisManager::GetNtupleFColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<float>* >::const_iterator it
    = fNtupleFColumnMap.find(id);
  if ( it == fNtupleFColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleFColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  


//_____________________________________________________________________________
tools::wcsv::ntuple::column<double>* 
G4CsvAnalysisManager::GetNtupleDColumn(G4int id) const
{
  std::map<G4int, tools::wcsv::ntuple::column<double>* >::const_iterator it
    = fNtupleDColumnMap.find(id);
  if ( it == fNtupleDColumnMap.end() ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::GetNtupleDColumn()",
                "Analysis_W009", JustWarning, description);
    return 0;
  }
  
  return it->second;
}  
 
// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFile(const G4String& fileName)
{
  // Add file extension .csv if no extension is given
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) name.append(".csv");

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("open", "analysis file", name);
#endif
  
  fFile = new std::ofstream(name);
  if ( fFile->fail() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4CsvAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("open", "analysis file", name);
#endif
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Write() 
{
  // nothing to be done here
  return true;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFile()
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
G4int G4CsvAnalysisManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::CreateH1()",
            "Analysis_W005", JustWarning, description);
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::CreateH2()",
            "Analysis_W005", JustWarning, description);
  return 0;
}                                         

//_____________________________________________________________________________
void G4CsvAnalysisManager::CreateNtuple(const G4String& name, 
                                        const G4String& title)
{
  if ( fNtuple ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple already exists. "
                << "(Only one ntuple is currently supported.)";
    G4Exception("G4CsvAnalysisManager::CreateNtuple()",
                "Analysis_W006", JustWarning, description);
    return;       
  }

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple", name);
#endif

  fNtuple = new tools::wcsv::ntuple(*fFile);
  fNtupleName = name;
  fNtupleTitle = title;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple", name);
#endif
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple I column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<int>* column = fNtuple->create_column<int>(name);  
  fNtupleIColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple I column", name);
#endif

  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple F column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<float>* column = fNtuple->create_column<float>(name);  
  fNtupleFColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple F column", name);
#endif

  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleDColumn(const G4String& name)   
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple D column", name);
#endif

  G4int index = fNtuple->columns().size();
  tools::wcsv::ntuple::column<double>* column = fNtuple->create_column<double>(name);  
  fNtupleDColumnMap[index] = column;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("create", "ntuple D column", name);
#endif

  return index + fFirstNtupleId;       
}                                         

//_____________________________________________________________________________
void G4CsvAnalysisManager::FinishNtuple()
{ 
  // nothing to be done here
}
   
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillH1(G4int /*id*/, 
                                    G4double /*value*/, G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::FillH1()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillH2(G4int /*id*/, 
                                    G4double /*xvalue*/, G4double /*yvalue*/,
                                    G4double /*weight*/)
{
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception("G4CsvAnalysisManager::FillH2()",
            "Analysis_W007", JustWarning, description);
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple I column", description);
  }  
#endif

  tools::wcsv::ntuple::column<int>* column = GetNtupleIColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleIColumn()",
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
G4bool G4CsvAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple F column", description);
  }  
#endif

  tools::wcsv::ntuple::column<float>* column = GetNtupleFColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleFColumn()",
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
G4bool G4CsvAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL2->Message("fill", "ntuple D column", description);
  }  
#endif

   tools::wcsv::ntuple::column<double>* column = GetNtupleDColumn(id);
  if ( ! column ) {
    G4ExceptionDescription description;
    description << "      " << "column " << id << " does not exist.";
    G4Exception("G4CsvAnalysisManager::FillNtupleDColumn()",
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
G4bool G4CsvAnalysisManager::AddNtupleRow()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL2 )
    fpVerboseL2->Message("add", "ntuple row", "");
#endif

  if ( ! fNtuple ) {
    G4ExceptionDescription description;
    description << "      " << "ntuple does not exist. ";
    G4Exception("G4CsvAnalysisManager::AddNtupleRow()",
                "Analysis_W008", JustWarning, description);
    return false;
  }  
  
  fNtuple->add_row();
#ifdef G4VERBOSE
  if ( fpVerboseL1 )
    fpVerboseL1->Message("add", "ntuple row", "");
#endif

  return true;
}
