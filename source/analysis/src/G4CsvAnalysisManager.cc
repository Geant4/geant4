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
   fNtuple(0),
   fNtupleBooking(0),
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
  delete fNtupleBooking;
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
 
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Reset()
{
  delete fNtuple;
  fNtuple = 0;
  
  return true;
}  
 
//_____________________________________________________________________________
void  G4CsvAnalysisManager::ExceptionForHistograms(
                                        const G4String& functionName) const
{
  G4String inFunction = "G4CsvAnalysisManager::";
  inFunction += functionName;
  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;
  G4Exception(inFunction, "Analysis_W005", JustWarning, description);
}  

//
// protected methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// To be added: ntuple

  return true;
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::OpenFile(const G4String& fileName)
{
  // Keep file name
  fFileName =  fileName;

  // Add file extension .csv if no extension is given
  G4String name(fileName);
  if ( name.find(".") == std::string::npos ) { 
    name.append(".");
    name.append(GetFileType());
  }  

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("open", "analysis file", name);
#endif
  
  // delete a previous file if it exists
  if ( fFile ) delete fFile; 
  
  fFile = new std::ofstream(name);
  if ( fFile->fail() ) {
    G4ExceptionDescription description;
    description << "      " << "Cannot open file " << fileName;
    G4Exception("G4CsvAnalysisManager::OpenFile()",
                "Analysis_W001", JustWarning, description);
    return false;
  }

  // Create ntuple if it was already booked
  if ( fNtupleBooking && ( ! fNtuple ) ) {
#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("create from booking", "ntuple", name);
#endif
    fNtuple = new tools::wcsv::ntuple(*fFile, G4cerr, *fNtupleBooking);
    if ( fNtupleBooking->m_columns.size() ) {
      // store ntuple columns in local maps
      const std::vector<tools::ntuple_booking::col_t>& columns 
        = fNtupleBooking->m_columns;
      std::vector<tools::ntuple_booking::col_t>::const_iterator it;
      G4int index = 0;
      for ( it = columns.begin(); it!=columns.end(); ++it) {
        if ( (*it).second == tools::_cid(int(0) ) ) {
          G4cout << "adding int " << fNtuple->find_column<int>((*it).first) << G4endl;
          fNtupleIColumnMap[index++] = fNtuple->find_column<int>((*it).first);
        }
        else if( (*it).second == tools::_cid(float(0) ) ) {
          fNtupleFColumnMap[index++] = fNtuple->find_column<float>((*it).first);
        } 
        else if((*it).second== tools::_cid(double(0))) {
          fNtupleDColumnMap[index++] = fNtuple->find_column<double>((*it).first);
        }
        else {
          G4ExceptionDescription description;
          description << "      " 
                      << "Unsupported column type " << (*it).first;
          G4Exception("G4CsvAnalysisManager::OpenFile()",
                      "Analysis_W004", JustWarning, description);
        }
      }
    }
  }   

  fLockFileName = true;
  
#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("open", "analysis file", name);
#endif
  
  return true;
}  
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::Write() 
{
  // nothing to be done for Csv file
  G4bool result = true;
  
  // Write ASCII if activated
  if ( IsAscii() ) {
    result = WriteAscii();
  }   

  return result;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::CloseFile()
{
  G4bool result = true;

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "file", GetFullFileName());
#endif

  // reset data
  result = Reset();
  if ( ! result ) {
      G4ExceptionDescription description;
      description << "      " << "Resetting data failed";
      G4Exception("G4CsvAnalysisManager::CloseFile()",
                "Analysis_W002", JustWarning, description);
      result = false;       
  } 

  // close file
  fFile->close(); 
  fLockFileName = false;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", GetFullFileName());
#endif

  return true; 
} 
   
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH1(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               const G4String& /*unitName*/, 
                               const G4String& /*fcnName*/)
{
  ExceptionForHistograms("CreateH1");
  return 0;
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateH2(const G4String& /*name*/, 
                               const G4String& /*title*/, 
                               G4int /*nxbins*/, 
                               G4double /*xmin*/, G4double /*xmax*/,
                               G4int /*nybins*/, 
                               G4double /*ymin*/, G4double /*ymax*/,
                               const G4String& /*xunitName*/, 
                               const G4String& /*yunitName*/, 
                               const G4String& /*xfcnName*/,
                               const G4String& /*yfcnName*/)
{
  ExceptionForHistograms("CreateH2");
  return 0;
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1(G4int /*id*/,
                                G4int /*nbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/,
                                const G4String& /*unitName*/, 
                                const G4String& /*fcnName*/)
{                                
  ExceptionForHistograms("SetH1");
  return false;
}
  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2(G4int /*id*/,
                                G4int /*nxbins*/, 
                                G4double /*xmin*/, G4double /*xmax*/, 
                                G4int /*nybins*/, 
                                G4double /*ymin*/, G4double /*ymax*/,
                                const G4String& /*xunitName*/, 
                                const G4String& /*yunitName*/, 
                                const G4String& /*xfcnName*/,
                                const G4String& /*yfcnName*/)
{                                
  ExceptionForHistograms("SetH2");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::ScaleH1(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH1");
  return false;
}
                                  
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::ScaleH2(G4int /*id*/, G4double /*factor*/)
{
  ExceptionForHistograms("ScaleH2");
  return false;
}

//_____________________________________________________________________________
void G4CsvAnalysisManager::CreateNtuple(const G4String& name, 
                                        const G4String& title)
{
  if ( fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple already exists. "
                << "(Only one ntuple is currently supported.)";
    G4Exception("G4CsvAnalysisManager::CreateNtuple()",
                "Analysis_W006", JustWarning, description);
    return;       
  }

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple", name);
#endif

  // Create ntuple booking
  fNtupleBooking = new tools::ntuple_booking();
  fNtupleBooking->m_name = name;
  fNtupleBooking->m_title = title;
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open
  if ( fFile ) {
    fNtuple = new tools::wcsv::ntuple(*fFile);
           // ntuple object is deleted when closing a file
  }  

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple", name);
#endif
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple I column", name);
#endif

  if ( ! fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleIColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::wcsv::ntuple::column<int>* column 
      = fNtuple->create_column<int>(name);  
    fNtupleIColumnMap[index] = column;
  }  

  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple I column", name);
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple F column", name);
#endif

  if ( ! fNtupleBooking )  {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleFColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::wcsv::ntuple::column<float>* column 
      = fNtuple->create_column<float>(name);  
    fNtupleFColumnMap[index] = column;
  }  

  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple F column", name);
#endif

  return index + fFirstNtupleColumnId;       
}                                         

//_____________________________________________________________________________
G4int G4CsvAnalysisManager::CreateNtupleDColumn(const G4String& name)   
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple D column", name);
#endif

  if ( ! fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4CsvAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::wcsv::ntuple::column<double>* column 
      = fNtuple->create_column<double>(name);  
    fNtupleDColumnMap[index] = column;
  }
    
  fLockFirstNtupleColumnId = true;

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("create", "ntuple D column", name);
#endif

  return index + fFirstNtupleColumnId;       
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple I column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::AddNtupleRow()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL4 )
    fpVerboseL4->Message("add", "ntuple row", "");
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
  if ( fpVerboseL4 )
    fpVerboseL4->Message("add", "ntuple row", "");
#endif

  return true;
}

//_____________________________________________________________________________
tools::wcsv::ntuple* G4CsvAnalysisManager::GetNtuple() const
{
  return fNtuple;
}  


//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH1Nbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Nbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xmax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH1Width(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Xwidth");
  return 0;
}
  
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH2Nxbins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NXbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Xmin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Xmax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Xmin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2XWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XWidth");
  return 0;
}
  
//_____________________________________________________________________________
G4int G4CsvAnalysisManager::GetH2Nybins(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2NYbins");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Ymin(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymin");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2Ymax(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Ymax");
  return 0;
}
  
//_____________________________________________________________________________
G4double G4CsvAnalysisManager::GetH2YWidth(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YWidth");
  return 0;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1Title(G4int /*id*/, 
                                        const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1XAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH1YAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH1YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2Title(G4int /*id*/, 
                                        const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2Title");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2XAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2XAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2YAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2YAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4bool G4CsvAnalysisManager::SetH2ZAxisTitle(G4int /*id*/, 
                                             const G4String& /*title*/)
{
  ExceptionForHistograms("SetH2ZAxisTitle");
  return false;
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH1YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH1YAxisTitle");
  return "";
}


//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2Title(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2Title");
  return "";
}
  
//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2XAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2XAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2YAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2YAxisTitle");
  return "";
}

//_____________________________________________________________________________
G4String G4CsvAnalysisManager::GetH2ZAxisTitle(G4int /*id*/) const
{
  ExceptionForHistograms("GetH2ZAxisTitle");
  return "";
}
  
