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
//
/// \file common/analysis/src/ExG4HbookAnalysisManager.cc
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
   fH2Vector(),
   fH1BookingVector(),
   fH2BookingVector(),
   fH1NameIdMap(),  
   fH2NameIdMap(),  
   fNtuple(0),
   fNtupleBooking(0),
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
  // Delete h1, h2, ntuple
  Reset();

  // Delete h1, h2 booking 
  std::vector<h1_booking*>::iterator it;
  for ( it = fH1BookingVector.begin(); it != fH1BookingVector.end(); it++ ) {
    delete *it;
  }  
  std::vector<h2_booking*>::iterator it2;
  for ( it2 = fH2BookingVector.begin(); it2 != fH2BookingVector.end(); it2++ ) {
    delete *it2;
  }  

  delete fNtupleBooking;
  delete fFile;  

  fgInstance = 0;
}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::SetH1HbookIdOffset()
{
// Set  fH1HbookIdOffset if needed

  if ( fH1HbookIdOffset == -1 ) {
    if ( fFirstHistoId > 0 ) 
      fH1HbookIdOffset = 0;
    else
      fH1HbookIdOffset = 1;
        
    if ( fH1HbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "H1 will be defined in HBOOK with ID = G4_firstHistoId + 1";
      G4Exception("ExG4HbookAnalysisManager::SetH1HbookIdOffset()",
                  "Analysis_W011", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::SetH2HbookIdOffset()
{
// Set  fH2HbookIdOffset if needed

  if ( fH2HbookIdOffset == -1 ) {
    if ( fFirstHistoId > 0 ) 
      fH2HbookIdOffset = 0;
    else
      fH2HbookIdOffset = 1;
        
    if ( fH2HbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "H2 will be defined in HBOOK with ID = G4_firstHistoId + 1";
      G4Exception("ExG4HbookAnalysisManager::SetH1HbookIdOffset",
                  "Analysis_W011", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::CreateH1FromBooking()
{
// Create h1 from h1_booking.

  // Do nothing if any h1 histogram already exists
  // or no h1 histograms are booked
  if ( fH1Vector.size() || ( fH1BookingVector.size() == 0 ) ) return;       

  // Go to histograms directory if defined
  if ( fHistoDirectoryName != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fHistoDirectoryName.data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  

  // Create histograms
  G4int index = 0;
  std::vector<h1_booking*>::const_iterator it;
  for ( it = fH1BookingVector.begin(); it != fH1BookingVector.end(); ++it) {
    // Get information
    G4int id = index + fFirstHistoId;    
    G4HnInformation* info = GetH1Information(id);
    // Hbook index
    G4int hbookIndex = fH1HbookIdOffset + index + fFirstHistoId;
    ++index;

#ifdef G4VERBOSE
    if ( fpVerboseL4 ) 
      fpVerboseL4->Message("create from booking", "h1", info->fName);
#endif

    // Create h1
    tools::hbook::h1* h1 
      = new tools::hbook::h1(hbookIndex, (*it)->fTitle, 
                             (*it)->fNbins, (*it)->fXmin, (*it)->fXmax);
    fH1Vector.push_back(h1);

#ifdef G4VERBOSE
    if ( fpVerboseL3 ) { 
      G4ExceptionDescription description;
      description << " name : " << info->fName << " hbook index : " << hbookIndex; 
      fpVerboseL3->Message("create from booking", "h1", description);
    }  
#endif
  } 
  
  if ( fHistoDirectoryName != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
}   

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::CreateH2FromBooking()
{
// Create h2 from h2_booking.

  // Do nothing if any h2 histogram already exists
  // or no h2 histograms are booked
  if ( fH2Vector.size() || ( fH2BookingVector.size() == 0 ) ) return;       

  // Go to histograms directory
  if ( fHistoDirectoryName != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fHistoDirectoryName.data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  
  
  // Create histograms
  G4int index = 0;
  std::vector<h2_booking*>::const_iterator it;
  for ( it = fH2BookingVector.begin(); it != fH2BookingVector.end(); ++it) {
    // Get information
    G4int id = index + fFirstHistoId;    
    G4HnInformation* info = GetH2Information(id);
    // Hbook index
    G4int hbookIndex = fH2HbookIdOffset + index + fFirstHistoId;
    ++index;

#ifdef G4VERBOSE
    if ( fpVerboseL3 ) 
      fpVerboseL3->Message("create from booking", "h2", info->fName);
#endif

    // Create h2
    tools::hbook::h2* h2 
      = new tools::hbook::h2(hbookIndex, (*it)->fTitle, 
                             (*it)->fNxbins, (*it)->fXmin, (*it)->fXmax,
                             (*it)->fNybins, (*it)->fYmin, (*it)->fYmax);
    fH2Vector.push_back(h2);

#ifdef G4VERBOSE
    if ( fpVerboseL3 ) { 
      G4ExceptionDescription description;
      description << " name : " << info->fName << " hbook index : " << hbookIndex; 
      fpVerboseL3->Message("create from booking", "h2", description);
    }  
#endif
  } 
  
  if ( fHistoDirectoryName != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
}   

//_____________________________________________________________________________
void ExG4HbookAnalysisManager::CreateNtupleFromBooking()
{
// Create ntuple from ntuple_booking.

  if ( fNtuple || (! fNtupleBooking) ) return;       

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create from booking", "ntuple", fNtupleBooking->m_name);
#endif

  // Create an "ntuple" directory both in memory and in the file
  fFile->cd_home();      //go under //PAWC/LUN1
  if ( fNtupleDirectoryName == "" )
    fFile->mkcd(fgkDefaultNtupleDirectoryName.data());
  else  
    fFile->mkcd(fNtupleDirectoryName.data());
  fLockNtupleDirectoryName = true;

  // Define ntuple ID in HBOOK
  if ( fNtupleHbookId == -1 ) fNtupleHbookId = fgkDefaultNtupleHbookId;
  
  // We should be under //PAWC/LUN1/ntuple
  fNtuple = new tools::hbook::wntuple(fNtupleHbookId, G4cout, *fNtupleBooking);
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
        G4Exception("G4HbookAnalysisManager::CreateNtupleFromBooking()",
                    "Analysis_W004", JustWarning, description);
      }
    }
  }
  FinishNtuple();

#ifdef G4VERBOSE
  if ( fpVerboseL3 ) 
    fpVerboseL3->Message("create from booking", "ntuple", fNtupleBooking->m_name);
#endif
}   

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
void ExG4HbookAnalysisManager::Reset()
{
// Reset histograms and ntuple  

  // Delete histograms
  std::vector<tools::hbook::h1*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  }  
  
  std::vector<tools::hbook::h2*>::iterator it2;
  for (it2 = fH2Vector.begin(); it2 != fH2Vector.end(); it2++ ) {
    delete *it2;
  }  

  // Clear vectors
  fH1Vector.clear();
  fH2Vector.clear();

  // Delete ntuple
  delete fNtuple;
  fNtuple = 0;
}  
 
//_____________________________________________________________________________
void ExG4HbookAnalysisManager::UpdateTitle(G4String& title, 
                                           const G4String& unitName, 
                                           const G4String& fcnName) const
{
// The units are not within [ ] as this causes strange font effect in 
// browsing the hbook file in PAW

  if ( fcnName != "none" )  { title += " "; title += fcnName; title += "("; }
  if ( unitName != "none" ) { title += " "; title += unitName; title += " ";}
  if ( fcnName != "none" )  { title += ")"; }
}  
                                                          
//_____________________________________________________________________________
h1_booking*  ExG4HbookAnalysisManager::GetH1Booking(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1BookingVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH1Booking()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }

  return fH1BookingVector[index];
}

//_____________________________________________________________________________
h2_booking*  ExG4HbookAnalysisManager::GetH2Booking(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH2BookingVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH2Booking()",
                  "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  return fH2BookingVector[index];
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// (Only H1 implemented by now)
// According to the implementation by Michel Maire, originally in
// extended examples.

  // h1 histograms
  for ( G4int i=0; i<G4int(fH1Vector.size()); ++i ) {
    G4int id = i + fFirstHistoId;
    G4HnInformation* info = GetH1Information(id); 
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->fAscii ) continue; 
    tools::hbook::h1* h1 = fH1Vector[i];

#ifdef G4VERBOSE
    if ( fpVerboseL3 ) 
      fpVerboseL3->Message("write on ascii", "h1", info->fName);
#endif
  
    output << "\n  1D histogram " << id << ": " << h1->title() 
           << "\n \n \t     X \t\t     Y" << G4endl;
    
    for (G4int j=0; j< G4int(h1->axis().bins()); ++j) {
       output << "  " << j << "\t" 
              << h1->axis().bin_center(j) << "\t"
              << h1->bin_height(j) << G4endl;
    } 
  }
  
  return true;
}  

//_____________________________________________________________________________
tools::hbook::h1*  ExG4HbookAnalysisManager::GetH1InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "ExG4HbookAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  // Do not return histogram if inactive 
  if ( fActivation && onlyIfActive && ( ! GetActivation(kH1, id) ) ) {
    return 0; 
  }  
  
  return fH1Vector[index];
}  
                                      
//_____________________________________________________________________________
tools::hbook::h2*  ExG4HbookAnalysisManager::GetH2InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstHistoId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "ExG4HbookAnalysisManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }

  // Do not return histogram if inactive 
  if ( fActivation  && onlyIfActive && ( ! GetActivation(kH2, id) ) ) {
    return 0; 
  }  
  
  return fH2Vector[index];
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
  // Keep file name
  fFileName =  fileName;

  // Add file extension .root if no extension is given
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
  // the five upper lines could have been done with :
  //fFile->cd_home();
  //fFile->mkcd("histo");

  // Create h1 histrograms if any is booked
  if ( ( fH1Vector.size() == 0 ) && ( fH1BookingVector.size() ) )  
    CreateH1FromBooking();

  // Create h2 histrograms if any is booked
  if ( ( fH2Vector.size() == 0 ) && ( fH2BookingVector.size() ) )  
    CreateH2FromBooking();

  // Create ntuple if it is booked
  if ( fNtupleBooking && ( ! fNtuple ) ) 
    CreateNtupleFromBooking();

  fLockFileName = true;
  fLockHistoDirectoryName = true;

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
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("write", "file", GetFullFileName());
#endif

  // ntuple 
  //if ( fNtuple ) fNtuple->add_row_end();

  // Return to //PAWC/LUN1 :
  //tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  G4bool result = fFile->write();  

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("write", "file", GetFullFileName(), result);
#endif

  // Write ASCII if activated
  if ( IsAscii() ) {
    G4bool result2 = WriteAscii();
    result = result && result2;
  }   

  return result;  
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::CloseFile()
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("close", "file", GetFullFileName());
#endif

  // reset data
  Reset();

  // close file
  G4bool result = fFile->close();  
  fLockFileName = false;

#ifdef G4VERBOSE
  if ( fpVerboseL1 ) 
    fpVerboseL1->Message("close", "file", GetFullFileName(), result);
#endif

  return result;
} 
   
//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateH1(const G4String& name, const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "H1", name);
#endif

  // Create h1 booking & information
  G4int index = fH1BookingVector.size();
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  G4String newTitle(title);
  UpdateTitle(newTitle, unitName, fcnName);  
  h1_booking* h1Booking = new h1_booking(nbins, fcn(xmin), fcn(xmax)); 
           // h1_booking object is deleted in destructor
  h1Booking->fTitle = newTitle;
  fH1BookingVector.push_back(h1Booking);
  AddH1Information(name, unitName, fcnName, unit, fcn);
  
  // Set  fH1HbookIdOffset if needed
  SetH1HbookIdOffset();
  
  // Hbook index
  G4int hbookIndex = fH1HbookIdOffset + index + fFirstHistoId;

  // Create h1 if the file is open
  if ( fFile) {
    // Go to histograms directory
    G4String histoPath = "//PAWC/LUN1/";
    if ( fHistoDirectoryName != "" ) {
      histoPath.append(fHistoDirectoryName.data());
      tools::hbook::CHCDIR(histoPath.data()," ");
    }  
    tools::hbook::CHCDIR(histoPath.data()," ");
    
    // Create histogram    
    tools::hbook::h1* h1 
      = new tools::hbook::h1(hbookIndex, newTitle, nbins, fcn(xmin), fcn(xmax));
            // h1 objects are deleted when closing a file.
    fH1Vector.push_back(h1);
 
    if ( fHistoDirectoryName != "" ) {
      // Return to //PAWC/LUN1 :
      tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    }
  }
  
  fLockFirstHistoId = true;

#ifdef G4VERBOSE
    if ( fpVerboseL2 ) { 
      G4ExceptionDescription description;
      description << " name : " << name << " hbook index : " << hbookIndex; 
      fpVerboseL2->Message("create", "H1", description);
    }  
#endif

  fH1NameIdMap[name] = index + fFirstHistoId;
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateH2(const G4String& name, const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "H2", name);
#endif

  // Create h2 booking & information
  G4int index = fH2BookingVector.size();
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4String newTitle(title);
  UpdateTitle(newTitle, xunitName, xfcnName);  
  UpdateTitle(newTitle, yunitName, yfcnName);  

  h2_booking* h2Booking = new h2_booking(nxbins, xfcn(xmin), xfcn(xmax), 
                                         nybins, yfcn(ymin), yfcn(ymax)); 
           // h2_booking object is deleted in destructor
  h2Booking->fTitle = newTitle;
  fH2BookingVector.push_back(h2Booking);
  AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                   xunit, yunit, xfcn, yfcn);
  
  // Set fH1HbookIdOffset if needed
  SetH2HbookIdOffset();
  
  // Hbook index
  G4int hbookIndex = fH2HbookIdOffset + index + fFirstHistoId;

  // Create h2 if the file is open
  if ( fFile) {
    // Go to histograms directory
    G4String histoPath = "//PAWC/LUN1/";
    if ( fHistoDirectoryName != "" ) {
      histoPath.append(fHistoDirectoryName.data());
    }  
    tools::hbook::CHCDIR(histoPath.data()," ");

    // Create histogram    
    tools::hbook::h2* h2 
      = new tools::hbook::h2(hbookIndex, title, 
                             nxbins, xfcn(xmin), xfcn(xmax), 
                             nybins, yfcn(ymin), yfcn(ymax));
            // h2 objects are deleted when closing a file.
    fH2Vector.push_back(h2);

    // Return to //PAWC/LUN1 
    if ( fHistoDirectoryName != "" ) {
      tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    }
  }    

  fLockFirstHistoId = true;

#ifdef G4VERBOSE
    if ( fpVerboseL2 ) {
      G4ExceptionDescription description;
      description << " name : " << name << " hbook index : " << hbookIndex; 
      fpVerboseL2->Message("create", "H2", description);
    }  
#endif

  fH2NameIdMap[name] = index + fFirstHistoId;
  return index + fFirstHistoId;
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH1(G4int id,
                                   G4int nbins, G4double xmin, G4double xmax,
                                   const G4String& unitName, 
                                   const G4String& fcnName)
{                                
  h1_booking* h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH1()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  G4HnInformation* info = GetH1Information(id);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("configure", "H1", info->fName);
#endif

  // Keep new parameters in booking & information
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  h1Booking->fNbins = nbins;
  h1Booking->fXmin = fcn(xmin);
  h1Booking->fXmax = fcn(xmax);
  info->fXUnitName = unitName;
  info->fYUnitName = unitName;
  info->fXFcnName = fcnName;
  info->fYFcnName = fcnName;
  info->fXUnit = unit;
  info->fYUnit = unit;
  info->fXFcn = fcn;
  info->fYFcn = fcn;
  SetActivation(kH1, id, true); 

  G4String newTitle(h1Booking->fTitle);
  UpdateTitle(newTitle, unitName, fcnName);  
  h1Booking->fTitle = newTitle;  
  
  // Re-configure histogram if it was already defined
  if ( fH1Vector.size() ) {
    tools::hbook::h1* h1 = GetH1(id);
    h1->configure(nbins, fcn(xmin), fcn(xmax));
    G4cout << " h1 = " << h1 << G4endl;   
  }  
  
  return true;
}
  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{                                
  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH2()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  G4HnInformation* info = GetH2Information(id);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("configure", "H2", info->fName);
#endif

  // Keep new parameters in booking & information
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  h2Booking->fNxbins = nxbins;
  h2Booking->fXmin = xfcn(xmin);
  h2Booking->fXmax = xfcn(xmax);
  h2Booking->fNybins = nybins;
  h2Booking->fYmin = yfcn(ymin);
  h2Booking->fYmax = yfcn(ymax);
  
  info->fXUnitName = xunitName;
  info->fYUnitName = yunitName;
  info->fXFcnName = xfcnName;
  info->fYFcnName = yfcnName;
  info->fXUnit = xunit;
  info->fYUnit = yunit;
  info->fXFcn = xfcn;
  info->fYFcn = yfcn;
  SetActivation(kH2, id, true); 

  G4String newTitle(h2Booking->fTitle);
  UpdateTitle(newTitle, xunitName, xfcnName);  
  UpdateTitle(newTitle, yunitName, yfcnName);  
  h2Booking->fTitle = newTitle;  
  
  // Re-configure histogram if it was already defined
  if ( fH2Vector.size() ) {
    tools::hbook::h2* h2 = GetH2(id);
    h2->configure(nxbins, xfcn(xmin), xfcn(xmax), 
                  nybins, yfcn(ymin), yfcn(ymax));
  }  
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::ScaleH1(G4int id, G4double factor)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "ScaleH1", false, false);
  if ( ! h1 ) return false;

  return h1->scale(factor);
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::ScaleH2(G4int id, G4double factor)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "ScaleH2", false, false);
  if ( ! h2 ) return false;

  return h2->scale(factor);
}  
                           
//_____________________________________________________________________________
void ExG4HbookAnalysisManager::CreateNtuple(const G4String& name, 
                                            const G4String& title)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple", name);
#endif

  if ( fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple already exists. "
                << "(Only one ntuple is currently supported.)";
    G4Exception("G4HbookAnalysisManager::CreateNtuple()",
                "Analysis_W006", JustWarning, description);
    return;       
  }

  // Create an "ntuple" directory both in memory and in the file
  if ( fFile ) {
    fFile->cd_home();      //go under //PAWC/LUN1
    if ( fNtupleDirectoryName == "" )
      fFile->mkcd(fgkDefaultNtupleDirectoryName.data());
    else  
      fFile->mkcd(fNtupleDirectoryName.data());
    fLockNtupleDirectoryName = true;
  }    

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple", name);
#endif

  // Define ntuple ID in HBOOK
  if ( fNtupleHbookId == -1 ) fNtupleHbookId = fgkDefaultNtupleHbookId;
  
  // Create ntuple booking
  fNtupleBooking = new tools::ntuple_booking();
  fNtupleBooking->m_name = name;
  fNtupleBooking->m_title = title;
           // ntuple booking object is deleted in destructor

  // Create ntuple if the file is open
  // We should be under //PAWC/LUN1/ntuple
  if ( fFile ) {
    fNtuple = new tools::hbook::wntuple(fNtupleHbookId, name);
           // ntuple object is deleted when closing a file
  }  

#ifdef G4VERBOSE
  if ( fpVerboseL2 ) {
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << fNtupleHbookId; 
    fpVerboseL2->Message("create", "ntuple", description);
  }  
#endif
}                                         

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::CreateNtupleIColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple I column", name);
#endif

  if ( ! fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleIColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<int>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::hbook::wntuple::column<int>* column 
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
G4int ExG4HbookAnalysisManager::CreateNtupleFColumn(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple F column", name);
#endif

  if ( ! fNtupleBooking )  {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleFColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<float>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::hbook::wntuple::column<float>* column 
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
G4int ExG4HbookAnalysisManager::CreateNtupleDColumn(const G4String& name) 
{
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("create", "ntuple D column", name);
#endif

  if ( ! fNtupleBooking ) {
    G4ExceptionDescription description;
    description << "      " 
                << "Ntuple has to be created first. ";
    G4Exception("G4HbookAnalysisManager::CreateNtupleDColumn()",
                "Analysis_W005", JustWarning, description);
    return -1;       
  }

  // Save column info in booking
  G4int index = fNtupleBooking->m_columns.size();
  fNtupleBooking->add_column<double>(name);  
 
  // Create column if ntuple already exists
  if ( fNtuple ) {
    tools::hbook::wntuple::column<double>* column 
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
void ExG4HbookAnalysisManager::FinishNtuple()
{ 
  if ( ! fNtuple ) return;

#ifdef G4VERBOSE
  if ( fpVerboseL4 ) 
    fpVerboseL4->Message("finish", "ntuple", fNtupleBooking->m_name);
#endif

  // Return to //PAWC/LUN1 :
  tools::hbook::CHCDIR("//PAWC/LUN1"," ");

  //fNtuple->add_row_beg();
#ifdef G4VERBOSE
  if ( fpVerboseL2 ) 
    fpVerboseL2->Message("finish", "ntuple", fNtupleBooking->m_name);
#endif
}
  
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillH1(G4int id, G4double value, G4double weight)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "FillH1", true, false);
  if ( ! h1 ) return false;

  if ( fActivation && ( ! GetActivation(kH1, id) ) ) {
    //G4cout << "Skipping FillH1 for " << id << G4endl; 
    return false; 
  }  

  G4HnInformation* info = GetInformation(kH1, id);
  h1->fill(info->fXFcn(value/info->fXUnit), weight);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillH2(G4int id, 
                                       G4double xvalue, G4double yvalue,
                                       G4double weight)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "FillH2", true, false);
  if ( ! h2 ) return false;

  if ( fActivation && ( ! GetActivation(kH2, id) ) ) return false; 

  G4HnInformation* info = GetInformation(kH2, id);
  h2->fill(info->fXFcn(xvalue/info->fXUnit), 
           info->fYFcn(yvalue/info->fYUnit), weight);
#ifdef G4VERBOSE
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
    fpVerboseL4->Message("fill", "H2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleIColumn(G4int id, G4int value)
{
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple I column", description);
  }  
#endif
 return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleFColumn(G4int id, G4float value)
{
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple F column", description);
  }  
#endif
  return true;       
}                                         
//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::FillNtupleDColumn(G4int id, G4double value)
{
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
  if ( fpVerboseL4 ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fpVerboseL4->Message("fill", "ntuple D column", description);
  }  
#endif
  return true;       
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::AddNtupleRow()
{ 
#ifdef G4VERBOSE
  if ( fpVerboseL4 )
    fpVerboseL4->Message("add", "ntuple row", "");
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
  if ( fpVerboseL3 )
    fpVerboseL3->Message("add", "ntuple row", "");
#endif
  return true;
}
 
//_____________________________________________________________________________
tools::hbook::h1*  ExG4HbookAnalysisManager::GetH1(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const 
{
  return GetH1InFunction(id, "GetH1", warn, onlyIfActive);
}

//_____________________________________________________________________________
tools::hbook::h2*  ExG4HbookAnalysisManager::GetH2(G4int id, G4bool warn,
                                                   G4bool onlyIfActive) const 
{
  return GetH2InFunction(id, "GetH2", warn, onlyIfActive);
}

//_____________________________________________________________________________
G4int  ExG4HbookAnalysisManager::GetH1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH1NameIdMap.find(name);
  if ( it ==  fH1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "ExG4HbookAnalysisManager::GetH1Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int  ExG4HbookAnalysisManager::GetH2Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH2NameIdMap.find(name);
  if ( it ==  fH2NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "ExG4HbookAnalysisManager::GetH2Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
tools::hbook::wntuple* ExG4HbookAnalysisManager::GetNtuple() const
{
  return fNtuple;
}

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::GetH1Nbins(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Nbins");
  if ( ! h1 ) return 0;
  
  return h1->axis().bins();
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH1Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Xmin");
  if ( ! h1 ) return 0;
  
  G4HnInformation* info = GetInformation(kH1, id);
  return info->fXFcn(h1->axis().lower_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH1Xmax(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Xmax");
  if ( ! h1 ) return 0;
  
  G4HnInformation* info = GetInformation(kH1, id);
  return info->fXFcn(h1->axis().upper_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH1Width(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1XWidth", true, false);
  if ( ! h1 ) return 0;
  
  G4int nbins = h1->axis().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH1Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH1, id);
  return ( info->fXFcn(h1->axis().upper_edge()) 
           - info->fXFcn(h1->axis().lower_edge()))*info->fXUnit/nbins;
}  

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::GetH2Nxbins(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2NXbins");
  if ( ! h2 ) return 0;
  
  return h2->axis_x().bins();
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Xmin");
  if ( ! h2 ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fXFcn(h2->axis_x().lower_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2Xmax(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Xmax");
  if ( ! h2 ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fXFcn(h2->axis_x().upper_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2XWidth(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2XWidth", true, false);
  if ( ! h2 ) return 0;
  
  G4int nbins = h2->axis_x().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH2, id);
  return ( info->fXFcn(h2->axis_x().upper_edge()) 
           - info->fXFcn(h2->axis_x().lower_edge()))*info->fXUnit/nbins;
}  

//_____________________________________________________________________________
G4int ExG4HbookAnalysisManager::GetH2Nybins(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2NYbins");
  if ( ! h2 ) return 0;
  
  return h2->axis_y().bins();
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Ymin");
  if ( ! h2 ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fYFcn(h2->axis_y().lower_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2Ymax(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Ymax");
  if ( ! h2 ) return 0;
  
  G4HnInformation* info = GetInformation(kH2, id);
  return info->fYFcn(h2->axis_y().upper_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double ExG4HbookAnalysisManager::GetH2YWidth(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2YWidth", true, false);
  if ( ! h2 ) return 0;
  
  G4int nbins = h2->axis_y().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = GetInformation(kH2, id);
  return ( info->fYFcn(h2->axis_y().upper_edge()) 
           - info->fYFcn(h2->axis_y().lower_edge()))*info->fYUnit/nbins;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH1Title(G4int id, const G4String& title)
{
  h1_booking* h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH1Title()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  h1Booking->fTitle = title;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH1XAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "SetH1XAxisTitle");
  if ( ! h1 ) return false;
  
  h1->add_annotation(tools::hbook::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH1YAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "SetH1YAxisTitle");
  if ( ! h1 ) return false;
  
  h1->add_annotation(tools::hbook::key_axis_y_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2Title(G4int id, const G4String& title)
{
  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH2Title()",
                "Analysis_W007", JustWarning, description);
    return false;
  }

  h2Booking->fTitle = title;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2XAxisTitle");
  if ( ! h2 ) return false;
  
  h2->add_annotation(tools::hbook::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2YAxisTitle");
  if ( ! h2 ) return false;
  
  h2->add_annotation(tools::hbook::key_axis_x_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4bool ExG4HbookAnalysisManager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2ZAxisTitle");
  if ( ! h2 ) return false;
  
  h2->add_annotation(tools::hbook::key_axis_z_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH1Title(G4int id) const
{
  h1_booking* h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetH1Title()",
                "Analysis_W007", JustWarning, description);
    return "";
  }
  
  return h1Booking->fTitle;
}  


//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH1XAxisTitle(G4int id) const 
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1XAxisTitle");
  if ( ! h1 ) return "";
  
  G4String title;
  G4bool result = h1->annotation(tools::hbook::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h1 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH1XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH1YAxisTitle(G4int id) const 
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1YAxisTitle");
  if ( ! h1 ) return "";
  
  G4String title;
  G4bool result = h1->annotation(tools::hbook::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h1 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH1YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH2Title(G4int id) const
{
  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetH2Title()",
                "Analysis_W007", JustWarning, description);
    return "";
  }
  
  return h2Booking->fTitle;
}  


//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH2XAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2XAxisTitle");
  if ( ! h2 ) return "";
  
  G4String title;
  G4bool result = h2->annotation(tools::hbook::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h2 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH2XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
} 

//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH2YAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2YAxisTitle");
  if ( ! h2 ) return "";
  
  G4String title;
  G4bool result = h2->annotation(tools::hbook::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h2 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH2YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String ExG4HbookAnalysisManager::GetH2ZAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2ZAxisTitle");
  if ( ! h2 ) return "";
  
  G4String title;
  G4bool result = h2->annotation(tools::hbook::key_axis_z_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get z_axis title for h2 id = " << id << ").";
    G4Exception("ExG4HbookAnalysisManager::GetH2ZAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
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
