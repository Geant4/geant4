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
/// \file hbook/src/ExG4HbookH2Manager.cc
/// \brief Implementation of the ExG4HbookH2Manager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookH2Manager.hh"
#include "ExG4HbookFileManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include <iostream>

using namespace G4Analysis;

const G4int ExG4HbookH2Manager::fgkDefaultH2HbookIdOffset = 100;

//_____________________________________________________________________________
ExG4HbookH2Manager::ExG4HbookH2Manager(const G4AnalysisManagerState& state)
 : G4VH2Manager(state),
   fBaseToolsManager("H2"),
   fFileManager(0),
   fH2HbookIdOffset(-1),
   fH2Vector(),
   fH2BookingVector(),
   fH2NameIdMap()
{
}

//_____________________________________________________________________________
ExG4HbookH2Manager::~ExG4HbookH2Manager()
{  
  // Delete h2
  Reset();

  // Delete h2 booking 
  std::vector<h2_booking*>::iterator it2;
  for ( it2 = fH2BookingVector.begin(); it2 != fH2BookingVector.end(); it2++ ) {
    delete *it2;
  }  
}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookH2Manager::SetH2HbookIdOffset()
{
// Set  fH2HbookIdOffset if needed

  if ( fH2HbookIdOffset == -1 ) {
    if ( fFirstId > 0 ) 
      fH2HbookIdOffset = 0;
    else
      fH2HbookIdOffset = 1;
        
    if ( fH2HbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "H2 will be defined in HBOOK with ID = G4_firstHistoId + 1";
      G4Exception("ExG4HbookH2Manager::SetH1HbookIdOffset",
                  "Analysis_W013", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookH2Manager::CreateH2sFromBooking()
{
// Create h2 from h2_booking.

  // Do nothing if any h2 histogram already exists
  // or no h2 histograms are booked
  if ( fH2Vector.size() || ( fH2BookingVector.size() == 0 ) ) return;       

  // Go to histograms directory
  if ( fFileManager->GetHistoDirectoryName() != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fFileManager->GetHistoDirectoryName().data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  
  
  // Create histograms
  G4int index = 0;
  std::vector<h2_booking*>::const_iterator it;
  for ( it = fH2BookingVector.begin(); it != fH2BookingVector.end(); ++it) {
    // Get information
    G4int id = index + fFirstId;    
    G4HnInformation* info = fHnManager->GetHnInformation(id, "CreateH2FromBooking");
    // Hbook index
    G4int hbookIndex = fH2HbookIdOffset + index + fFirstId;
    ++index;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("create from booking", "h2", info->GetName());
#endif

    // Create h2
    tools::hbook::h2* h2 
      = new tools::hbook::h2(hbookIndex, (*it)->fTitle, 
                             (*it)->fNxbins, (*it)->fXmin, (*it)->fXmax,
                             (*it)->fNybins, (*it)->fYmin, (*it)->fYmax);
    fH2Vector.push_back(h2);

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) { 
      G4ExceptionDescription description;
      description << " name : " << info->GetName() << " hbook index : " << hbookIndex; 
      fState.GetVerboseL3()->Message("create from booking", "h2", description);
    }  
#endif
  } 
  
  if ( fFileManager->GetHistoDirectoryName() != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
}   

//_____________________________________________________________________________
void ExG4HbookH2Manager::Reset()
{
// Reset histograms and ntuple  

  // Delete histograms
  std::vector<tools::hbook::h2*>::iterator it2;
  for (it2 = fH2Vector.begin(); it2 != fH2Vector.end(); it2++ ) {
    delete *it2;
  }  

  // Clear vectors
  fH2Vector.clear();
}  
 
//_____________________________________________________________________________
h2_booking*  ExG4HbookH2Manager::GetH2Booking(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH2BookingVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH2Booking()",
                  "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }
  
  return fH2BookingVector[index];
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.
// Not yet available for H2

  return false;
}  

//_____________________________________________________________________________
tools::hbook::h2*  ExG4HbookH2Manager::GetH2InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "ExG4HbookH2Manager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }

  // Do not return histogram if inactive 
  if ( fState.GetIsActivation()  && onlyIfActive && ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fH2Vector[index];
}
  
// 
// public methods
//

//_____________________________________________________________________________
G4int ExG4HbookH2Manager::CreateH2(const G4String& name, const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName,
                               const G4String& xbinSchemeName,
                               const G4String& ybinSchemeName)
                               
{
  // HBook does not support user defined binning for H2
  if ( xbinSchemeName != "linear" || ybinSchemeName != "linear" ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Logarithimc binning is not supported for H2.";
    G4Exception("ExG4HbookH2Manager::CreateH2",
                "Analysis_F015", FatalException, description);
  }              

#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H2", name);
#endif

  // Create h2 booking & information
  G4int index = fH2BookingVector.size();
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  G4String newTitle(title);
  UpdateTitle(newTitle, xunitName, xfcnName);  
  UpdateTitle(newTitle, yunitName, yfcnName);  

  h2_booking* h2Booking = new h2_booking(nxbins, xfcn(xmin), xfcn(xmax), 
                                         nybins, yfcn(ymin), yfcn(ymax)); 
           // h2_booking object is deleted in destructor
  h2Booking->fTitle = newTitle;
  fH2BookingVector.push_back(h2Booking);
  fHnManager
    ->AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                       xunit, yunit, xfcn, yfcn, xbinScheme, ybinScheme);
  
  // Set fH1HbookIdOffset if needed
  SetH2HbookIdOffset();
  
  // Hbook index
  G4int hbookIndex = fH2HbookIdOffset + index + fFirstId;

  // Create h2 if the file is open
  if ( fFileManager->IsFile() ) {
    // Go to histograms directory
    G4String histoPath = "//PAWC/LUN1/";
    if ( fFileManager->GetHistoDirectoryName() != "" ) {
      histoPath.append(fFileManager->GetHistoDirectoryName().data());
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
    if ( fFileManager->GetHistoDirectoryName() != "" ) {
      tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    }
  }    

  fLockFirstId = true;

#ifdef G4VERBOSE
    if ( fState.GetVerboseL2() ) {
      G4ExceptionDescription description;
      description << " name : " << name << " hbook index : " << hbookIndex; 
      fState.GetVerboseL2()->Message("create", "H2", description);
    }  
#endif

  fH2NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

//_____________________________________________________________________________
G4int ExG4HbookH2Manager::CreateH2(const G4String& /*name*/,  const G4String& /*title*/,
                          const std::vector<G4double>& /*xedges*/,
                          const std::vector<G4double>& /*yedges*/,
                          const G4String& /*xunitName*/, const G4String& /*yunitName*/,
                          const G4String& /*xfcnName*/, const G4String& /*yfcnName*/)
{                          
  // HBook does not support user defined binning for H2
  G4ExceptionDescription description;
  description 
    << "      " 
    << "User defined binning is not supported for H2.";
  G4Exception("ExG4HbookH2Manager::CreateH2",
              "Analysis_F015", FatalException, description);
  return 0;              
}              
                               
//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName,
                                const G4String& xbinScheme, const G4String& ybinScheme)
{                                
  // HBook does not support user defined binning for H2
  if ( xbinScheme != "linear" || ybinScheme != "linear" ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "Logarithimc binning is not supported for H2.";
    G4Exception("ExG4HbookH2Manager::CreateH2",
                "Analysis_F015", FatalException, description);
  }              

  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH2()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  G4HnInformation* hnInfo = fHnManager->GetHnInformation(id, "SetH2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H2", hnInfo->GetName());
#endif

  // Keep new parameters in booking
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  h2Booking->fNxbins = nxbins;
  h2Booking->fXmin = xfcn(xmin/xunit);
  h2Booking->fXmax = xfcn(xmax/xunit);
  h2Booking->fNybins = nybins;
  h2Booking->fYmin = yfcn(ymin/yunit);
  h2Booking->fYmax = yfcn(ymax/yunit);
  
  // Keep new parameters in information
  G4HnDimensionInformation* xInfo 
    = hnInfo->GetHnDimensionInformation(G4HnInformation::kX);
  xInfo->fUnitName = xunitName;
  xInfo->fFcnName = xfcnName;
  xInfo->fUnit = xunit;
  xInfo->fFcn = xfcn;
    
  G4HnDimensionInformation* yInfo 
    = hnInfo->GetHnDimensionInformation(G4HnInformation::kY);
  yInfo->fUnitName = yunitName;
  yInfo->fFcnName = yfcnName;
  yInfo->fUnit = yunit;
  yInfo->fFcn = yfcn;
  fHnManager->SetActivation(id, true); 

  G4String newTitle(h2Booking->fTitle);
  UpdateTitle(newTitle, xunitName, xfcnName);  
  UpdateTitle(newTitle, yunitName, yfcnName);  
  h2Booking->fTitle = newTitle;  
  
  // Re-configure histogram if it was already defined
  if ( fH2Vector.size() ) {
    tools::hbook::h2* h2 = GetH2(id);
    h2->configure(nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                  nybins, yfcn(ymin/yunit), yfcn(ymax/yunit));
  }  
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2(G4int /*id*/,
                            const std::vector<G4double>& /*xedges*/,
                            const std::vector<G4double>& /*yedges*/,
                            const G4String& /*xunitName*/, const G4String& /*yunitName*/,
                            const G4String& /*xfcnName*/, const G4String& /*yfcnName*/)
{                          
  // HBook does not support user defined binning for H2
  G4ExceptionDescription description;
  description 
    << "      " 
    << "User defined binning is not supported for H2.";
  G4Exception("ExG4HbookH2Manager::CreateH2",
              "Analysis_F015", FatalException, description);

  return false;              
}              
                            
//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::ScaleH2(G4int id, G4double factor)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "ScaleH2", false, false);
  if ( ! h2 ) return false;

  return h2->scale(factor);
}  
                           
//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::FillH2(G4int id, 
                                       G4double xvalue, G4double yvalue,
                                       G4double weight)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "FillH2", true, false);
  if ( ! h2 ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) return false; 

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillH2");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillH2");
  h2->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
           yInfo->fFcn(yvalue/yInfo->fUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue 
                << " xfcn(xvalue/xunit) " <<  xInfo->fFcn(xvalue/xInfo->fUnit) 
                << " yvalue " << yvalue
                << " yfcn(yvalue/yunit) " <<  yInfo->fFcn(yvalue/yInfo->fUnit) 
                << " weight " << weight;
    fState.GetVerboseL4()->Message("fill", "H2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
tools::hbook::h2*  ExG4HbookH2Manager::GetH2(G4int id, G4bool warn,
                                                   G4bool onlyIfActive) const 
{
  return GetH2InFunction(id, "GetH2", warn, onlyIfActive);
}

//_____________________________________________________________________________
G4int  ExG4HbookH2Manager::GetH2Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH2NameIdMap.find(name);
  if ( it ==  fH2NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "ExG4HbookH2Manager::GetH2Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int ExG4HbookH2Manager::GetH2Nxbins(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2NXbins");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetNbins(h2->axis_x());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Xmin");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetMin(h2->axis_x());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2Xmax(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Xmax");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetMin(h2->axis_x());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2XWidth(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2XWidth", true, false);
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetWidth(h2->axis_x());
}  

//_____________________________________________________________________________
G4int ExG4HbookH2Manager::GetH2Nybins(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2NYbins");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetNbins(h2->axis_y());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Ymin");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetMin(h2->axis_y());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2Ymax(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2Ymax");
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetMax(h2->axis_y());
}  

//_____________________________________________________________________________
G4double ExG4HbookH2Manager::GetH2YWidth(G4int id) const
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2YWidth", true, false);
  if ( ! h2 ) return 0;
  
  return fBaseToolsManager.GetWidth(h2->axis_y());
}  

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2Title(G4int id, const G4String& title)
{
  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH2Title()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  h2Booking->fTitle = title;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2XAxisTitle");
  if ( ! h2 ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2, ExG4HbookBaseHnManager::kX, title);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2YAxisTitle");
  if ( ! h2 ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2, ExG4HbookBaseHnManager::kY, title);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "SetH2ZAxisTitle");
  if ( ! h2 ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2, ExG4HbookBaseHnManager::kZ, title);
}  

//_____________________________________________________________________________
G4String ExG4HbookH2Manager::GetH2Title(G4int id) const
{
  h2_booking* h2Booking = GetH2Booking(id, false);
  if ( ! h2Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetH2Title()",
                "Analysis_W011", JustWarning, description);
    return "";
  }
  
  return h2Booking->fTitle;
}  


//_____________________________________________________________________________
G4String ExG4HbookH2Manager::GetH2XAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2XAxisTitle");
  if ( ! h2 ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2, ExG4HbookBaseHnManager::kX);
} 

//_____________________________________________________________________________
G4String ExG4HbookH2Manager::GetH2YAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2YAxisTitle");
  if ( ! h2 ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2, ExG4HbookBaseHnManager::kY);
}  

//_____________________________________________________________________________
G4String ExG4HbookH2Manager::GetH2ZAxisTitle(G4int id) const 
{
  tools::hbook::h2* h2 = GetH2InFunction(id, "GetH2ZAxisTitle");
  if ( ! h2 ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2, ExG4HbookBaseHnManager::kZ);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH2Manager::SetH2HbookIdOffset(G4int offset) 
{
  if ( fH2Vector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set H2HbookIdOffset as some H2 histogramms already exist.";
    G4Exception("G4HbookAnalysisManager::SetH2HbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }

  if ( fFirstId + offset < 1 ) {
    G4ExceptionDescription description;
    description << "The first histogram HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetH1HbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }
  
  fH2HbookIdOffset = offset;
  return true;
}  

#endif
