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
/// \file hbook/src/ExG4HbookH1Manager.cc
/// \brief Implementation of the ExG4HbookH1Manager class

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookH1Manager.hh"
#include "ExG4HbookFileManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
ExG4HbookH1Manager::ExG4HbookH1Manager(const G4AnalysisManagerState& state)
 : G4VH1Manager(state),
   fBaseToolsManager("H1"),
   fFileManager(0),
   fH1HbookIdOffset(-1),
   fH1Vector(),
   fH1BookingVector(),
   fH1NameIdMap()
{
}

//_____________________________________________________________________________
ExG4HbookH1Manager::~ExG4HbookH1Manager()
{  
  // Delete h1
  Reset();

  // Delete h1 booking 
  std::vector<h1_booking*>::iterator it;
  for ( it = fH1BookingVector.begin(); it != fH1BookingVector.end(); it++ ) {
    delete *it;
  }  
}

//
// utility functions
//

namespace {

//_____________________________________________________________________________
void ConvertToFloat(const std::vector<G4double>& doubleVector,
                    std::vector<float>& floatVector)
{
  for (G4int i=0; i<G4int(doubleVector.size()); ++i) 
    floatVector.push_back((float)doubleVector[i]);
}                        

//_____________________________________________________________________________
void UpdateH1Information(G4HnInformation* hnInformation,
                          const G4String& unitName, 
                          const G4String& fcnName,
                          G4BinScheme binScheme)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  G4HnDimensionInformation* information 
    = hnInformation->GetHnDimensionInformation(G4HnInformation::kX);
  information->fUnitName = unitName;
  information->fFcnName = fcnName;
  information->fUnit = unit;
  information->fFcn = fcn;
  information->fBinScheme = binScheme;
}  

//_____________________________________________________________________________
h1_booking* CreateH1Booking(const G4String& title,
                   G4int nbins, G4double xmin, G4double xmax,
                   const G4String& unitName,
                   const G4String& fcnName,
                   G4BinScheme binScheme)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);

  h1_booking* h1Booking = 0; 
  if ( binScheme != kLogBinScheme ) {
    if ( binScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookH1Manager::CreateH1",
                "Analysis_W013", JustWarning, description);
    }              
    h1Booking = new h1_booking(nbins, fcn(xmin/unit), fcn(xmax/unit)); 
                    // h1_booking object is deleted in destructor
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, unit, fcn, binScheme, edges);
    h1Booking = new h1_booking(edges); 
                    // h1_booking object is deleted in destructor
  }

  h1Booking->fTitle = title;
  UpdateTitle(h1Booking->fTitle, unitName, fcnName);  

  return h1Booking;
}

//_____________________________________________________________________________
h1_booking* CreateH1Booking(const G4String& title,
                   const std::vector<G4double>& edges,
                   const G4String& unitName,
                   const G4String& fcnName)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);

  // Apply function
  std::vector<G4double> newEdges;
  ComputeEdges(edges, unit, fcn, newEdges);
  
  h1_booking* h1Booking = new h1_booking(newEdges); 
                    // h1_booking object is deleted in destructor

  h1Booking->fTitle = title;
  UpdateTitle(h1Booking->fTitle, unitName, fcnName);  
  
  return h1Booking;
}

//_____________________________________________________________________________
void UpdateH1Booking(h1_booking* h1Booking,
                     G4int nbins, G4double xmin, G4double xmax,  
                     const G4String& unitName,
                     const G4String& fcnName,
                     const G4String& binSchemeName)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  G4BinScheme binScheme = GetBinScheme(binSchemeName);

  if ( binScheme != kLogBinScheme ) {
    if ( binScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookH1Manager::SetH1",
                "Analysis_W013", JustWarning, description);
    }              
    h1Booking->fNbins = nbins;
    h1Booking->fXmin = fcn(xmin/unit);
    h1Booking->fXmax = fcn(xmax/unit);
  }
  else {
    // Compute edges
    ComputeEdges(nbins, xmin, xmax, unit, fcn, binScheme, h1Booking->fEdges);
  }

  UpdateTitle(h1Booking->fTitle, unitName, fcnName);  
}     

//_____________________________________________________________________________
void UpdateH1Booking(h1_booking* h1Booking,
                     const std::vector<G4double>& edges,
                     const G4String& unitName,
                     const G4String& fcnName)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);

  // Apply function
  ComputeEdges(edges, unit, fcn, h1Booking->fEdges);

  UpdateTitle(h1Booking->fTitle, unitName, fcnName);  
}     

//_____________________________________________________________________________
void ConfigureHbookH1(tools::hbook::h1* h1,
                      G4int nbins, G4double xmin, G4double xmax,  
                      const G4String& unitName,
                      const G4String& fcnName,
                      G4BinScheme binScheme)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);

  if ( binScheme != kLogBinScheme ) {
    if ( binScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookH1Manager::SetH1",
                "Analysis_W013", JustWarning, description);
    }              
    h1->configure(nbins, fcn(xmin/unit), fcn(xmax/unit));
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, unit, fcn, binScheme, edges);
    // Convert to float
    std::vector<float> fedges;
    ConvertToFloat(edges, fedges); 

    h1->configure(fedges);
  }
}     

//_____________________________________________________________________________
void ConfigureHbookH1(tools::hbook::h1* h1,
                      const std::vector<G4double>& edges,
                      const G4String& unitName,
                      const G4String& fcnName)
{
  // Apply function to edges
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, unit, fcn, newEdges);
  
  // Convert to float
  std::vector<float> newFEdges;
  ConvertToFloat(newEdges, newFEdges); 

  h1->configure(newFEdges);
}

}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookH1Manager::SetH1HbookIdOffset()
{
// Set  fH1HbookIdOffset if needed

  if ( fH1HbookIdOffset == -1 ) {
    if ( fFirstId > 0 ) 
      fH1HbookIdOffset = 0;
    else
      fH1HbookIdOffset = 1;
        
    if ( fH1HbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "H1 will be defined in HBOOK with ID = G4_firstHistoId + 1";
      G4Exception("ExG4HbookH1Manager::SetH1HbookIdOffset()",
                  "Analysis_W013", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookH1Manager::AddH1Information(const G4String& name,  
                                          const G4String& unitName, 
                                          const G4String& fcnName,
                                          G4BinScheme binScheme) const
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  fHnManager->AddH1Information(name, unitName, fcnName, unit, fcn, binScheme);
}  

//_____________________________________________________________________________
G4int ExG4HbookH1Manager::CreateH1FromBooking(h1_booking* h1Booking, 
                                              G4bool chDir)
{
// Create h1 from h1_booking.

  if ( chDir ) {
    // Go to histograms directory if defined
    if ( fFileManager->GetHistoDirectoryName() != "" ) {
      G4String histoPath = "//PAWC/LUN1/";
      histoPath.append(fFileManager->GetHistoDirectoryName().data());
      tools::hbook::CHCDIR(histoPath.data()," ");
    }
  }    

  G4int index = fH1Vector.size();
  G4int id = index + fFirstId;    
  G4HnInformation* 
    info = fHnManager->GetHnInformation(id, "CreateH1FromBooking");
  // Hbook index
  G4int hbookIndex = fH1HbookIdOffset + index + fFirstId;
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create from booking", "h1", info->GetName());
#endif

  // Create h1
  tools::hbook::h1* h1 = 0; 
  if ( ! h1Booking->fEdges.size() ) {
    h1 = new tools::hbook::h1(
               hbookIndex, h1Booking->fTitle, 
               h1Booking->fNbins, h1Booking->fXmin, h1Booking->fXmax);
  }
  else {               
    // Convert to float
    std::vector<float> newEdges;
    ConvertToFloat(h1Booking->fEdges, newEdges); 

    h1 = new tools::hbook::h1(hbookIndex, h1Booking->fTitle, newEdges);
  }
                           
  fH1Vector.push_back(h1);
  
  if ( chDir ) {
    if ( fFileManager->GetHistoDirectoryName() != "" ) {
      // Return to //PAWC/LUN1 :
      tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    }  
  }
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) { 
    G4ExceptionDescription description;
    description << " name : " << info->GetName() << " hbook index : " << hbookIndex; 
    fState.GetVerboseL3()->Message("create from booking", "h1", description);
  }  
#endif
  
  return id;
}  

//_____________________________________________________________________________
G4int ExG4HbookH1Manager::RegisterH1Booking(const G4String& name, 
                                            h1_booking* h1Booking)
{
  // Register h1
  G4int index = fH1BookingVector.size();  
  fH1BookingVector.push_back(h1Booking);
  fH1NameIdMap[name] = index + fFirstId;

  // Lock id
  fLockFirstId = true;

  return index + fFirstId;
}  

//_____________________________________________________________________________
void ExG4HbookH1Manager::BeginCreateH1(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H1", name);
#endif

  // Set  fH1HbookIdOffset if needed
  SetH1HbookIdOffset();
}

//_____________________________________________________________________________
G4int ExG4HbookH1Manager::FinishCreateH1(
                               const G4String& name, h1_booking* h1Booking,
                               const G4String& unitName, const G4String& fcnName,
                               G4BinScheme binScheme)
{
  // Register h1 booking
  G4int id = RegisterH1Booking(name, h1Booking);
  
  // Save H1 information
  AddH1Information(name, unitName, fcnName, binScheme);

  // Create h1 if the file is open
  if ( fFileManager->IsFile() ) {
    CreateH1FromBooking(h1Booking);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) { 
    G4int hbookIndex = fH1HbookIdOffset + id;
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << hbookIndex; 
    fState.GetVerboseL2()->Message("create", "H1", description);
  }  
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::BeginSetH1(
                               G4int id,
                               h1_booking* h1Booking,
                               G4HnInformation* info)
{                                
  h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH1()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  info = fHnManager->GetHnInformation(id,"SetH1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H1", info->GetName());
#endif

  return true;
}
  
//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::FinishSetH1(
                               G4int id,
                               G4HnInformation* info,
                               const G4String& unitName, 
                               const G4String& fcnName,
                               G4BinScheme binScheme)
{                                
  // Update information
  UpdateH1Information(info, unitName, fcnName, binScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
  
                                        
//_____________________________________________________________________________
void ExG4HbookH1Manager::CreateH1sFromBooking()
{
// Create all h1 from h1_booking.

  // Do nothing if any h1 histogram already exists
  // or no h1 histograms are booked
  if ( fH1Vector.size() || ( fH1BookingVector.size() == 0 ) ) return;       

  // Go to histograms directory if defined
  if ( fFileManager->GetHistoDirectoryName() != "" ) {
    G4String histoPath = "//PAWC/LUN1/";
    histoPath.append(fFileManager->GetHistoDirectoryName().data());
    tools::hbook::CHCDIR(histoPath.data()," ");
  }  

  // Create histograms
  std::vector<h1_booking*>::const_iterator it;
  for ( it = fH1BookingVector.begin(); it != fH1BookingVector.end(); ++it) {
    CreateH1FromBooking(*it, false);
  }  
  
  // Return backi from histograms directory if defined
  if ( fFileManager->GetHistoDirectoryName() != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
}   

//_____________________________________________________________________________
void ExG4HbookH1Manager::Reset()
{
// Reset histograms and ntuple  

  // Delete histograms
  std::vector<tools::hbook::h1*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete *it;
  }  

  // Clear vectors
  fH1Vector.clear();
}  
 
//_____________________________________________________________________________
h1_booking*  ExG4HbookH1Manager::GetH1Booking(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH1BookingVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "histo " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetH1Booking()",
                  "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }

  return fH1BookingVector[index];
}

//
// protected methods
//

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// (Only H1 implemented by now)
// According to the implementation by Michel Maire, originally in
// extended examples.

  // h1 histograms
  for ( G4int i=0; i<G4int(fH1Vector.size()); ++i ) {
    G4int id = i + fFirstId;
    G4HnInformation* info 
      = fHnManager->GetHnInformation(id, "WriteOnAscii"); 
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->GetAscii() ) continue; 
    tools::hbook::h1* h1 = fH1Vector[i];

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("write on ascii", "h1", info->GetName());
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
tools::hbook::h1*  ExG4HbookH1Manager::GetH1InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "ExG4HbookH1Manager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }
  
  // Do not return histogram if inactive 
  if ( fState.GetIsActivation() && onlyIfActive && ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fH1Vector[index];
}  
                                      
// 
// public methods
//

//_____________________________________________________________________________
G4int ExG4HbookH1Manager::CreateH1(
                               const G4String& name, const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName,
                               const G4String& binSchemeName)
{
  BeginCreateH1(name);

  G4BinScheme binScheme = GetBinScheme(binSchemeName);

  // Create h1 booking
  h1_booking* h1Booking 
    = CreateH1Booking(title, nbins, xmin, xmax, unitName, fcnName, binScheme);
    
  return FinishCreateH1(name, h1Booking, unitName, fcnName, binScheme); 
}                                         

//_____________________________________________________________________________
G4int ExG4HbookH1Manager::CreateH1(
                               const G4String& name, const G4String& title,
                               const std::vector<G4double>& edges,
                               const G4String& unitName,
                               const G4String& fcnName)
{                       
  BeginCreateH1(name);

  // Create h1 booking
  h1_booking* h1Booking 
    = CreateH1Booking(title, edges, unitName, fcnName);
    
  return FinishCreateH1(name, h1Booking, unitName, fcnName, kUserBinScheme); 
}                                         


//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1(G4int id,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, 
                               const G4String& fcnName,
                               const G4String& binSchemeName)
{                                
  h1_booking* h1Booking = 0;
  G4HnInformation* info = 0;

  if ( ! BeginSetH1(id, h1Booking, info) ) return false; 

  G4BinScheme binScheme = GetBinScheme(binSchemeName);

  // Update H1 booking
  UpdateH1Booking(h1Booking, 
                  nbins, xmin, xmax, unitName, fcnName, binScheme);

  // Re-configure histogram if it was already defined
  if ( fH1Vector.size() ) {
    tools::hbook::h1* h1 = GetH1(id);
    ConfigureHbookH1(h1, nbins, xmin, xmax, unitName, fcnName, binScheme);
  }  
  
  return FinishSetH1(id, info, unitName, fcnName, binScheme);
}
  
//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1(G4int id,
                               const std::vector<G4double>& edges,
                               const G4String& unitName, 
                               const G4String& fcnName)
{                                
  h1_booking* h1Booking = 0;
  G4HnInformation* info = 0;

  if ( ! BeginSetH1(id, h1Booking, info) ) return false; 

  // Update H1 booking
  UpdateH1Booking(h1Booking, edges, unitName, fcnName);

  // Re-configure histogram if it was already defined
  if ( fH1Vector.size() ) {
    tools::hbook::h1* h1 = GetH1(id);
    ConfigureHbookH1(h1, edges, unitName, fcnName);
  }  
  
  return 
    FinishSetH1(id, info, unitName, fcnName, kUserBinScheme);
}
  
//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::ScaleH1(G4int id, G4double factor)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "ScaleH1", false, false);
  if ( ! h1 ) return false;

  return h1->scale(factor);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::FillH1(G4int id, G4double value, G4double weight)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "FillH1", true, false);
  if ( ! h1 ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillH1 for " << id << G4endl; 
    return false; 
  }  

  G4HnDimensionInformation* info 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillH1");
  h1->fill(info->fFcn(value/info->fUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value 
                << " fcn(value/unit) " << info->fFcn(value/info->fUnit)  
                << " weight " << weight;
    fState.GetVerboseL4()->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
tools::hbook::h1*  ExG4HbookH1Manager::GetH1(G4int id, G4bool warn,
                                             G4bool onlyIfActive) const 
{
  return GetH1InFunction(id, "GetH1", warn, onlyIfActive);
}

//_____________________________________________________________________________
G4int  ExG4HbookH1Manager::GetH1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH1NameIdMap.find(name);
  if ( it ==  fH1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "ExG4HbookH1Manager::GetH1Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int ExG4HbookH1Manager::GetH1Nbins(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Nbins");
  if ( ! h1 ) return 0;
  
  return fBaseToolsManager.GetNbins(h1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookH1Manager::GetH1Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Xmin");
  if ( ! h1 ) return 0;
  
  return fBaseToolsManager.GetMin(h1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookH1Manager::GetH1Xmax(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1Xmax");
  if ( ! h1 ) return 0;
  
  return fBaseToolsManager.GetMax(h1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookH1Manager::GetH1Width(G4int id) const
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1XWidth", true, false);
  if ( ! h1 ) return 0;
  
  return fBaseToolsManager.GetWidth(h1->axis());
}  

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1Title(G4int id, const G4String& title)
{
  h1_booking* h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetH1Title()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  h1Booking->fTitle = title;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1XAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "SetH1XAxisTitle");
  if ( ! h1 ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h1, ExG4HbookBaseHnManager::kX, title);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1YAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "SetH1YAxisTitle");
  if ( ! h1 ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h1, ExG4HbookBaseHnManager::kY, title);
}  

//_____________________________________________________________________________
G4String ExG4HbookH1Manager::GetH1Title(G4int id) const
{
  h1_booking* h1Booking = GetH1Booking(id, false);
  if ( ! h1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "histogram " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetH1Title()",
                "Analysis_W011", JustWarning, description);
    return "";
  }
  
  return h1Booking->fTitle;
}  


//_____________________________________________________________________________
G4String ExG4HbookH1Manager::GetH1XAxisTitle(G4int id) const 
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1XAxisTitle");
  if ( ! h1 ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h1, ExG4HbookBaseHnManager::kX);
}  

//_____________________________________________________________________________
G4String ExG4HbookH1Manager::GetH1YAxisTitle(G4int id) const 
{
  tools::hbook::h1* h1 = GetH1InFunction(id, "GetH1YAxisTitle");
  if ( ! h1 ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h1, ExG4HbookBaseHnManager::kY);
}  

//_____________________________________________________________________________
G4bool ExG4HbookH1Manager::SetH1HbookIdOffset(G4int offset) 
{
  if ( fH1Vector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set H1HbookIdOffset as some H1 histogramms already exist.";
    G4Exception("G4HbookAnalysisManager::SetH1HbookIdOffset()",
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
  
  fH1HbookIdOffset = offset;
  return true;
}  

#endif
