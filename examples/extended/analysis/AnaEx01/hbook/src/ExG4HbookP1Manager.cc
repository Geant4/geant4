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
/// \file hbook/src/ExG4HbookP1Manager.cc
/// \brief Implementation of the ExG4HbookP1Manager class

// Author: Ivana Hrivnacova, 03/11/2014  (ivana@ipno.in2p3.fr)

#ifdef G4_USE_HBOOK

#include "ExG4HbookP1Manager.hh"
#include "ExG4HbookFileManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"
#include "G4BinScheme.hh"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
ExG4HbookP1Manager::ExG4HbookP1Manager(const G4AnalysisManagerState& state)
 : G4VP1Manager(state),
   fBaseToolsManager("P1"),
   fFileManager(0),
   fP1HbookIdOffset(-1),
   fP1Vector(),
   fP1BookingVector(),
   fP1NameIdMap()
{
}

//_____________________________________________________________________________
ExG4HbookP1Manager::~ExG4HbookP1Manager()
{  
  // Delete p1
  Reset();

  // Delete p1 booking 
  std::vector<p1_booking*>::iterator it;
  for ( it = fP1BookingVector.begin(); it != fP1BookingVector.end(); it++ ) {
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
void UpdateP1Information(G4HnInformation* hnInformation,
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  
  G4HnDimensionInformation* xInformation 
    = hnInformation->GetHnDimensionInformation(G4HnInformation::kX);
  xInformation->fUnitName = xunitName;
  xInformation->fFcnName = xfcnName;
  xInformation->fUnit = xunit;
  xInformation->fFcn = xfcn;
  xInformation->fBinScheme = xbinScheme;

  G4HnDimensionInformation* yInformation 
    = hnInformation->GetHnDimensionInformation(G4HnInformation::kY);
  yInformation->fUnitName = yunitName;
  yInformation->fFcnName = yfcnName;
  yInformation->fUnit = yunit;
  yInformation->fFcn = yfcn;
  yInformation->fBinScheme = kLinearBinScheme;
}  

//_____________________________________________________________________________
p1_booking* CreateP1Booking(const G4String& title,
                   G4int nbins, G4double xmin, G4double xmax,
                   G4double ymin, G4double ymax,
                   const G4String& xunitName,
                   const G4String& yunitName,
                   const G4String& xfcnName,
                   const G4String& yfcnName,
                   G4BinScheme xbinScheme)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  p1_booking* p1Booking = 0; 
  if ( xbinScheme != kLogBinScheme ) {
    if ( xbinScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookP1Manager::CreateP1",
                "Analysis_W013", JustWarning, description);
    }              
    p1Booking = new p1_booking(nbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                               yfcn(ymin/yunit), yfcn(ymax/yunit)); 
                    // p1_booking object is deleted in destructor
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    p1Booking = new p1_booking(edges, yfcn(ymin/yunit), yfcn(ymax/yunit)); 
                    // p1_booking object is deleted in destructor
  }

  p1Booking->fTitle = title;
  UpdateTitle(p1Booking->fTitle, xunitName, xfcnName);  

  return p1Booking;
}

//_____________________________________________________________________________
p1_booking* CreateP1Booking(const G4String& title,
                   const std::vector<G4double>& edges,
                   G4double ymin, G4double ymax,
                   const G4String& xunitName,
                   const G4String& yunitName,
                   const G4String& xfcnName,
                   const G4String& yfcnName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  // Apply function
  std::vector<G4double> newEdges;
  ComputeEdges(edges, xunit, xfcn, newEdges);
  
  p1_booking* p1Booking = new p1_booking(newEdges, yfcn(ymin/yunit), yfcn(ymax/yunit)); 
                    // p1_booking object is deleted in destructor

  p1Booking->fTitle = title;
  UpdateTitle(p1Booking->fTitle, xunitName, xfcnName);  
  
  return p1Booking;
}

//_____________________________________________________________________________
void UpdateP1Booking(p1_booking* p1Booking,
                     G4int nbins, G4double xmin, G4double xmax,  
                     G4double ymin, G4double ymax,  
                     const G4String& xunitName,
                     const G4String& yunitName,
                     const G4String& xfcnName,
                     const G4String& yfcnName,
                     const G4String& xbinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);

  if ( xbinScheme != kLogBinScheme ) {
    if ( xbinScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookP1Manager::SetP1",
                "Analysis_W013", JustWarning, description);
    }              
    p1Booking->fNbins = nbins;
    p1Booking->fXmin = xfcn(xmin/xunit);
    p1Booking->fXmax = xfcn(xmax/xunit);
    p1Booking->fYmin = yfcn(ymin/yunit);
    p1Booking->fYmax = yfcn(ymax/yunit);
  }
  else {
    // Compute edges
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, p1Booking->fEdges);
    p1Booking->fYmin = yfcn(ymin/yunit);
    p1Booking->fYmax = yfcn(ymax/yunit);
  }

  UpdateTitle(p1Booking->fTitle, xunitName, xfcnName);  
}     

//_____________________________________________________________________________
void UpdateP1Booking(p1_booking* p1Booking,
                     const std::vector<G4double>& edges,
                     G4double ymin, G4double ymax,  
                     const G4String& xunitName,
                     const G4String& yunitName,
                     const G4String& xfcnName,
                     const G4String& yfcnName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  // Apply function
  ComputeEdges(edges, xunit, xfcn, p1Booking->fEdges);
  p1Booking->fYmin = yfcn(ymin/yunit);
  p1Booking->fYmax = yfcn(ymax/yunit);

  UpdateTitle(p1Booking->fTitle, xunitName, xfcnName);  
}     

//_____________________________________________________________________________
void ConfigureHbookP1(tools::hbook::p1* p1,
                      G4int nbins, G4double xmin, G4double xmax,  
                      G4double ymin, G4double ymax,  
                      const G4String& xunitName,
                      const G4String& yunitName,
                      const G4String& xfcnName,
                      const G4String& yfcnName,
                      G4BinScheme xbinScheme)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);

  if ( xbinScheme != kLogBinScheme ) {
    if ( xbinScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("ExG4HbookP1Manager::SetP1",
                "Analysis_W013", JustWarning, description);
    }              
    // not available !!              
    p1->configure(nbins, xfcn(xmin/xunit), xfcn(xmax/xunit),
                  yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    // Convert to float
    std::vector<float> fedges;
    ConvertToFloat(edges, fedges); 

    // not available !!              
    //p1->configure(edges,  yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
}     

//_____________________________________________________________________________
void ConfigureHbookP1(tools::hbook::p1* /*p1*/,
                      const std::vector<G4double>& edges,
                       G4double /*ymin*/, G4double /*ymax*/,  
                      const G4String& xunitName,
                      const G4String& /*yunitName*/,
                      const G4String& xfcnName,
                      const G4String& /*yfcnName*/)
{
  // Apply function to edges
  G4double xunit = GetUnitValue(xunitName);
  //G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  //G4Fcn yfcn = GetFunction(yfcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, xunit, xfcn, newEdges);
  
  // Convert to float
  std::vector<float> newFEdges;
  ConvertToFloat(newEdges, newFEdges); 

  // not available !!              
  //p1->configure(newFEdges);
}

}

// 
// private methods
//

//_____________________________________________________________________________
void ExG4HbookP1Manager::SetP1HbookIdOffset()
{
// Set  fP1HbookIdOffset if needed

  if ( fP1HbookIdOffset == -1 ) {
    if ( fFirstId > 0 ) 
      fP1HbookIdOffset = 0;
    else
      fP1HbookIdOffset = 1;
        
    if ( fP1HbookIdOffset > 0 ) {
      G4ExceptionDescription description;
      description << "P1 will be defined in HBOOK with ID = G4_firstProfileId + 1";
      G4Exception("ExG4HbookP1Manager::SetP1HbookIdOffset()",
                  "Analysis_W013", JustWarning, description);
    }              
  }
}  

//_____________________________________________________________________________
void ExG4HbookP1Manager::AddP1Information(const G4String& name,  
                            const G4String& xunitName, 
                            const G4String& yunitName, 
                            const G4String& xfcnName,
                            const G4String& yfcnName,
                            G4BinScheme xbinScheme) const
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  fHnManager
    ->AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                       xunit, yunit, xfcn, yfcn, 
                       xbinScheme, xbinScheme);
}  

//_____________________________________________________________________________
G4int ExG4HbookP1Manager::CreateP1FromBooking(p1_booking* p1Booking, 
                                              G4bool chDir)
{
// Create p1 from p1_booking.

  if ( chDir ) {
    // Go to profiles directory if defined
    if ( fFileManager->GetProfileDirectoryName() != "" ) {
      G4String profilePath = "//PAWC/LUN1/";
      profilePath.append(fFileManager->GetProfileDirectoryName().data());
      tools::hbook::CHCDIR(profilePath.data()," ");
    }
  }    

  G4int index = fP1Vector.size();
  G4int id = index + fFirstId;    
  G4HnInformation* 
    info = fHnManager->GetHnInformation(id, "CreateP1FromBooking");
  // Hbook index
  G4int hbookIndex = fP1HbookIdOffset + index + fFirstId;
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create from booking", "p1", info->GetName());
#endif

  // Create p1
  tools::hbook::p1* p1 = 0; 
  if ( ! p1Booking->fEdges.size() ) {
    p1 = new tools::hbook::p1(
               hbookIndex, p1Booking->fTitle, 
               p1Booking->fNbins, p1Booking->fXmin, p1Booking->fXmax, 
               p1Booking->fYmin, p1Booking->fYmax);
  }
  else {               
    // Convert to float
    std::vector<float> newEdges;
    ConvertToFloat(p1Booking->fEdges, newEdges); 

    // not supported
    //p1 = new tools::hbook::p1(hbookIndex, p1Booking->fTitle, newEdges,
    //                          p1Booking->fYmin, p1Booking->fYmax);
  }
                           
  fP1Vector.push_back(p1);
  
  if ( chDir ) {
    if ( fFileManager->GetProfileDirectoryName() != "" ) {
      // Return to //PAWC/LUN1 :
      tools::hbook::CHCDIR("//PAWC/LUN1"," ");
    }  
  }
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL3() ) { 
    G4ExceptionDescription description;
    description << " name : " << info->GetName() << " hbook index : " << hbookIndex; 
    fState.GetVerboseL3()->Message("create from booking", "p1", description);
  }  
#endif
  
  return id;
}  

//_____________________________________________________________________________
G4int ExG4HbookP1Manager::RegisterP1Booking(const G4String& name, 
                                            p1_booking* p1Booking)
{
  // Register p1
  G4int index = fP1BookingVector.size();  
  fP1BookingVector.push_back(p1Booking);
  fP1NameIdMap[name] = index + fFirstId;

  // Lock id
  fLockFirstId = true;

  return index + fFirstId;
}  

//_____________________________________________________________________________
void ExG4HbookP1Manager::BeginCreateP1(const G4String& name)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "P1", name);
#endif

  // Set  fP1HbookIdOffset if needed
  SetP1HbookIdOffset();
}

//_____________________________________________________________________________
G4int ExG4HbookP1Manager::FinishCreateP1(
                               const G4String& name, p1_booking* p1Booking,
                               const G4String& xunitName, 
                               const G4String& yunitName, 
                               const G4String& xfcnName,
                               const G4String& yfcnName,
                               G4BinScheme xbinScheme)
{
  // Register p1 booking
  G4int id = RegisterP1Booking(name, p1Booking);
  
  // Save P1 information
  AddP1Information(name, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Create p1 if the file is open
  if ( fFileManager->IsFile() ) {
    CreateP1FromBooking(p1Booking);
  }

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) { 
    G4int hbookIndex = fP1HbookIdOffset + id;
    G4ExceptionDescription description;
    description << " name : " << name << " hbook index : " << hbookIndex; 
    fState.GetVerboseL2()->Message("create", "P1", description);
  }  
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::BeginSetP1(
                               G4int id,
                               p1_booking* p1Booking,
                               G4HnInformation* info)
{                                
  p1Booking = GetP1Booking(id, false);
  if ( ! p1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "profile " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetP1()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  info = fHnManager->GetHnInformation(id,"SetP1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "P1", info->GetName());
#endif

  return true;
}
  
//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::FinishSetP1(
                               G4int id,
                               G4HnInformation* info,
                               const G4String& xunitName, 
                               const G4String& yunitName, 
                               const G4String& xfcnName,
                               const G4String& yfcnName,
                               G4BinScheme xbinScheme)
{                                
  // Update information
  UpdateP1Information(info, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
  
                                        
//_____________________________________________________________________________
void ExG4HbookP1Manager::CreateP1sFromBooking()
{
// Create all p1 from p1_booking.

  // Do nothing if any p1 profile already exists
  // or no p1 profiles are booked
  if ( fP1Vector.size() || ( fP1BookingVector.size() == 0 ) ) return;       

  // Go to profiles directory if defined
  if ( fFileManager->GetProfileDirectoryName() != "" ) {
    G4String profilePath = "//PAWC/LUN1/";
    profilePath.append(fFileManager->GetProfileDirectoryName().data());
    tools::hbook::CHCDIR(profilePath.data()," ");
  }  

  // Create profiles
  std::vector<p1_booking*>::const_iterator it;
  for ( it = fP1BookingVector.begin(); it != fP1BookingVector.end(); ++it) {
    CreateP1FromBooking(*it, false);
  }  
  
  // Return backi from profiles directory if defined
  if ( fFileManager->GetProfileDirectoryName() != "" ) {
    // Return to //PAWC/LUN1 :
    tools::hbook::CHCDIR("//PAWC/LUN1"," ");
  }  
}   

//_____________________________________________________________________________
void ExG4HbookP1Manager::Reset()
{
// Reset profiles and ntuple  

  // Delete profiles
  std::vector<tools::hbook::p1*>::iterator it;
  for (it = fP1Vector.begin(); it != fP1Vector.end(); it++ ) {
    delete *it;
  }  

  // Clear vectors
  fP1Vector.clear();
}  
 
//_____________________________________________________________________________
p1_booking*  ExG4HbookP1Manager::GetP1Booking(G4int id, G4bool warn) const 
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fP1BookingVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "profile " << id << " does not exist.";
      G4Exception("G4HbookAnalysisManager::GetP1Booking()",
                  "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }

  return fP1BookingVector[index];
}

//
// protected methods
//
/*
//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// (Only P1 implemented by now)
// According to the implementation by Michel Maire, originally in
// extended examples.

  // p1 profiles
  for ( G4int i=0; i<G4int(fP1Vector.size()); ++i ) {
    G4int id = i + fFirstId;
    G4HnInformation* info 
      = fHnManager->GetHnInformation(id, "WriteOnAscii"); 
    // skip writing if activation is enabled and P1 is inactivated
    if ( ! info->fAscii ) continue; 
    tools::hbook::p1* p1 = fP1Vector[i];

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("write on ascii", "p1", info->GetName());
#endif
  
    output << "\n  1D profile " << id << ": " << p1->title() 
           << "\n \n \t     X \t\t     Y" << G4endl;
    
    for (G4int j=0; j< G4int(p1->axis().bins()); ++j) {
       output << "  " << j << "\t" 
              << p1->axis().bin_center(j) << "\t"
              << p1->bin_height(j) << G4endl;
    } 
  }
  
  return true;
}  
*/
//_____________________________________________________________________________
tools::hbook::p1*  ExG4HbookP1Manager::GetP1InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fP1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "ExG4HbookP1Manager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "profile " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }
  
  // Do not return profile if inactive 
  if ( fState.GetIsActivation() && onlyIfActive && ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fP1Vector[index];
}  
                                      
// 
// public methods
//

//_____________________________________________________________________________
G4int ExG4HbookP1Manager::CreateP1(
                               const G4String& name, const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, 
                               const G4String& yunitName, 
                               const G4String& xfcnName,
                               const G4String& yfcnName,
                               const G4String& xbinSchemeName)
{
  BeginCreateP1(name);

  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);

  // Create p1 booking
  p1_booking* p1Booking 
    = CreateP1Booking(title, nbins, xmin, xmax, ymin, ymax, 
                      xunitName, yunitName, xfcnName, yfcnName, xbinScheme);
    
  return FinishCreateP1(name, p1Booking, 
                        xunitName, yunitName, xfcnName, yfcnName, xbinScheme); 
}                                         

//_____________________________________________________________________________
G4int ExG4HbookP1Manager::CreateP1(
                               const G4String& name, const G4String& title,
                               const std::vector<G4double>& edges,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, 
                               const G4String& yunitName, 
                               const G4String& xfcnName,
                               const G4String& yfcnName)
{                       
  BeginCreateP1(name);

  // Create p1 booking
  p1_booking* p1Booking 
    = CreateP1Booking(title, edges, ymin, ymax,
                      xunitName, yunitName, xfcnName, yfcnName);
    
  return FinishCreateP1(name, p1Booking,
                        xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme); 
}                                         


//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1(G4int id,
                               G4int nbins, G4double xmin, G4double xmax,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName,
                               const G4String& xbinSchemeName)
{                                
  p1_booking* p1Booking = 0;
  G4HnInformation* info = 0;

  if ( ! BeginSetP1(id, p1Booking, info) ) return false; 

  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);

  // Update P1 booking
  UpdateP1Booking(p1Booking, nbins, xmin, xmax, ymin, ymax, 
                  xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Re-configure profile if it was already defined
  if ( fP1Vector.size() ) {
    tools::hbook::p1* p1 = GetP1(id);
    ConfigureHbookP1(p1, nbins, xmin, xmax, ymin, ymax, 
                     xunitName, yunitName, xfcnName, yfcnName, xbinScheme);
  }  
  
  return FinishSetP1(id, info, 
                     xunitName, yunitName, xfcnName, yfcnName, xbinScheme);
}
  
//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1(G4int id,
                               const std::vector<G4double>& edges,
                               G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
{                                
  p1_booking* p1Booking = 0;
  G4HnInformation* info = 0;

  if ( ! BeginSetP1(id, p1Booking, info) ) return false; 

  // Update P1 booking
  UpdateP1Booking(p1Booking, edges, ymin, ymax, 
                  xunitName, yunitName, xfcnName, yfcnName);

  // Re-configure profile if it was already defined
  if ( fP1Vector.size() ) {
    tools::hbook::p1* p1 = GetP1(id);
    ConfigureHbookP1(p1, edges, ymin, ymax,
                     xunitName, yunitName, xfcnName, yfcnName);
  }  
  
  return 
    FinishSetP1(id, info, 
                xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme);
}
  
//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::ScaleP1(G4int id, G4double factor)
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "ScaleP1", false, false);
  if ( ! p1 ) return false;

  return p1->scale(factor);
}  

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::FillP1(G4int id, G4double xvalue, G4double yvalue,
                                  G4double weight)
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "FillP1", true, false);
  if ( ! p1 ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillP1 for " << id << G4endl; 
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillP1");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillP1");

  p1->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
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
    fState.GetVerboseL4()->Message("fill", "P1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
tools::hbook::p1*  ExG4HbookP1Manager::GetP1(G4int id, G4bool warn,
                                             G4bool onlyIfActive) const 
{
  return GetP1InFunction(id, "GetP1", warn, onlyIfActive);
}

//_____________________________________________________________________________
G4int  ExG4HbookP1Manager::GetP1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fP1NameIdMap.find(name);
  if ( it ==  fP1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "ExG4HbookP1Manager::GetP1Id";
      G4ExceptionDescription description;
      description << "      " << "profile " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int ExG4HbookP1Manager::GetP1Nbins(G4int id) const
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1Nbins");
  if ( ! p1 ) return 0;
  
  return fBaseToolsManager.GetNbins(p1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookP1Manager::GetP1Xmin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1Xmin");
  if ( ! p1 ) return 0;
  
  return fBaseToolsManager.GetMin(p1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookP1Manager::GetP1Xmax(G4int id) const
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1Xmax");
  if ( ! p1 ) return 0;
  
  return fBaseToolsManager.GetMax(p1->axis());
}  

//_____________________________________________________________________________
G4double ExG4HbookP1Manager::GetP1XWidth(G4int id) const
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1XWidth", true, false);
  if ( ! p1 ) return 0;
  
  G4int nbins = p1->axis().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for p1 id = " << id << ").";
    G4Exception("ExG4HbookP1Manager::GetP1Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  return fBaseToolsManager.GetWidth(p1->axis());
}

//_____________________________________________________________________________
G4double ExG4HbookP1Manager::GetP1Ymin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1Ymin", true, false);
  if ( ! p1 ) return 0;
  
  // not available
  //return p1->min_v();

  G4String inFunction = "ExG4HbookP1Manager::GetP1Ymin";
  G4ExceptionDescription description;
  description << "Get function not available.";
  G4Exception(inFunction, "Analysis_W011", JustWarning, description);
  return 0;
}

//_____________________________________________________________________________
G4double ExG4HbookP1Manager::GetP1Ymax(G4int id) const
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1Ymax", true, false);
  if ( ! p1 ) return 0;
  
  // not available
  //return p1->max_v();

  G4String inFunction = "ExG4HbookP1Manager::GetP1Ymax";
  G4ExceptionDescription description;
  description << "Get function not available.";
  G4Exception(inFunction, "Analysis_W011", JustWarning, description);
  return 0;
}

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1Title(G4int id, const G4String& title)
{
  p1_booking* p1Booking = GetP1Booking(id, false);
  if ( ! p1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "profile " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::SetP1Title()",
                "Analysis_W011", JustWarning, description);
    return false;
  }

  p1Booking->fTitle = title;
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1XAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "SetP1XAxisTitle");
  if ( ! p1 ) return false;
  
  p1->add_annotation(tools::hbook::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1YAxisTitle(G4int id, const G4String& title)
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "SetP1YAxisTitle");
  if ( ! p1 ) return false;
  
  p1->add_annotation(tools::hbook::key_axis_y_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4String ExG4HbookP1Manager::GetP1Title(G4int id) const
{
  p1_booking* p1Booking = GetP1Booking(id, false);
  if ( ! p1Booking ) {
    G4ExceptionDescription description;
    description << "      " << "profile " << id << " does not exist.";
    G4Exception("G4HbookAnalysisManager::GetP1Title()",
                "Analysis_W011", JustWarning, description);
    return "";
  }
  
  return p1Booking->fTitle;
}  


//_____________________________________________________________________________
G4String ExG4HbookP1Manager::GetP1XAxisTitle(G4int id) const 
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1XAxisTitle");
  if ( ! p1 ) return "";
  
  G4String title;
  G4bool result = p1->annotation(tools::hbook::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for p1 id = " << id << ").";
    G4Exception("ExG4HbookP1Manager::GetP1XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String ExG4HbookP1Manager::GetP1YAxisTitle(G4int id) const 
{
  tools::hbook::p1* p1 = GetP1InFunction(id, "GetP1YAxisTitle");
  if ( ! p1 ) return "";
  
  G4String title;
  G4bool result = p1->annotation(tools::hbook::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for p1 id = " << id << ").";
    G4Exception("ExG4HbookP1Manager::GetP1YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4bool ExG4HbookP1Manager::SetP1HbookIdOffset(G4int offset) 
{
  if ( fP1Vector.size() ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set P1HbookIdOffset as some P1 profiles already exist.";
    G4Exception("G4HbookAnalysisManager::SetP1HbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }
  
  if ( fFirstId + offset < 1 ) {
    G4ExceptionDescription description;
    description << "The first profile HBOOK id must be >= 1.";
    G4Exception("G4HbookAnalysisManager::SetP1HbookIdOffset()",
                 "Analysis_W013", JustWarning, description);
    return false;             
  }
  
  fP1HbookIdOffset = offset;
  return true;
}  

#endif
