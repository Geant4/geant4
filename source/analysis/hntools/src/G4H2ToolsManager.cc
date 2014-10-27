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
// $Id: G4H2ToolsManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013 (ivana@ipno.in2p3.fr)

#include "G4H2ToolsManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h2d"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4H2ToolsManager::G4H2ToolsManager(const G4AnalysisManagerState& state)
 : G4VH2Manager(state),
   fBaseToolsManager("H2"),
   fH2Vector(), 
   fH2NameIdMap()
{
}

//_____________________________________________________________________________
G4H2ToolsManager::~G4H2ToolsManager()
{  
  std::vector<tools::histo::h2d*>::iterator it;
  for (it = fH2Vector.begin(); it != fH2Vector.end(); it++ ) {
    delete (*it);
  }  
}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void UpdateH2Information(G4HnInformation* hnInformation,
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme)
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
  yInformation->fBinScheme = ybinScheme;
}  
                           
//_____________________________________________________________________________
void AddH2Annotation(tools::histo::h2d* h2d,
                     const G4String& xunitName, 
                     const G4String& yunitName, 
                     const G4String& xfcnName,
                     const G4String& yfcnName)
{                          
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  h2d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
}               
                          
//_____________________________________________________________________________
tools::histo::h2d* CreateToolsH2(
                         const G4String& title,
                         G4int nxbins, G4double xmin, G4double xmax,
                         G4int nybins, G4double ymin, G4double ymax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& xfcnName, 
                         const G4String& yfcnName,
                         const G4String& xbinSchemeName,
                         const G4String& ybinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  
  if ( xbinScheme != kLogBinScheme && ybinScheme !=  kLogBinScheme) {
    if ( xbinScheme == kUserBinScheme || ybinScheme == kUserBinScheme) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H2ToolsManager::CreateH2",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::h2d(title, 
                                 nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                                 nybins, yfcn(ymin/yunit), yfcn(ymax/yunit));
               // h2 objects are deleted in destructor and reset when 
               // closing a file.
  }
  else {
    // Compute edges
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    return new tools::histo::h2d(title, xedges, yedges); 
  }
}     

//_____________________________________________________________________________
tools::histo::h2d* CreateToolsH2(
                         const G4String& title,
                         const std::vector<G4double>& xedges,
                         const std::vector<G4double>& yedges,
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
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);
  
  return new tools::histo::h2d(title, xnewEdges, ynewEdges); 
             // h2 objects are deleted in destructor and reset when 
             // closing a file.
}  

//_____________________________________________________________________________
void  ConfigureToolsH2(tools::histo::h2d* h2d,
                       G4int nxbins, G4double xmin, G4double xmax,
                       G4int nybins, G4double ymin, G4double ymax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& xfcnName, 
                       const G4String& yfcnName,
                       const G4String& xbinSchemeName,
                       const G4String& ybinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  
  if ( xbinScheme != kLogBinScheme && ybinScheme !=  kLogBinScheme) {
    if ( xbinScheme == kUserBinScheme || ybinScheme == kUserBinScheme) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H2ToolsManager::CreateH2",
                "Analysis_W013", JustWarning, description);
    }              
    h2d->configure(nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                   nybins, yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
  else {
    // Compute bins
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    h2d->configure(xedges, yedges);
  }
}     

//_____________________________________________________________________________
void  ConfigureToolsH2(tools::histo::h2d* h2d,
                       const std::vector<G4double>& xedges,
                       const std::vector<G4double>& yedges,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);

  G4double yunit = GetUnitValue(yunitName);
  G4Fcn yfcn = GetFunction(yfcnName);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);

  h2d->configure(xnewEdges, ynewEdges);
}

}


//
// private methods
//

//_____________________________________________________________________________
tools::histo::h2d*  G4H2ToolsManager::GetH2InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH2Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4H2ToolsManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }

  // Do not return histogram if inactive 
  if ( fState.GetIsActivation()  && onlyIfActive && 
       ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fH2Vector[index];
}
  

//_____________________________________________________________________________
void G4H2ToolsManager::AddH2Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme) const
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  fHnManager
    ->AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                       xunit, yunit, xfcn, yfcn, 
                       xbinScheme, ybinScheme);
}  

//_____________________________________________________________________________
G4int G4H2ToolsManager::RegisterToolsH2(tools::histo::h2d* h2d, 
                          const G4String& name)
{
  G4int index = fH2Vector.size();
  fH2Vector.push_back(h2d);
  
  fLockFirstId = true;
  fH2NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H2ToolsManager::CreateH2(const G4String& name,  const G4String& title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& xbinSchemeName, 
                          const G4String& ybinSchemeName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H2", name);
#endif
  tools::histo::h2d* h2d
    = CreateToolsH2(title, nxbins, xmin, xmax, nybins, ymin, ymax, 
                    xunitName, yunitName, xfcnName, yfcnName, 
                    xbinSchemeName, ybinSchemeName);
    
  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Save H2 information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  AddH2Information(
    name, xunitName, yunitName, xfcnName, yfcnName, xbinScheme, ybinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH2(h2d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H2", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4int G4H2ToolsManager::CreateH2(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& xedges,
                          const std::vector<G4double>& yedges,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H2", name);
#endif
  tools::histo::h2d* h2d
    = CreateToolsH2(title, xedges, yedges, 
                    xunitName, yunitName, xfcnName, yfcnName); 
    
  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Save H2 information
  AddH2Information(
    name, xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme, kUserBinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH2(h2d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H2", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2(G4int id,
                            G4int nxbins, G4double xmin, G4double xmax, 
                            G4int nybins, G4double ymin, G4double ymax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& xbinSchemeName, 
                            const G4String& ybinSchemeName)
{                                
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetH2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H2", info->GetName());
#endif

  // Configure tools h2
  ConfigureToolsH2(
    h2d, nxbins, xmin, xmax, nybins, ymin, ymax, 
    xunitName, yunitName, xfcnName, yfcnName, xbinSchemeName, ybinSchemeName);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Update information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  UpdateH2Information(
    info, xunitName, yunitName, xfcnName, yfcnName, xbinScheme, ybinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2(G4int id,
                            const std::vector<G4double>& xedges,
                            const std::vector<G4double>& yedges,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName)
{                                
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetH2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H2", info->GetName());
#endif

  // Configure tools h2
  ConfigureToolsH2(h2d, xedges, yedges, xunitName, yunitName, xfcnName, yfcnName);

  // Add annotation
  AddH2Annotation(h2d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Update information
  UpdateH2Information(
    info, xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme, kUserBinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4H2ToolsManager::ScaleH2(G4int id, G4double factor)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "ScaleH2", false, false);
  if ( ! h2d ) return false;

  return h2d->scale(factor);
}  
                           
//_____________________________________________________________________________
G4bool G4H2ToolsManager::FillH2(G4int id, 
                                     G4double xvalue, G4double yvalue, 
                                     G4double weight)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "FillH2", true, false);
  if ( ! h2d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillH2");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillH2");

  h2d->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
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
G4int  G4H2ToolsManager::GetH2Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH2NameIdMap.find(name);
  if ( it ==  fH2NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4H2ToolsManager::GetH2Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return kInvalidId;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nxbins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NXbins");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetNbins(*h2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmin");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetMin(*h2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmax");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetMax(*h2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2XWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XWidth", true, false);
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetWidth(*h2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nybins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NYbins");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetNbins(*h2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymin");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetMin(*h2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymax");
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetMax(*h2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2YWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YWidth", true, false);
  if ( ! h2d ) return 0;
  
  return fBaseToolsManager.GetWidth(*h2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2Title(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2Title");
  if ( ! h2d ) return false;
  
  return fBaseToolsManager.SetTitle(*h2d, title);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2XAxisTitle");
  if ( ! h2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2d, G4BaseToolsManager::kX, title);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2YAxisTitle");
  if ( ! h2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2d, G4BaseToolsManager::kY, title);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2ZAxisTitle");
  if ( ! h2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h2d, G4BaseToolsManager::kZ, title);
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2Title(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Title");
  if ( ! h2d ) return "";
  
  return fBaseToolsManager.GetTitle(*h2d);
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2XAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XAxisTitle");
  if ( ! h2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2d, G4BaseToolsManager::kX);
} 

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2YAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YAxisTitle");
  if ( ! h2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2ZAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2ZAxisTitle");
  if ( ! h2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h2d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.
// Not yet available for H2

  return ! fHnManager->IsAscii();
} 

//
// public methods
// 

//_____________________________________________________________________________
G4int G4H2ToolsManager::AddH2(const G4String& name, tools::histo::h2d* h2d)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("add", "H2", name);
#endif

  // Add annotation
  AddH2Annotation(h2d, "none", "none", "none", "none");        
  // Add information
  AddH2Information(name, "none", "none", "none", "none", 
                   kLinearBinScheme, kLinearBinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH2(h2d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("add", "H2", name);
#endif
  return id;
}  

//_____________________________________________________________________________
void G4H2ToolsManager::AddH2Vector(
                          const std::vector<tools::histo::h2d*>& h2Vector)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("merge", "all h2", "");
#endif
  std::vector<tools::histo::h2d*>::const_iterator itw = h2Vector.begin();
  std::vector<tools::histo::h2d*>::iterator it;
  for (it = fH2Vector.begin(); it != fH2Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("merge", "all h2", "");
#endif
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::h2d*>::iterator it;
  for (it = fH2Vector.begin(); it != fH2Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::IsEmpty() const
{
  return ! fH2Vector.size();
}  
 
//_____________________________________________________________________________
tools::histo::h2d*  G4H2ToolsManager::GetH2(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const 
{
  return GetH2InFunction(id, "GetH2", warn, onlyIfActive);
}

