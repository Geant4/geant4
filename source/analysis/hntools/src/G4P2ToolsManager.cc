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

// Author: Ivana Hrivnacova, 24/07/2014  (ivana@ipno.in2p3.fr)

#include "G4P2ToolsManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/p2d"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4P2ToolsManager::G4P2ToolsManager(const G4AnalysisManagerState& state)
 : G4VP2Manager(state),
   fBaseToolsManager("P2"),
   fP2Vector(), 
   fP2NameIdMap()
{
}

//_____________________________________________________________________________
G4P2ToolsManager::~G4P2ToolsManager()
{  
  std::vector<tools::histo::p2d*>::iterator it;
  for (it = fP2Vector.begin(); it != fP2Vector.end(); it++ ) {
    delete (*it);
  }  
}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void UpdateP2Information(G4HnInformation* hnInformation,
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);

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

  G4HnDimensionInformation* zInformation 
    = hnInformation->GetHnDimensionInformation(G4HnInformation::kZ);
  zInformation->fUnitName = zunitName;
  zInformation->fFcnName = zfcnName;
  zInformation->fUnit = zunit;
  zInformation->fFcn = zfcn;
  zInformation->fBinScheme = kLinearBinScheme;
}  
                           
//_____________________________________________________________________________
void AddP2Annotation(tools::histo::p2d* p2d,
                     const G4String& xunitName, 
                     const G4String& yunitName, 
                     const G4String& zunitName, 
                     const G4String& xfcnName,
                     const G4String& yfcnName,
                     const G4String& zfcnName)
{                          
  G4String xaxisTitle;
  G4String yaxisTitle;
  G4String zaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  UpdateTitle(zaxisTitle, zunitName, zfcnName);        
  p2d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  p2d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
  p2d->add_annotation(tools::histo::key_axis_z_title(), zaxisTitle);
}               
                          
//_____________________________________________________________________________
tools::histo::p2d* CreateToolsP2(
                         const G4String& title,
                         G4int nxbins, G4double xmin, G4double xmax,
                         G4int nybins, G4double ymin, G4double ymax,
                         G4double zmin, G4double zmax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& zunitName,
                         const G4String& xfcnName, 
                         const G4String& zfcnName,
                         const G4String& yfcnName,
                         const G4String& xbinSchemeName,
                         const G4String& ybinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  
  if ( xbinScheme != kLogBinScheme && ybinScheme !=  kLogBinScheme) {
    if ( xbinScheme == kUserBinScheme || ybinScheme == kUserBinScheme ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4P2ToolsManager::CreateP2",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::p2d(title, 
                                 nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                                 nybins, yfcn(ymin/yunit), yfcn(ymax/yunit),
                                 zfcn(zmin/zunit), zfcn(zmax/zunit));
               // p2 objects are deleted in destructor and reset when 
               // closing a file.
  }
  else {
    // Compute edges
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    return new tools::histo::p2d(title, xedges, yedges, 
                                 zfcn(zmin/zunit), zfcn(zmax/zunit)); 
  }
}     

//_____________________________________________________________________________
tools::histo::p2d* CreateToolsP2(
                         const G4String& title,
                         const std::vector<G4double>& xedges,
                         const std::vector<G4double>& yedges,
                         G4double zmin, G4double zmax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& zunitName,
                         const G4String& xfcnName,
                         const G4String& yfcnName,
                         const G4String& zfcnName)
{                          
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);

  // Apply function 
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);
  
  return new tools::histo::p2d(title, xnewEdges, ynewEdges, 
                               zfcn(zmin/zunit), zfcn(zmax/zunit)); 
             // p2 objects are deleted in destructor and reset when 
             // closing a file.
}  

//_____________________________________________________________________________
void  ConfigureToolsP2(tools::histo::p2d* p2d,
                       G4int nxbins, G4double xmin, G4double xmax,
                       G4int nybins, G4double ymin, G4double ymax,
                       G4double zmin, G4double zmax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& zunitName,
                       const G4String& xfcnName, 
                       const G4String& yfcnName,
                       const G4String& zfcnName,
                       const G4String& xbinSchemeName,
                       const G4String& ybinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);
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
      G4Exception("G4P2ToolsManager::CreateP2",
                "Analysis_W013", JustWarning, description);
    }              
    p2d->configure(nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                   nybins, yfcn(ymin/yunit), yfcn(ymax/yunit),
                   zfcn(zmin/zunit), zfcn(zmax/zunit));
  }
  else {
    // Compute bins
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    p2d->configure(xedges, yedges, zfcn(zmin/zunit), zfcn(zmax/zunit));
  }
}     

//_____________________________________________________________________________
void  ConfigureToolsP2(tools::histo::p2d* p2d,
                       const std::vector<G4double>& xedges,
                       const std::vector<G4double>& yedges,
                       G4double zmin, G4double zmax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& zunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName,
                       const G4String& zfcnName)
{
  // Apply function to edges
  G4double xunit = GetUnitValue(xunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);

  G4double yunit = GetUnitValue(yunitName);
  G4Fcn yfcn = GetFunction(yfcnName);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);

  G4double zunit = GetUnitValue(zunitName);
  G4Fcn zfcn = GetFunction(zfcnName);
  p2d->configure(xnewEdges, ynewEdges, zfcn(zmin/zunit), zfcn(zmax/zunit));
}

}


//
// private methods
//

//_____________________________________________________________________________
tools::histo::p2d*  G4P2ToolsManager::GetP2InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fP2Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4P2ToolsManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "profile " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return 0;         
  }

  // Do not return profile if inactive 
  if ( fState.GetIsActivation()  && onlyIfActive && 
       ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fP2Vector[index];
}
  

//_____________________________________________________________________________
void G4P2ToolsManager::AddP2Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme) const
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);
  fHnManager
    ->AddH3Information(name, xunitName, yunitName, zunitName,
                       xfcnName, yfcnName, zfcnName,
                       xunit, yunit, zunit, xfcn, yfcn, zfcn,
                       xbinScheme, ybinScheme, kLinearBinScheme);
}  

//_____________________________________________________________________________
G4int G4P2ToolsManager::RegisterToolsP2(tools::histo::p2d* p2d, 
                          const G4String& name)
{
  G4int index = fP2Vector.size();
  fP2Vector.push_back(p2d);
  
  fLockFirstId = true;
  fP2NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4P2ToolsManager::CreateP2(const G4String& name,  const G4String& title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax,
                          G4double zmin, G4double zmax,
                          const G4String& xunitName, const G4String& yunitName, 
                          const G4String& zunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& zfcnName, 
                          const G4String& xbinSchemeName, 
                          const G4String& ybinSchemeName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "P2", name);
#endif
  tools::histo::p2d* p2d
    = CreateToolsP2(title, 
                    nxbins, xmin, xmax, nybins, ymin, ymax, zmin, zmax,
                    xunitName, yunitName, zunitName, 
                    xfcnName, yfcnName, zfcnName, 
                    xbinSchemeName, ybinSchemeName);
    
  // Add annotation
  AddP2Annotation(p2d, xunitName, yunitName, zunitName,
                  xfcnName, yfcnName, zfcnName);        
    
  // Save P2 information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  AddP2Information(
    name, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    xbinScheme, ybinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP2(p2d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "P2", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4int G4P2ToolsManager::CreateP2(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& xedges,
                          const std::vector<G4double>& yedges,
                          G4double zmin, G4double zmax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& zunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& zfcnName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "P2", name);
#endif
  tools::histo::p2d* p2d
    = CreateToolsP2(title, xedges, yedges, zmin, zmax, 
        xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName); 
    
  // Add annotation
  AddP2Annotation(
    p2d, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName);        
    
  // Save P2 information
  AddP2Information(
    name, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    kUserBinScheme, kUserBinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP2(p2d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "P2", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2(G4int id,
                            G4int nxbins, G4double xmin, G4double xmax, 
                            G4int nybins, G4double ymin, G4double ymax,
                            G4double zmin, G4double zmax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& zunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& zfcnName,
                            const G4String& xbinSchemeName, 
                            const G4String& ybinSchemeName)
{                                
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2", false, false);
  if ( ! p2d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetP2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "P2", info->GetName());
#endif

  // Configure tools p2
  ConfigureToolsP2(
    p2d, nxbins, xmin, xmax, nybins, ymin, ymax, zmin, zmax,
    xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    xbinSchemeName, ybinSchemeName);

  // Add annotation
  AddP2Annotation(p2d, xunitName, yunitName, zunitName, 
                  xfcnName, yfcnName, zfcnName);        
    
  // Update information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  UpdateP2Information(
    info, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName,
    xbinScheme, ybinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2(G4int id,
                            const std::vector<G4double>& xedges,
                            const std::vector<G4double>& yedges,
                            G4double zmin, G4double zmax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& zunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& zfcnName)
{                                
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2", false, false);
  if ( ! p2d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetP2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "P2", info->GetName());
#endif

  // Configure tools p2
  ConfigureToolsP2(p2d, xedges, yedges, zmin, zmax, 
    xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName);

  // Add annotation
  AddP2Annotation(p2d, xunitName, yunitName, zunitName, 
                  xfcnName, yfcnName, zfcnName);        
    
  // Update information
  UpdateP2Information(
    info, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName,
    kUserBinScheme, kUserBinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4P2ToolsManager::ScaleP2(G4int id, G4double factor)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "ScaleP2", false, false);
  if ( ! p2d ) return false;

  return p2d->scale(factor);
}  
                           
//_____________________________________________________________________________
G4bool G4P2ToolsManager::FillP2(G4int id, 
                             G4double xvalue, G4double yvalue, G4double zvalue,
                             G4double weight)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "FillP2", true, false);
  if ( ! p2d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillP2");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillP2");
  G4HnDimensionInformation* zInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kZ, "FillP2");

  p2d->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
            yInfo->fFcn(yvalue/yInfo->fUnit), 
            zInfo->fFcn(zvalue/zInfo->fUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    //description << " id " << id 
    //  << " xvalue " << xvalue << " yvalue " << yvalue << " zvalue " << zvalue;
    description << " id " << id 
                << " xvalue " << xvalue 
                << " xfcn(xvalue/xunit) " <<  xInfo->fFcn(xvalue/xInfo->fUnit) 
                << " yvalue " << yvalue
                << " yfcn(yvalue/yunit) " <<  yInfo->fFcn(yvalue/yInfo->fUnit)
                << " zvalue " << zvalue
                << " zfcn(zvalue/zunit) " <<  zInfo->fFcn(zvalue/zInfo->fUnit)
                << " weight " << weight;
    fState.GetVerboseL4()->Message("fill", "P2", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4int  G4P2ToolsManager::GetP2Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fP2NameIdMap.find(name);
  if ( it ==  fP2NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4P2ToolsManager::GetP2Id";
      G4ExceptionDescription description;
      description << "      " << "profile " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return kInvalidId;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int G4P2ToolsManager::GetP2Nxbins(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2NXbins");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetNbins(*p2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Xmin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Xmin");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMin(*p2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Xmax(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Xmax");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMax(*p2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2XWidth(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2XWidth", true, false);
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetWidth(*p2d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4int G4P2ToolsManager::GetP2Nybins(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2NYbins");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetNbins(*p2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Ymin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Ymin");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMin(*p2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Ymax(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Ymax");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMax(*p2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2YWidth(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2YWidth", true, false);
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetWidth(*p2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Zmin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Zmin");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMin(*p2d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4double G4P2ToolsManager::GetP2Zmax(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Zmax");
  if ( ! p2d ) return 0;
  
  return fBaseToolsManager.GetMax(*p2d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2Title(G4int id, const G4String& title)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2Title");
  if ( ! p2d ) return false;
  
  return fBaseToolsManager.SetTitle(*p2d, title);
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2XAxisTitle");
  if ( ! p2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*p2d, G4BaseToolsManager::kX, title);
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2YAxisTitle");
  if ( ! p2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*p2d, G4BaseToolsManager::kY, title);
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::SetP2ZAxisTitle(G4int id, const G4String& title)
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "SetP2ZAxisTitle");
  if ( ! p2d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*p2d, G4BaseToolsManager::kZ, title);
}  

//_____________________________________________________________________________
G4String G4P2ToolsManager::GetP2Title(G4int id) const
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2Title");
  if ( ! p2d ) return "";
  
  return fBaseToolsManager.GetTitle(*p2d);
}  

//_____________________________________________________________________________
G4String G4P2ToolsManager::GetP2XAxisTitle(G4int id) const 
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2XAxisTitle");
  if ( ! p2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*p2d, G4BaseToolsManager::kX);
} 

//_____________________________________________________________________________
G4String G4P2ToolsManager::GetP2YAxisTitle(G4int id) const 
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2YAxisTitle");
  if ( ! p2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*p2d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4String G4P2ToolsManager::GetP2ZAxisTitle(G4int id) const 
{
  tools::histo::p2d* p2d = GetP2InFunction(id, "GetP2ZAxisTitle");
  if ( ! p2d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*p2d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.
// Not yet available for P2

  return ! fHnManager->IsAscii();
} 

//
// public methods
// 

//_____________________________________________________________________________
G4int G4P2ToolsManager::AddP2(const G4String& name, tools::histo::p2d* p2d)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("add", "P2", name);
#endif

  // Add annotation
  AddP2Annotation(p2d, "none", "none", "none", "none", "none", "none");        
  // Add information
  AddP2Information(name, "none", "none", "none", "none", "none", "none", 
                   kLinearBinScheme, kLinearBinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP2(p2d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("add", "P2", name);
#endif
  return id;
}  

//_____________________________________________________________________________
void G4P2ToolsManager::AddP2Vector(
                          const std::vector<tools::histo::p2d*>& p2Vector)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("merge", "all p2", "");
#endif
  std::vector<tools::histo::p2d*>::const_iterator itw = p2Vector.begin();
  std::vector<tools::histo::p2d*>::iterator it;
  for (it = fP2Vector.begin(); it != fP2Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("merge", "all p2", "");
#endif
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::Reset()
{
// Reset profiles and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::p2d*>::iterator it;
  for (it = fP2Vector.begin(); it != fP2Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4P2ToolsManager::IsEmpty() const
{
  return ! fP2Vector.size();
}  
 
//_____________________________________________________________________________
tools::histo::p2d*  G4P2ToolsManager::GetP2(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const 
{
  return GetP2InFunction(id, "GetP2", warn, onlyIfActive);
}

