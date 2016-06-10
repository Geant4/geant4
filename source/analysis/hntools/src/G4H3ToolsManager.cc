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

#include "G4H3ToolsManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h3d"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4H3ToolsManager::G4H3ToolsManager(const G4AnalysisManagerState& state)
 : G4VH3Manager(state),
   fBaseToolsManager("H3"),
   fH3Vector(), 
   fH3NameIdMap()
{
}

//_____________________________________________________________________________
G4H3ToolsManager::~G4H3ToolsManager()
{  
  std::vector<tools::histo::h3d*>::iterator it;
  for (it = fH3Vector.begin(); it != fH3Vector.end(); it++ ) {
    delete (*it);
  }  
}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void UpdateH3Information(G4HnInformation* hnInformation,
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme,
                          G4BinScheme zbinScheme)
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
  zInformation->fBinScheme = zbinScheme;
}  
                           
//_____________________________________________________________________________
void AddH3Annotation(tools::histo::h3d* h3d,
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
  h3d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h3d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
  h3d->add_annotation(tools::histo::key_axis_z_title(), zaxisTitle);
}               
                          
//_____________________________________________________________________________
tools::histo::h3d* CreateToolsH3(
                         const G4String& title,
                         G4int nxbins, G4double xmin, G4double xmax,
                         G4int nybins, G4double ymin, G4double ymax,
                         G4int nzbins, G4double zmin, G4double zmax,
                         const G4String& xunitName,
                         const G4String& yunitName,
                         const G4String& zunitName,
                         const G4String& xfcnName, 
                         const G4String& zfcnName,
                         const G4String& yfcnName,
                         const G4String& xbinSchemeName,
                         const G4String& ybinSchemeName,
                         const G4String& zbinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  G4BinScheme zbinScheme = GetBinScheme(zbinSchemeName);
  
  if ( xbinScheme != kLogBinScheme && ybinScheme !=  kLogBinScheme && zbinScheme != kLogBinScheme) {
    if ( xbinScheme == kUserBinScheme || ybinScheme == kUserBinScheme || zbinScheme == kUserBinScheme) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H3ToolsManager::CreateH3",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::h3d(title, 
                                 nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                                 nybins, yfcn(ymin/yunit), yfcn(ymax/yunit),
                                 nzbins, zfcn(zmin/zunit), zfcn(zmax/zunit));
               // h3 objects are deleted in destructor and reset when 
               // closing a file.
  }
  else {
    // Compute edges
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    std::vector<G4double> zedges;
    ComputeEdges(nzbins, zmin, zmax, zunit, zfcn, zbinScheme, zedges);
    return new tools::histo::h3d(title, xedges, yedges, zedges); 
  }
}     

//_____________________________________________________________________________
tools::histo::h3d* CreateToolsH3(
                         const G4String& title,
                         const std::vector<G4double>& xedges,
                         const std::vector<G4double>& yedges,
                         const std::vector<G4double>& zedges,
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
  std::vector<G4double> znewEdges;
  ComputeEdges(zedges, zunit, zfcn, znewEdges);
  
  return new tools::histo::h3d(title, xnewEdges, ynewEdges, znewEdges); 
             // h3 objects are deleted in destructor and reset when 
             // closing a file.
}  

//_____________________________________________________________________________
void  ConfigureToolsH3(tools::histo::h3d* h3d,
                       G4int nxbins, G4double xmin, G4double xmax,
                       G4int nybins, G4double ymin, G4double ymax,
                       G4int nzbins, G4double zmin, G4double zmax,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& zunitName,
                       const G4String& xfcnName, 
                       const G4String& yfcnName,
                       const G4String& zfcnName,
                       const G4String& xbinSchemeName,
                       const G4String& ybinSchemeName,
                       const G4String& zbinSchemeName)
{
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4double zunit = GetUnitValue(zunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  G4Fcn zfcn = GetFunction(zfcnName);
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  G4BinScheme zbinScheme = GetBinScheme(zbinSchemeName);
  
  if ( xbinScheme != kLogBinScheme && ybinScheme !=  kLogBinScheme && zbinScheme !=  kLogBinScheme) {
    if ( xbinScheme == kUserBinScheme || ybinScheme == kUserBinScheme || zbinScheme == kUserBinScheme) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H3ToolsManager::CreateH3",
                "Analysis_W013", JustWarning, description);
    }              
    h3d->configure(nxbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                   nybins, yfcn(ymin/yunit), yfcn(ymax/yunit),
                   nzbins, zfcn(zmin/zunit), zfcn(zmax/zunit));
  }
  else {
    // Compute bins
    std::vector<G4double> xedges;
    ComputeEdges(nxbins, xmin, xmax, xunit, xfcn, xbinScheme, xedges);
    std::vector<G4double> yedges;
    ComputeEdges(nybins, ymin, ymax, yunit, yfcn, ybinScheme, yedges);
    std::vector<G4double> zedges;
    ComputeEdges(nzbins, zmin, zmax, zunit, zfcn, zbinScheme, zedges);
    h3d->configure(xedges, yedges, zedges);
  }
}     

//_____________________________________________________________________________
void  ConfigureToolsH3(tools::histo::h3d* h3d,
                       const std::vector<G4double>& xedges,
                       const std::vector<G4double>& yedges,
                       const std::vector<G4double>& zedges,
                       const G4String& xunitName,
                       const G4String& yunitName,
                       const G4String& zunitName,
                       const G4String& xfcnName,
                       const G4String& yfcnName,
                       const G4String& zfcnName)
{
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
  std::vector<G4double> znewEdges;
  ComputeEdges(zedges, zunit, zfcn, znewEdges);

  h3d->configure(xnewEdges, ynewEdges, znewEdges);
}

}


//
// private methods
//

//_____________________________________________________________________________
tools::histo::h3d*  G4H3ToolsManager::GetH3InFunction(G4int id, 
                                      G4String functionName, G4bool warn,
                                      G4bool onlyIfActive) const
{                                      
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH3Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4H3ToolsManager::";
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
  
  return fH3Vector[index];
}
  

//_____________________________________________________________________________
void G4H3ToolsManager::AddH3Information(const G4String& name,  
                          const G4String& xunitName, 
                          const G4String& yunitName, 
                          const G4String& zunitName, 
                          const G4String& xfcnName,
                          const G4String& yfcnName,
                          const G4String& zfcnName,
                          G4BinScheme xbinScheme,
                          G4BinScheme ybinScheme,
                          G4BinScheme zbinScheme) const
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
                       xbinScheme, ybinScheme, zbinScheme);
}  

//_____________________________________________________________________________
G4int G4H3ToolsManager::RegisterToolsH3(tools::histo::h3d* h3d, 
                          const G4String& name)
{
  G4int index = fH3Vector.size();
  fH3Vector.push_back(h3d);
  
  fLockFirstId = true;
  fH3NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H3ToolsManager::CreateH3(const G4String& name,  const G4String& title,
                          G4int nxbins, G4double xmin, G4double xmax,
                          G4int nybins, G4double ymin, G4double ymax,
                          G4int nzbins, G4double zmin, G4double zmax,
                          const G4String& xunitName, const G4String& yunitName, 
                          const G4String& zunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& zfcnName, 
                          const G4String& xbinSchemeName, 
                          const G4String& ybinSchemeName,
                          const G4String& zbinSchemeName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H3", name);
#endif
  tools::histo::h3d* h3d
    = CreateToolsH3(title, 
                    nxbins, xmin, xmax, nybins, ymin, ymax, nzbins, zmin, zmax,
                    xunitName, yunitName, zunitName, 
                    xfcnName, yfcnName, zfcnName, 
                    xbinSchemeName, ybinSchemeName, zbinSchemeName);
    
  // Add annotation
  AddH3Annotation(h3d, xunitName, yunitName, zunitName,
                  xfcnName, yfcnName, zfcnName);        
    
  // Save H3 information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  G4BinScheme zbinScheme = GetBinScheme(zbinSchemeName);
  AddH3Information(
    name, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    xbinScheme, ybinScheme, zbinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH3(h3d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H3", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4int G4H3ToolsManager::CreateH3(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& xedges,
                          const std::vector<G4double>& yedges,
                          const std::vector<G4double>& zedges,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& zunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& zfcnName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H3", name);
#endif
  tools::histo::h3d* h3d
    = CreateToolsH3(title, xedges, yedges, zedges, 
        xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName); 
    
  // Add annotation
  AddH3Annotation(h3d, xunitName, yunitName, zunitName, 
                  xfcnName, yfcnName, zfcnName);        
    
  // Save H3 information
  AddH3Information(
    name, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    kUserBinScheme, kUserBinScheme, kUserBinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH3(h3d, name); 

#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H3", name);
#endif

  return id;
}                                         

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3(G4int id,
                            G4int nxbins, G4double xmin, G4double xmax, 
                            G4int nybins, G4double ymin, G4double ymax,
                            G4int nzbins, G4double zmin, G4double zmax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& zunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& zfcnName,
                            const G4String& xbinSchemeName, 
                            const G4String& ybinSchemeName,
                            const G4String& zbinSchemeName)
{                                
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3", false, false);
  if ( ! h3d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetH3");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H3", info->GetName());
#endif

  // Configure tools h3
  ConfigureToolsH3(
    h3d, nxbins, xmin, xmax, nybins, ymin, ymax, nzbins, zmin, zmax,
    xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    xbinSchemeName, ybinSchemeName, zbinSchemeName);

  // Add annotation
  AddH3Annotation(h3d, xunitName, yunitName, zunitName, 
                  xfcnName, yfcnName, zfcnName);        
    
  // Update information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  G4BinScheme ybinScheme = GetBinScheme(ybinSchemeName);
  G4BinScheme zbinScheme = GetBinScheme(zbinSchemeName);
  UpdateH3Information(
    info, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName,
    xbinScheme, ybinScheme, zbinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3(G4int id,
                            const std::vector<G4double>& xedges,
                            const std::vector<G4double>& yedges,
                            const std::vector<G4double>& zedges,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& zunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& zfcnName)
{                                
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3", false, false);
  if ( ! h3d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetH3");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H3", info->GetName());
#endif

  // Configure tools h3
  ConfigureToolsH3(h3d, xedges, yedges, zedges, 
    xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName);

  // Add annotation
  AddH3Annotation(h3d, xunitName, yunitName, zunitName, 
                  xfcnName, yfcnName, zfcnName);        
    
  // Update information
  UpdateH3Information(
    info, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName,
    kUserBinScheme, kUserBinScheme, kUserBinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4H3ToolsManager::ScaleH3(G4int id, G4double factor)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "ScaleH3", false, false);
  if ( ! h3d ) return false;

  return h3d->scale(factor);
}  
                           
//_____________________________________________________________________________
G4bool G4H3ToolsManager::FillH3(G4int id, 
                             G4double xvalue, G4double yvalue, G4double zvalue,
                             G4double weight)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "FillH3", true, false);
  if ( ! h3d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillH3");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillH3");
  G4HnDimensionInformation* zInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kZ, "FillH3");

  h3d->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
            yInfo->fFcn(yvalue/yInfo->fUnit), 
            zInfo->fFcn(zvalue/zInfo->fUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue 
                << " xfcn(xvalue/xunit) " <<  xInfo->fFcn(xvalue/xInfo->fUnit) 
                << " yvalue " << yvalue
                << " yfcn(yvalue/yunit) " <<  yInfo->fFcn(yvalue/yInfo->fUnit)
                << " zvalue " << zvalue
                << " zfcn(zvalue/zunit) " <<  zInfo->fFcn(zvalue/zInfo->fUnit)
                << " weight " << weight;
    fState.GetVerboseL4()->Message("fill", "H3", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4int  G4H3ToolsManager::GetH3Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH3NameIdMap.find(name);
  if ( it ==  fH3NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4H3ToolsManager::GetH3Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return kInvalidId;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nxbins(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3NXbins");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetNbins(*h3d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Xmin");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMin(*h3d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Xmax(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Xmax");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMax(*h3d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3XWidth(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3XWidth", true, false);
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetWidth(*h3d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nybins(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3NYbins");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetNbins(*h3d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Ymin");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMin(*h3d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Ymax(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Ymax");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMax(*h3d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3YWidth(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3YWidth", true, false);
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetWidth(*h3d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nzbins(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3NZbins");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetNbins(*h3d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Zmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Zmin");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMin(*h3d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Zmax(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Zmax");
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetMax(*h3d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3ZWidth(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3ZWidth", true, false);
  if ( ! h3d ) return 0;
  
  return fBaseToolsManager.GetWidth(*h3d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3Title(G4int id, const G4String& title)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3Title");
  if ( ! h3d ) return false;
  
  return fBaseToolsManager.SetTitle(*h3d, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3XAxisTitle");
  if ( ! h3d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h3d, G4BaseToolsManager::kX, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3YAxisTitle");
  if ( ! h3d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h3d, G4BaseToolsManager::kY, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3ZAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "SetH3ZAxisTitle");
  if ( ! h3d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*h3d, G4BaseToolsManager::kZ, title);
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3Title(G4int id) const
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3Title");
  if ( ! h3d ) return "";
  
  return fBaseToolsManager.GetTitle(*h3d);
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3XAxisTitle(G4int id) const 
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3XAxisTitle");
  if ( ! h3d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h3d, G4BaseToolsManager::kX);
} 

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3YAxisTitle(G4int id) const 
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3YAxisTitle");
  if ( ! h3d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h3d, G4BaseToolsManager::kY);
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3ZAxisTitle(G4int id) const 
{
  tools::histo::h3d* h3d = GetH3InFunction(id, "GetH3ZAxisTitle");
  if ( ! h3d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*h3d, G4BaseToolsManager::kZ);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.
// Not yet available for H3

  return ! fHnManager->IsAscii();
} 

//
// public methods
// 

//_____________________________________________________________________________
G4int G4H3ToolsManager::AddH3(const G4String& name, tools::histo::h3d* h3d)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("add", "H3", name);
#endif

  // Add annotation
  AddH3Annotation(h3d, "none", "none", "none", "none", "none", "none");        
  // Add information
  AddH3Information(name, "none", "none", "none", "none", "none", "none", 
                   kLinearBinScheme, kLinearBinScheme, kLinearBinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH3(h3d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("add", "H3", name);
#endif
  return id;
}  

//_____________________________________________________________________________
void G4H3ToolsManager::AddH3Vector(
                          const std::vector<tools::histo::h3d*>& h3Vector)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("merge", "all h3", "");
#endif
  std::vector<tools::histo::h3d*>::const_iterator itw = h3Vector.begin();
  std::vector<tools::histo::h3d*>::iterator it;
  for (it = fH3Vector.begin(); it != fH3Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("merge", "all h3", "");
#endif
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::h3d*>::iterator it;
  for (it = fH3Vector.begin(); it != fH3Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  

  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::IsEmpty() const
{
  return ! fH3Vector.size();
}  
 
//_____________________________________________________________________________
tools::histo::h3d*  G4H3ToolsManager::GetH3(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const 
{
  return GetH3InFunction(id, "GetH3", warn, onlyIfActive);
}

