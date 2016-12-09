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
#include "G4AnalysisManagerState.hh"
#include "G4BaseHistoUtilities.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h3d"

#include <fstream>

using namespace G4Analysis;

// static definitions
const G4int G4H3ToolsManager::kDimension = 3;

//_____________________________________________________________________________
G4H3ToolsManager::G4H3ToolsManager(const G4AnalysisManagerState& state)
 : G4VH3Manager(),
   G4THnManager<tools::histo::h3d>(state, "H3")
{}

//_____________________________________________________________________________
G4H3ToolsManager::~G4H3ToolsManager()
{}

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
  hnInformation->SetDimension(kX, xunitName, xfcnName, xbinScheme);
  hnInformation->SetDimension(kY, yunitName, yfcnName, ybinScheme);
  hnInformation->SetDimension(kZ, zunitName, zfcnName, zbinScheme);
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
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto zunit = GetUnitValue(zunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto zfcn = GetFunction(zfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  auto zbinScheme = GetBinScheme(zbinSchemeName);
  
  if ( xbinScheme != G4BinScheme::kLog && ybinScheme !=  G4BinScheme::kLog && zbinScheme != G4BinScheme::kLog) {
    if ( xbinScheme == G4BinScheme::kUser || ybinScheme == G4BinScheme::kUser || zbinScheme == G4BinScheme::kUser) {
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
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto zunit = GetUnitValue(zunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto zfcn = GetFunction(zfcnName);

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
  auto xunit = GetUnitValue(xunitName);
  auto yunit = GetUnitValue(yunitName);
  auto zunit = GetUnitValue(zunitName);
  auto xfcn = GetFunction(xfcnName);
  auto yfcn = GetFunction(yfcnName);
  auto zfcn = GetFunction(zfcnName);
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  auto zbinScheme = GetBinScheme(zbinSchemeName);
  
  if ( xbinScheme != G4BinScheme::kLog && ybinScheme !=  G4BinScheme::kLog && zbinScheme !=  G4BinScheme::kLog) {
    if ( xbinScheme == G4BinScheme::kUser || ybinScheme == G4BinScheme::kUser || zbinScheme == G4BinScheme::kUser) {
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
  auto xunit = GetUnitValue(xunitName);
  auto xfcn = GetFunction(xfcnName);
  std::vector<G4double> xnewEdges;
  ComputeEdges(xedges, xunit, xfcn, xnewEdges);

  auto yunit = GetUnitValue(yunitName);
  auto yfcn = GetFunction(yfcnName);
  std::vector<G4double> ynewEdges;
  ComputeEdges(yedges, yunit, yfcn, ynewEdges);

  auto zunit = GetUnitValue(zunitName);
  auto zfcn = GetFunction(zfcnName);
  std::vector<G4double> znewEdges;
  ComputeEdges(zedges, zunit, zfcn, znewEdges);

  h3d->configure(xnewEdges, ynewEdges, znewEdges);
}

}


//
// private methods
//

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
  auto hnInformation = fHnManager->AddHnInformation(name, 3);
  hnInformation->AddDimension(xunitName, xfcnName, xbinScheme);
  hnInformation->AddDimension(yunitName, yfcnName, ybinScheme);
  hnInformation->AddDimension(zunitName, zfcnName, zbinScheme);
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
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  auto zbinScheme = GetBinScheme(zbinSchemeName);
  AddH3Information(
    name, xunitName, yunitName, zunitName, xfcnName, yfcnName, zfcnName, 
    xbinScheme, ybinScheme, zbinScheme);
    
  // Register histogram 
  G4int id = RegisterT(h3d, name); 

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
    G4BinScheme::kUser, G4BinScheme::kUser, G4BinScheme::kUser);
    
  // Register histogram 
  G4int id = RegisterT(h3d, name); 

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
  auto h3d = GetTInFunction(id, "SetH3", false, false);
  if ( ! h3d ) return false;

  auto info = fHnManager->GetHnInformation(id, "SetH3");
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
  auto xbinScheme = GetBinScheme(xbinSchemeName);
  auto ybinScheme = GetBinScheme(ybinSchemeName);
  auto zbinScheme = GetBinScheme(zbinSchemeName);
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
  auto h3d = GetTInFunction(id, "SetH3", false, false);
  if ( ! h3d ) return false;

  auto info = fHnManager->GetHnInformation(id, "SetH3");
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
    G4BinScheme::kUser, G4BinScheme::kUser, G4BinScheme::kUser);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                                  
//_____________________________________________________________________________
G4bool G4H3ToolsManager::ScaleH3(G4int id, G4double factor)
{
  auto h3d = GetTInFunction(id, "ScaleH3", false, false);
  if ( ! h3d ) return false;

  return h3d->scale(factor);
}  
                           
//_____________________________________________________________________________
G4bool G4H3ToolsManager::FillH3(G4int id, 
                             G4double xvalue, G4double yvalue, G4double zvalue,
                             G4double weight)
{
  auto h3d = GetTInFunction(id, "FillH3", true, false);
  if ( ! h3d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, kX, "FillH3");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, kY, "FillH3");
  G4HnDimensionInformation* zInfo 
    = fHnManager->GetHnDimensionInformation(id, kZ, "FillH3");

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
  return GetTId(name, warn);
}  
                                      
//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nxbins(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3NXbins");
  if ( ! h3d ) return 0;
  
  return GetNbins(*h3d, kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h3d = GetTInFunction(id, "GetH3Xmin");
  if ( ! h3d ) return 0.;
  
  return GetMin(*h3d, kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Xmax(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3Xmax");
  if ( ! h3d ) return 0.;
  
  return GetMax(*h3d, kX);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3XWidth(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3XWidth", true, false);
  if ( ! h3d ) return 0.;
  
  return GetWidth(*h3d, kX, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nybins(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3NYbins");
  if ( ! h3d ) return 0;
  
  return GetNbins(*h3d, kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h3d = GetTInFunction(id, "GetH3Ymin");
  if ( ! h3d ) return 0.;
  
  return GetMin(*h3d, kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Ymax(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3Ymax");
  if ( ! h3d ) return 0.;
  
  return GetMax(*h3d, kY);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3YWidth(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3YWidth", true, false);
  if ( ! h3d ) return 0.;
  
  return GetWidth(*h3d, kY, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4int G4H3ToolsManager::GetH3Nzbins(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3NZbins");
  if ( ! h3d ) return 0;
  
  return GetNbins(*h3d, kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Zmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h3d = GetTInFunction(id, "GetH3Zmin");
  if ( ! h3d ) return 0.;
  
  return GetMin(*h3d, kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3Zmax(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3Zmax");
  if ( ! h3d ) return 0.;
  
  return GetMax(*h3d, kZ);
}  

//_____________________________________________________________________________
G4double G4H3ToolsManager::GetH3ZWidth(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3ZWidth", true, false);
  if ( ! h3d ) return 0.;
  
  return GetWidth(*h3d, kZ, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3Title(G4int id, const G4String& title)
{
  auto h3d = GetTInFunction(id, "SetH3Title");
  if ( ! h3d ) return false;
  
  return SetTitle(*h3d, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3XAxisTitle(G4int id, const G4String& title)
{
  auto h3d = GetTInFunction(id, "SetH3XAxisTitle");
  if ( ! h3d ) return false;
  
  return SetAxisTitle(*h3d, kX, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3YAxisTitle(G4int id, const G4String& title)
{
  auto h3d = GetTInFunction(id, "SetH3YAxisTitle");
  if ( ! h3d ) return false;
  
  return SetAxisTitle(*h3d, kY, title);
}  

//_____________________________________________________________________________
G4bool G4H3ToolsManager::SetH3ZAxisTitle(G4int id, const G4String& title)
{
  auto h3d = GetTInFunction(id, "SetH3ZAxisTitle");
  if ( ! h3d ) return false;
  
  return SetAxisTitle(*h3d, kZ, title);
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3Title(G4int id) const
{
  auto h3d = GetTInFunction(id, "GetH3Title");
  if ( ! h3d ) return "";
  
  return GetTitle(*h3d);
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3XAxisTitle(G4int id) const 
{
  auto h3d = GetTInFunction(id, "GetH3XAxisTitle");
  if ( ! h3d ) return "";
  
  return GetAxisTitle(*h3d, kX, fHnManager->GetHnType());
} 

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3YAxisTitle(G4int id) const 
{
  auto h3d = GetTInFunction(id, "GetH3YAxisTitle");
  if ( ! h3d ) return "";
  
  return GetAxisTitle(*h3d, kY, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4String G4H3ToolsManager::GetH3ZAxisTitle(G4int id) const 
{
  auto h3d = GetTInFunction(id, "GetH3ZAxisTitle");
  if ( ! h3d ) return "";
  
  return GetAxisTitle(*h3d, kZ, fHnManager->GetHnType());
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
                   G4BinScheme::kLinear, G4BinScheme::kLinear, G4BinScheme::kLinear);
    
  // Register histogram 
  G4int id = RegisterT(h3d, name); 
  
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
  AddTVector(h3Vector);
}  

//_____________________________________________________________________________
tools::histo::h3d*  G4H3ToolsManager::GetH3(G4int id, G4bool warn,
                                                 G4bool onlyIfActive) const 
{
  return GetTInFunction(id, "GetH3", warn, onlyIfActive);
}

