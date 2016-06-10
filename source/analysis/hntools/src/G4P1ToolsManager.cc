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

#include "G4P1ToolsManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/p1d"

#include <fstream>

using namespace G4Analysis;

//
// Constructors, destructor
//

//_____________________________________________________________________________
G4P1ToolsManager::G4P1ToolsManager(const G4AnalysisManagerState& state)
 : G4VP1Manager(state),
   fBaseToolsManager("P1"),
   fP1Vector(),
   fP1NameIdMap()
{
}

//_____________________________________________________________________________
G4P1ToolsManager::~G4P1ToolsManager()
{  
  std::vector<tools::histo::p1d*>::iterator it;
  for (it = fP1Vector.begin(); it != fP1Vector.end(); it++ ) {
    delete (*it);
  }  
}

//
// Utility functions
//

namespace {

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
void AddP1Annotation(tools::histo::p1d* p1d,
                     const G4String& xunitName, 
                     const G4String& yunitName, 
                     const G4String& xfcnName,
                     const G4String& yfcnName)
{                                   
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  p1d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  p1d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
}  

//_____________________________________________________________________________
tools::histo::p1d* CreateToolsP1(const G4String& title,
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
      G4Exception("G4P1ToolsManager::CreateP1",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::p1d(title, 
                                 nbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                                 yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    return new tools::histo::p1d(title, edges, yfcn(ymin/yunit), yfcn(ymax/yunit)); 
  }
}     

//_____________________________________________________________________________
tools::histo::p1d* CreateToolsP1(const G4String& title,
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
  
  return new tools::histo::p1d(title, newEdges, yfcn(ymin/yunit), yfcn(ymax/yunit)); 
}  

//_____________________________________________________________________________
void ConfigureToolsP1(tools::histo::p1d* p1d,
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
      G4Exception("G4P1ToolsManager::SetP1",
                "Analysis_W013", JustWarning, description);
    }              
    p1d->configure(nbins, xfcn(xmin/xunit), xfcn(xmax/xunit), 
                   yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, xunit, xfcn, xbinScheme, edges);
    p1d->configure(edges,  yfcn(ymin/yunit), yfcn(ymax/yunit));
  }
}     

//_____________________________________________________________________________
void ConfigureToolsP1(tools::histo::p1d* p1d,
                      const std::vector<G4double>& edges,
                      G4double ymin, G4double ymax,  
                      const G4String& xunitName,
                      const G4String& yunitName,
                      const G4String& xfcnName,
                      const G4String& yfcnName)
{
  // Apply function to edges
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, xunit, xfcn, newEdges);

  p1d->configure(newEdges, yfcn(ymin/yunit), yfcn(ymax/yunit));
}
}

// 
// private methods
//

//_____________________________________________________________________________
tools::histo::p1d*  G4P1ToolsManager::GetP1InFunction(G4int id, 
                                         G4String functionName, G4bool warn,
                                         G4bool onlyIfActive) const
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fP1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4P1ToolsManager::";
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

//_____________________________________________________________________________
void G4P1ToolsManager::AddP1Information(const G4String& name,  
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
G4int G4P1ToolsManager::RegisterToolsP1(tools::histo::p1d* p1d, 
                                        const G4String& name)
{
  G4int index = fP1Vector.size();
  fP1Vector.push_back(p1d);
  
  fLockFirstId = true;
  fP1NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4P1ToolsManager::CreateP1(const G4String& name,  const G4String& title,
                          G4int nbins, G4double xmin, G4double xmax,
                          G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName,
                          const G4String& xbinSchemeName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "P1", name);
#endif
  tools::histo::p1d* p1d
    = CreateToolsP1(title, nbins, xmin, xmax, ymin, ymax, 
                    xunitName, yunitName, xfcnName, yfcnName, 
                    xbinSchemeName);
    
  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Save P1 information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  AddP1Information(
    name, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP1(p1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "P1", name);
#endif
  return id;
}                                         

//_____________________________________________________________________________
G4int G4P1ToolsManager::CreateP1(const G4String& name,  const G4String& title,
                          const std::vector<G4double>& edges,
                          G4double ymin, G4double ymax,
                          const G4String& xunitName, const G4String& yunitName,
                          const G4String& xfcnName, const G4String& yfcnName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "P1", name);
#endif
  tools::histo::p1d* p1d 
    = CreateToolsP1(title, edges, ymin, ymax, 
                    xunitName, yunitName, xfcnName, yfcnName);
    
  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Save P1 information
  AddP1Information(
    name, xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP1(p1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "P1", name);
#endif
  return id;
}                                         

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1(G4int id,
                            G4int nbins, G4double xmin, G4double xmax,
                            G4double ymin, G4double ymax,
                            const G4String& xunitName, const G4String& yunitName,
                            const G4String& xfcnName, const G4String& yfcnName,
                            const G4String& xbinSchemeName)
{                                
  tools::histo::p1d* p1d = GetP1InFunction(id, "SetP1", false, false);
  if ( ! p1d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id,"SetP1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "P1", info->GetName());
#endif

  // Configure tools p1
  ConfigureToolsP1(
    p1d, nbins, xmin, xmax, ymin, ymax, 
    xunitName, yunitName, xfcnName, yfcnName, xbinSchemeName);

  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);        
    
  // Update information
  G4BinScheme xbinScheme = GetBinScheme(xbinSchemeName);
  UpdateP1Information(
    info, xunitName, yunitName, xfcnName, yfcnName, xbinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1(G4int id,
                           const std::vector<G4double>& edges,
                           G4double ymin, G4double ymax,
                           const G4String& xunitName, const G4String& yunitName,
                           const G4String& xfcnName, const G4String& yfcnName)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "SetP1", false, false);
  if ( ! p1d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id,"SetP1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "P1", info->GetName());
#endif

  // Configure tools p1
  ConfigureToolsP1(p1d, edges, ymin, ymax, 
                   xunitName, yunitName, xfcnName, yfcnName);
      
  // Add annotation
  AddP1Annotation(p1d, xunitName, yunitName, xfcnName, yfcnName);        

  // Update information 
  UpdateP1Information(
    info, xunitName, yunitName, xfcnName, yfcnName, kUserBinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                           
  
//_____________________________________________________________________________
G4bool G4P1ToolsManager::ScaleP1(G4int id, G4double factor)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "ScaleP1", false, false);
  if ( ! p1d ) return false;

  return p1d->scale(factor);
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::FillP1(G4int id, G4double xvalue, G4double yvalue, 
                                G4double weight)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "FillP1", true, false);
  if ( ! p1d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillP1 for " << id << G4endl; 
    return false; 
  }  

  G4HnDimensionInformation* xInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kX, "FillP1");
  G4HnDimensionInformation* yInfo 
    = fHnManager->GetHnDimensionInformation(id, G4HnInformation::kY, "FillP1");

  p1d->fill(xInfo->fFcn(xvalue/xInfo->fUnit), 
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
G4int  G4P1ToolsManager::GetP1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fP1NameIdMap.find(name);
  if ( it ==  fP1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4P1ToolsManager::GetP1Id";
      G4ExceptionDescription description;
      description << "      " << "profile " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }
    return kInvalidId;         
  }
  return it->second;
}  

//_____________________________________________________________________________
G4int G4P1ToolsManager::GetP1Nbins(G4int id) const
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Nbins");
  if ( ! p1d ) return 0;
  
  return fBaseToolsManager.GetNbins(*p1d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Xmin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Xmin");
  if ( ! p1d ) return 0;
  
  return fBaseToolsManager.GetMin(*p1d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Xmax(G4int id) const
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Xmax");
  if ( ! p1d ) return 0;
  
  return fBaseToolsManager.GetMax(*p1d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1XWidth(G4int id) const
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1XWidth", true, false);
  if ( ! p1d ) return 0;
  
  return fBaseToolsManager.GetWidth(*p1d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Ymin(G4int id) const
{
// Returns xmin value with applied unit and profile function

  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Ymin");
  if ( ! p1d ) return 0;
  
  return p1d->min_v();
}  

//_____________________________________________________________________________
G4double G4P1ToolsManager::GetP1Ymax(G4int id) const
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Ymax");
  if ( ! p1d ) return 0;
  
  return p1d->max_v();
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1Title(G4int id, const G4String& title)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "SetP1Title");
  if ( ! p1d ) return false;
  
  return fBaseToolsManager.SetTitle(*p1d, title);
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "SetP1XAxisTitle");
  if ( ! p1d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*p1d, G4BaseToolsManager::kX, title);
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::SetP1YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "SetP1YAxisTitle");
  if ( ! p1d ) return false;
  
  return fBaseToolsManager.SetAxisTitle(*p1d, G4BaseToolsManager::kY, title);
}  

//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1Title(G4int id) const
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1Title");
  if ( ! p1d ) return "";
  
  return fBaseToolsManager.GetTitle(*p1d);
}  


//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1XAxisTitle(G4int id) const 
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1XAxisTitle");
  if ( ! p1d ) return "";

  return fBaseToolsManager.GetAxisTitle(*p1d, G4BaseToolsManager::kX);
}  

//_____________________________________________________________________________
G4String G4P1ToolsManager::GetP1YAxisTitle(G4int id) const 
{
  tools::histo::p1d* p1d = GetP1InFunction(id, "GetP1YAxisTitle");
  if ( ! p1d ) return "";
  
  return fBaseToolsManager.GetAxisTitle(*p1d, G4BaseToolsManager::kY);
}  

// 
// public methods
//

//_____________________________________________________________________________
G4int G4P1ToolsManager::AddP1(const G4String& name, tools::histo::p1d* p1d)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("add", "P1", name);
#endif
    
  // Add annotation
  AddP1Annotation(p1d, "none", "none", "none", "none");        
  // Add information
  AddP1Information(name, "none", "none", "none", "none", kLinearBinScheme);
    
  // Register profile 
  G4int id = RegisterToolsP1(p1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("add", "P1", name);
#endif
  return id;
}  

//_____________________________________________________________________________
void G4P1ToolsManager::AddP1Vector(
                          const std::vector<tools::histo::p1d*>& p1Vector)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("merge", "all p1", "");
#endif
  std::vector<tools::histo::p1d*>::const_iterator itw = p1Vector.begin();
  std::vector<tools::histo::p1d*>::iterator it;
  for (it = fP1Vector.begin(); it != fP1Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("merge", "all p1", "");
#endif
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::Reset()
{
// Reset profiles and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::p1d*>::iterator it;
  for (it = fP1Vector.begin(); it != fP1Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  
  
  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4P1ToolsManager::IsEmpty() const
{
  return ! fP1Vector.size();
}  
 
//_____________________________________________________________________________
tools::histo::p1d*  G4P1ToolsManager::GetP1(G4int id, G4bool warn,
                                            G4bool onlyIfActive) const 
{
  return GetP1InFunction(id, "GetP1", warn, onlyIfActive);
}

