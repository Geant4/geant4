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
// $Id: G4H1ToolsManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4H1ToolsManager.hh"
#include "G4HnManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h1d"

#include <fstream>

using namespace G4Analysis;

//
// Constructors, destructor
//

//_____________________________________________________________________________
G4H1ToolsManager::G4H1ToolsManager(const G4AnalysisManagerState& state)
 : G4VH1Manager(state),
   fH1Vector(),
   fH1NameIdMap()
{
}

//_____________________________________________________________________________
G4H1ToolsManager::~G4H1ToolsManager()
{  
  std::vector<tools::histo::h1d*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    delete (*it);
  }  
}

//
// Utility functions
//

namespace {

//_____________________________________________________________________________
void  UpdateH1Information(G4HnInformation* information,
                          const G4String& unitName, 
                          const G4String& fcnName,
                          G4BinScheme binScheme)
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  information->fXUnitName = unitName;
  information->fYUnitName = unitName;
  information->fXFcnName = fcnName;
  information->fYFcnName = fcnName;
  information->fXUnit = unit;
  information->fYUnit = unit;
  information->fXFcn = fcn;
  information->fYFcn = fcn;
  information->fXBinScheme = binScheme;
  information->fYBinScheme = binScheme;
}  

//_____________________________________________________________________________
void AddH1Annotation(tools::histo::h1d* h1d,
                     const G4String& unitName, 
                     const G4String& fcnName)
{                                   
  G4String axisTitle;
  UpdateTitle(axisTitle, unitName, fcnName);        
  h1d->add_annotation(tools::histo::key_axis_x_title(), axisTitle);
}  

//_____________________________________________________________________________
tools::histo::h1d* CreateToolsH1(const G4String& title,
                                 G4int nbins, G4double xmin, G4double xmax, 
                                 const G4String& fcnName,
                                 const G4String& binSchemeName)
{
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
      G4Exception("G4H1ToolsManager::CreateH1",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::h1d(title, nbins, fcn(xmin), fcn(xmax));
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, fcn, binScheme, edges);
    return new tools::histo::h1d(title, edges); 
  }
}     

//_____________________________________________________________________________
tools::histo::h1d* CreateToolsH1(const G4String& title,
                                 const std::vector<G4double>& edges,
                                 const G4String& fcnName)
{
  G4Fcn fcn = GetFunction(fcnName);

  // Apply function 
  std::vector<G4double> newEdges;
  ComputeEdges(edges, fcn, newEdges);
  
  return new tools::histo::h1d(title, newEdges); 
}  

//_____________________________________________________________________________
void ConfigureToolsH1(tools::histo::h1d* h1d,
                       G4int nbins, G4double xmin, G4double xmax,  
                       const G4String& fcnName,
                       const G4String& binSchemeName)
{
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
      G4Exception("G4H1ToolsManager::SetH1",
                "Analysis_W013", JustWarning, description);
    }              
    h1d->configure(nbins, fcn(xmin), fcn(xmax));
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, fcn, binScheme, edges);
    h1d->configure(edges);
  }
}     

//_____________________________________________________________________________
void ConfigureToolsH1(tools::histo::h1d* h1d,
                      const std::vector<G4double>& edges,
                      const G4String& fcnName)
{
  // Apply function to edges
  G4Fcn fcn = GetFunction(fcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, fcn, newEdges);

  h1d->configure(newEdges);
}

}

// 
// private methods
//

//_____________________________________________________________________________
tools::histo::h1d*  G4H1ToolsManager::GetH1InFunction(G4int id, 
                                         G4String functionName, G4bool warn,
                                         G4bool onlyIfActive) const
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fH1Vector.size()) ) {
    if ( warn) {
      G4String inFunction = "G4H1ToolsManager::";
      inFunction += functionName;
      G4ExceptionDescription description;
      description << "      " << "histogram " << id << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return 0;         
  }
  
  // Do not return histogram if inactive 
  if ( fState.GetIsActivation() && onlyIfActive && ( ! fHnManager->GetActivation(id) ) ) {
    return 0; 
  }  
  
  return fH1Vector[index];
}  

//_____________________________________________________________________________
void G4H1ToolsManager::AddH1Information(const G4String& name,  
                                        const G4String& unitName, 
                                        const G4String& fcnName,
                                        G4BinScheme binScheme) const
{
  G4double unit = GetUnitValue(unitName);
  G4Fcn fcn = GetFunction(fcnName);
  fHnManager->AddH1Information(name, unitName, fcnName, unit, fcn, binScheme);
}  
                                        
//_____________________________________________________________________________
G4int G4H1ToolsManager::RegisterToolsH1(tools::histo::h1d* h1d, 
                                        const G4String& name)
{
  G4int index = fH1Vector.size();
  fH1Vector.push_back(h1d);
  
  fLockFirstId = true;
  fH1NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H1ToolsManager::CreateH1(const G4String& name,  const G4String& title,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName,
                               const G4String& binSchemeName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H1", name);
#endif
  tools::histo::h1d* h1d
    = CreateToolsH1(title, nbins, xmin, xmax, fcnName, binSchemeName);
    
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Save H1 information
  G4BinScheme binScheme = GetBinScheme(binSchemeName);
  AddH1Information(name, unitName, fcnName, binScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH1(h1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H1", name);
#endif
  return id;
}                                         

//_____________________________________________________________________________
G4int G4H1ToolsManager::CreateH1(const G4String& name,  const G4String& title,
                                 const std::vector<G4double>& edges,
                                 const G4String& unitName, const G4String& fcnName)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H1", name);
#endif
  tools::histo::h1d* h1d 
    = CreateToolsH1(title, edges, fcnName);
    
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Save H1 information
  AddH1Information(name, unitName, fcnName, kUserBinScheme);
    
  // Register histogram 
  G4int id = RegisterToolsH1(h1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H1", name);
#endif
  return id;
}                                         

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1(G4int id,
                               G4int nbins, G4double xmin, G4double xmax,
                               const G4String& unitName, const G4String& fcnName,
                               const G4String& binSchemeName)
{                                
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1", false, false);
  if ( ! h1d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id,"SetH1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H1", info->fName);
#endif

  // Configure tools h1
  ConfigureToolsH1(h1d, nbins, xmin, xmax, fcnName, binSchemeName);

  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Update information
  G4BinScheme binScheme = GetBinScheme(binSchemeName);
  UpdateH1Information(info, unitName, fcnName, binScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1(G4int id,
                           const std::vector<G4double>& edges,
                           const G4String& unitName,
                           const G4String& fcnName)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1", false, false);
  if ( ! h1d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id,"SetH1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H1", info->fName);
#endif

  // Configure tools h1
  ConfigureToolsH1(h1d, edges, fcnName);
      
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        

  // Update information 
  UpdateH1Information(info, unitName, fcnName, kUserBinScheme);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                           
  
//_____________________________________________________________________________
G4bool G4H1ToolsManager::ScaleH1(G4int id, G4double factor)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "ScaleH1", false, false);
  if ( ! h1d ) return false;

  return h1d->scale(factor);
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::FillH1(G4int id, G4double value, G4double weight)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "FillH1", true, false);
  if ( ! h1d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillH1 for " << id << G4endl; 
    return false; 
  }  

  G4HnInformation* info = fHnManager->GetHnInformation(id, "FillId");
  h1d->fill(info->fXFcn(value/info->fXUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " id " << id << " value " << value;
    fState.GetVerboseL4()->Message("fill", "H1", description);
  }  
#endif
  return true;
}

//_____________________________________________________________________________
G4int  G4H1ToolsManager::GetH1Id(const G4String& name, G4bool warn) const
{
  std::map<G4String, G4int>::const_iterator it = fH1NameIdMap.find(name);
  if ( it ==  fH1NameIdMap.end() ) {  
    if ( warn) {
      G4String inFunction = "G4H1ToolsManager::GetH1Id";
      G4ExceptionDescription description;
      description << "      " << "histogram " << name << " does not exist.";
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  

//_____________________________________________________________________________
G4int G4H1ToolsManager::GetH1Nbins(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Nbins");
  if ( ! h1d ) return 0;
  
  return h1d->axis().bins();
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Xmin");
  if ( ! h1d ) return 0;
  
  return h1d->axis().lower_edge();
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Xmax(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Xmax");
  if ( ! h1d ) return 0;
  
  return h1d->axis().upper_edge();
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Width(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1XWidth", true, false);
  if ( ! h1d ) return 0;
  
  G4int nbins = h1d->axis().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h1 id = " << id << ").";
    G4Exception("G4H1ToolsManager::GetH1Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  return ( h1d->axis().upper_edge() - h1d->axis().lower_edge())/nbins;
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1Title(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1Title");
  if ( ! h1d ) return false;
  
  return h1d->set_title(title);
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1XAxisTitle");
  if ( ! h1d ) return false;
  
  h1d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "SetH1YAxisTitle");
  if ( ! h1d ) return false;
  
  h1d->add_annotation(tools::histo::key_axis_y_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1Title(G4int id) const
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1Title");
  if ( ! h1d ) return "";
  
  return h1d->title();
}  


//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1XAxisTitle(G4int id) const 
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1XAxisTitle");
  if ( ! h1d ) return "";
  
  G4String title;
  G4bool result = h1d->annotation(tools::histo::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h1 id = " << id << ").";
    G4Exception("G4H1ToolsManager::GetH1XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1YAxisTitle(G4int id) const 
{
  tools::histo::h1d* h1d = GetH1InFunction(id, "GetH1YAxisTitle");
  if ( ! h1d ) return "";
  
  G4String title;
  G4bool result = h1d->annotation(tools::histo::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h1 id = " << id << ").";
    G4Exception("G4H1ToolsManager::GetH1YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.

  // h1 histograms
  for ( G4int i=0; i<G4int(fH1Vector.size()); ++i ) {
    G4int id = i + fFirstId;
    G4HnInformation* info = fHnManager->GetHnInformation(id,"WriteOnAscii"); 
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->fAscii ) continue; 
    tools::histo::h1d* h1 = fH1Vector[i];

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("write on ascii", "h1d", info->fName);
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

// 
// public methods
//

//_____________________________________________________________________________
void G4H1ToolsManager::AddH1Vector(
                          const std::vector<tools::histo::h1d*>& h1Vector)
{
#ifdef G4VERBOSE
    if ( fState.GetVerboseL4() ) 
      fState.GetVerboseL4()->Message("merge", "all h1", "");
#endif
  std::vector<tools::histo::h1d*>::const_iterator itw = h1Vector.begin();
  std::vector<tools::histo::h1d*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    (*it)->add(*(*itw++));
  }  
#ifdef G4VERBOSE
    if ( fState.GetVerboseL1() ) 
      fState.GetVerboseL1()->Message("merge", "all h1", "");
#endif
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::Reset()
{
// Reset histograms and ntuple

  G4bool finalResult = true;

  std::vector<tools::histo::h1d*>::iterator it;
  for (it = fH1Vector.begin(); it != fH1Vector.end(); it++ ) {
    G4bool result = (*it)->reset();
    if ( ! result ) finalResult = false;
  }  
  
  return finalResult;
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::IsEmpty() const
{
  return ! fH1Vector.size();
}  
 
//_____________________________________________________________________________
tools::histo::h1d*  G4H1ToolsManager::GetH1(G4int id, G4bool warn,
                                            G4bool onlyIfActive) const 
{
  return GetH1InFunction(id, "GetH1", warn, onlyIfActive);
}

