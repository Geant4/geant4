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
#include "G4AnalysisManagerState.hh"
#include "G4BaseHistoUtilities.hh"
#include "G4AnalysisUtilities.hh"

#include "tools/histo/h1d"

#include <fstream>

using namespace G4Analysis;

// static definitions
const G4int G4H1ToolsManager::kDimension = 1;

//
// Constructors, destructor
//

//_____________________________________________________________________________
G4H1ToolsManager::G4H1ToolsManager(const G4AnalysisManagerState& state)
 : G4VH1Manager(),
   G4THnManager<tools::histo::h1d>(state, "H1")
{}

//_____________________________________________________________________________
G4H1ToolsManager::~G4H1ToolsManager()
{}

//
// Utility functions
//


namespace {

//_____________________________________________________________________________
void  UpdateH1Information(G4HnInformation* hnInformation,
                          const G4String& unitName, 
                          const G4String& fcnName,
                          G4BinScheme binScheme)
{
  hnInformation->SetDimension(kX, unitName, fcnName, binScheme);
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
                                 const G4String& unitName,
                                 const G4String& fcnName,
                                 const G4String& binSchemeName)
{
  auto unit = GetUnitValue(unitName);
  auto fcn = GetFunction(fcnName);
  auto binScheme = GetBinScheme(binSchemeName);
  
  if ( binScheme != G4BinScheme::kLog ) {
    if ( binScheme == G4BinScheme::kUser ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H1ToolsManager::CreateH1",
                "Analysis_W013", JustWarning, description);
    }              
    return new tools::histo::h1d(title, nbins, fcn(xmin/unit), fcn(xmax/unit));
  }
  else {
    // Compute edges
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, unit, fcn, binScheme, edges);
    return new tools::histo::h1d(title, edges); 
  }
}     

//_____________________________________________________________________________
tools::histo::h1d* CreateToolsH1(const G4String& title,
                                 const std::vector<G4double>& edges,
                                 const G4String& unitName,
                                 const G4String& fcnName)
{
  auto unit = GetUnitValue(unitName);
  auto fcn = GetFunction(fcnName);

  // Apply function 
  std::vector<G4double> newEdges;
  ComputeEdges(edges, unit, fcn, newEdges);
  
  return new tools::histo::h1d(title, newEdges); 
}  

//_____________________________________________________________________________
void ConfigureToolsH1(tools::histo::h1d* h1d,
                       G4int nbins, G4double xmin, G4double xmax,  
                       const G4String& unitName,
                       const G4String& fcnName,
                       const G4String& binSchemeName)
{
  auto unit = GetUnitValue(unitName);
  auto fcn = GetFunction(fcnName);
  auto binScheme = GetBinScheme(binSchemeName);

  if ( binScheme != G4BinScheme::kLog ) {
    if ( binScheme == G4BinScheme::kUser ) {
      // This should never happen, but let's make sure about it
      // by issuing a warning
      G4ExceptionDescription description;
      description 
        << "    User binning scheme setting was ignored." << G4endl
        << "    Linear binning will be applied with given (nbins, xmin, xmax) values";
      G4Exception("G4H1ToolsManager::SetH1",
                "Analysis_W013", JustWarning, description);
    }              
    h1d->configure(nbins, fcn(xmin/unit), fcn(xmax/unit));
  }
  else {
    // Compute bins
    std::vector<G4double> edges;
    ComputeEdges(nbins, xmin, xmax, unit, fcn, binScheme, edges);
    h1d->configure(edges);
  }
}     

//_____________________________________________________________________________
void ConfigureToolsH1(tools::histo::h1d* h1d,
                      const std::vector<G4double>& edges,
                      const G4String& unitName,
                      const G4String& fcnName)
{
  // Apply function to edges
  auto unit = GetUnitValue(unitName);
  auto fcn = GetFunction(fcnName);
  std::vector<G4double> newEdges;
  ComputeEdges(edges, unit, fcn, newEdges);

  h1d->configure(newEdges);
}

}

// 
// private methods
//

//_____________________________________________________________________________
void G4H1ToolsManager::AddH1Information(const G4String& name,  
                                        const G4String& unitName, 
                                        const G4String& fcnName,
                                        G4BinScheme binScheme) const
{
  auto hnInformation = fHnManager->AddHnInformation(name, kDimension);
  hnInformation->AddDimension(unitName, fcnName, binScheme);
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
  
  // Create H1
  auto h1d
    = CreateToolsH1(title, nbins, xmin, xmax, unitName, fcnName, binSchemeName);
    
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Save H1 information
  auto binScheme = GetBinScheme(binSchemeName);
  AddH1Information(name, unitName, fcnName, binScheme);
    
  // Register histogram 
  auto id = RegisterT(h1d, name); 
  
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
  auto h1d 
    = CreateToolsH1(title, edges, unitName, fcnName);
    
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Save H1 information
  AddH1Information(name, unitName, fcnName, G4BinScheme::kUser);
    
  // Register histogram 
  auto id = RegisterT(h1d, name); 
  
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
  auto h1d = GetTInFunction(id, "SetH1", false, false);
  if ( ! h1d ) return false;

  auto info = fHnManager->GetHnInformation(id,"SetH1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H1", info->GetName());
#endif

  // Configure tools h1
  ConfigureToolsH1(h1d, nbins, xmin, xmax, unitName, fcnName, binSchemeName);

  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        
    
  // Update information
  auto binScheme = GetBinScheme(binSchemeName);
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
  auto h1d = GetTInFunction(id, "SetH1", false, false);
  if ( ! h1d ) return false;

  auto info = fHnManager->GetHnInformation(id,"SetH1");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H1", info->GetName());
#endif

  // Configure tools h1
  ConfigureToolsH1(h1d, edges, unitName, fcnName);
      
  // Add annotation
  AddH1Annotation(h1d, unitName, fcnName);        

  // Update information 
  UpdateH1Information(info, unitName, fcnName, G4BinScheme::kUser);

  // Set activation
  fHnManager->SetActivation(id, true); 
  
  return true;
}
                           
  
//_____________________________________________________________________________
G4bool G4H1ToolsManager::ScaleH1(G4int id, G4double factor)
{
  auto h1d = GetTInFunction(id, "ScaleH1", false, false);
  if ( ! h1d ) return false;

  return h1d->scale(factor);
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::FillH1(G4int id, G4double value, G4double weight)
{
  auto h1d = GetTInFunction(id, "FillH1", true, false);
  if ( ! h1d ) return false;

  if ( fState.GetIsActivation() && ( ! fHnManager->GetActivation(id) ) ) {
    //G4cout << "Skipping FillH1 for " << id << G4endl; 
    return false; 
  }  

  auto info 
    = fHnManager->GetHnDimensionInformation(id, kX, "FillH1");
  h1d->fill(info->fFcn(value/info->fUnit), weight);
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
G4int  G4H1ToolsManager::GetH1Id(const G4String& name, G4bool warn) const
{
  return GetTId(name, warn);
}  

//_____________________________________________________________________________
G4int G4H1ToolsManager::GetH1Nbins(G4int id) const
{
  auto h1d = GetTInFunction(id, "GetH1Nbins");
  if ( ! h1d ) return 0;
  
  return GetNbins(*h1d, kX);
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  auto h1d = GetTInFunction(id, "GetH1Xmin");
  if ( ! h1d ) return 0.;
  
  return GetMin(*h1d, kX);
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Xmax(G4int id) const
{
  auto h1d = GetTInFunction(id, "GetH1Xmax");
  if ( ! h1d ) return 0.;
  
  return GetMax(*h1d, kX);
}  

//_____________________________________________________________________________
G4double G4H1ToolsManager::GetH1Width(G4int id) const
{
  auto h1d = GetTInFunction(id, "GetH1XWidth", true, false);
  if ( ! h1d ) return 0.;
  
  return GetWidth(*h1d, kX, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1Title(G4int id, const G4String& title)
{
  auto h1d = GetTInFunction(id, "SetH1Title");
  if ( ! h1d ) return false;
  
  return SetTitle(*h1d, title);
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1XAxisTitle(G4int id, const G4String& title)
{
  auto h1d = GetTInFunction(id, "SetH1XAxisTitle");
  if ( ! h1d ) return false;
  
  return SetAxisTitle(*h1d, kX, title);
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::SetH1YAxisTitle(G4int id, const G4String& title)
{
  auto h1d = GetTInFunction(id, "SetH1YAxisTitle");
  if ( ! h1d ) return false;
  
  return SetAxisTitle(*h1d, kY, title);
}  

//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1Title(G4int id) const
{
  auto h1d = GetTInFunction(id, "GetH1Title");
  if ( ! h1d ) return "";
  
  return GetTitle(*h1d);
}  


//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1XAxisTitle(G4int id) const 
{
  auto h1d = GetTInFunction(id, "GetH1XAxisTitle");
  if ( ! h1d ) return "";
  
  return GetAxisTitle(*h1d, kX, fHnManager->GetHnType());
}

//_____________________________________________________________________________
G4String G4H1ToolsManager::GetH1YAxisTitle(G4int id) const 
{
  auto h1d = GetTInFunction(id, "GetH1YAxisTitle");
  if ( ! h1d ) return "";
  
  return GetAxisTitle(*h1d, kY, fHnManager->GetHnType());
}  

//_____________________________________________________________________________
G4bool G4H1ToolsManager::WriteOnAscii(std::ofstream& output)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.

  // h1 histograms
  for ( G4int i=0; i<G4int(fTVector.size()); ++i ) {
    auto id = i + fHnManager->GetFirstId();
    auto info = fHnManager->GetHnInformation(id,"WriteOnAscii"); 
    // skip writing if activation is enabled and H1 is inactivated
    if ( ! info->GetAscii() ) continue; 
    auto h1 = fTVector[i];

#ifdef G4VERBOSE
    if ( fState.GetVerboseL3() ) 
      fState.GetVerboseL3()->Message("write on ascii", "h1d", info->GetName());
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
G4int G4H1ToolsManager::AddH1(const G4String& name, tools::histo::h1d* h1d)
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("add", "H1", name);
#endif
    
  // Add annotation
  AddH1Annotation(h1d, "none", "none");        
  // Add information
  AddH1Information(name, "none", "none", G4BinScheme::kLinear);
    
  // Register histogram 
  auto id = RegisterT(h1d, name); 
  
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("add", "H1", name);
#endif
  return id;
}  


//_____________________________________________________________________________
void G4H1ToolsManager::AddH1Vector(
                          const std::vector<tools::histo::h1d*>& h1Vector)
{
  AddTVector(h1Vector);
}  

//_____________________________________________________________________________
tools::histo::h1d*  G4H1ToolsManager::GetH1(G4int id, G4bool warn,
                                            G4bool onlyIfActive) const 
{
  return GetTInFunction(id, "GetH1", warn, onlyIfActive);
}

