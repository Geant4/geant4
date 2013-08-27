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

#include "tools/histo/h2d"

#include <fstream>

using namespace G4Analysis;

//_____________________________________________________________________________
G4H2ToolsManager::G4H2ToolsManager(const G4AnalysisManagerState& state)
 : G4VH2Manager(state),
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
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
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
  
// 
// protected methods
//

//_____________________________________________________________________________
G4int G4H2ToolsManager::CreateH2(const G4String& name,  const G4String& title,
                               G4int nxbins, G4double xmin, G4double xmax,
                               G4int nybins, G4double ymin, G4double ymax,
                               const G4String& xunitName, const G4String& yunitName,
                               const G4String& xfcnName, const G4String& yfcnName)
                               
{
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("create", "H2", name);
#endif
  G4int index = fH2Vector.size();
  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  tools::histo::h2d* h2 
    = new tools::histo::h2d(title, 
                            nxbins, xfcn(xmin), xfcn(xmax), 
                            nybins, yfcn(ymin), yfcn(ymax));
            // h2 objects are deleted in destructor and reset when 
            // closing a file.

  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  h2->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);
             
  fH2Vector.push_back(h2);
  fHnManager->AddH2Information(name, xunitName, yunitName, xfcnName, yfcnName, 
                               xunit, yunit, xfcn, yfcn, 
                               kLinearBinScheme, kLinearBinScheme);
  
  fLockFirstId = true;
#ifdef G4VERBOSE
  if ( fState.GetVerboseL2() ) 
    fState.GetVerboseL2()->Message("create", "H2", name);
#endif
  fH2NameIdMap[name] = index + fFirstId;
  return index + fFirstId;
}                                         

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2(G4int id,
                                G4int nxbins, G4double xmin, G4double xmax, 
                                G4int nybins, G4double ymin, G4double ymax,
                                const G4String& xunitName, const G4String& yunitName,
                                const G4String& xfcnName, const G4String& yfcnName)
{                                
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2", false, false);
  if ( ! h2d ) return false;

  G4HnInformation* info = fHnManager->GetHnInformation(id, "SetH2");
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) 
    fState.GetVerboseL4()->Message("configure", "H2", info->fName);
#endif

  G4double xunit = GetUnitValue(xunitName);
  G4double yunit = GetUnitValue(yunitName);
  G4Fcn xfcn = GetFunction(xfcnName);
  G4Fcn yfcn = GetFunction(yfcnName);
  h2d->configure(nxbins, xfcn(xmin), xfcn(xmax), 
                 nybins, yfcn(ymin), yfcn(ymax));

  info->fXUnitName = xunitName;
  info->fYUnitName = yunitName;
  info->fXFcnName = xfcnName;
  info->fYFcnName = yfcnName;
  info->fXUnit = xunit;
  info->fYUnit = yunit;
  info->fXFcn = xfcn;
  info->fYFcn = yfcn;
  fHnManager->SetActivation(id, true); 
  
  G4String xaxisTitle;
  G4String yaxisTitle;
  UpdateTitle(xaxisTitle, xunitName, xfcnName);        
  UpdateTitle(yaxisTitle, yunitName, yfcnName);        
  h2d->add_annotation(tools::histo::key_axis_x_title(), xaxisTitle);
  h2d->add_annotation(tools::histo::key_axis_y_title(), yaxisTitle);

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

  G4HnInformation* info = fHnManager->GetHnInformation(id, "FillH2");
  h2d->fill(info->fXFcn(xvalue/info->fXUnit), 
            info->fYFcn(yvalue/info->fYUnit), weight);
#ifdef G4VERBOSE
  if ( fState.GetVerboseL4() ) {
    G4ExceptionDescription description;
    description << " id " << id 
                << " xvalue " << xvalue << " yvalue " << yvalue;
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
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }
    return -1;         
  }
  return it->second;
}  
                                      
//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nxbins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NXbins");
  if ( ! h2d ) return 0;
  
  return h2d->axis_x().bins();
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmin");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2Xmin");
  return info->fXFcn(h2d->axis_x().lower_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Xmax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Xmax");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2Xmax");
  return info->fXFcn(h2d->axis_x().upper_edge()*info->fXUnit);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2XWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XWidth", true, false);
  if ( ! h2d ) return 0;
  
  G4int nbins = h2d->axis_x().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h2 id = " << id << ").";
    G4Exception("G4H2ToolsManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2XWidth");
  return ( info->fXFcn(h2d->axis_x().upper_edge()) 
           - info->fXFcn(h2d->axis_x().lower_edge()))*info->fXUnit/nbins;
}  

//_____________________________________________________________________________
G4int G4H2ToolsManager::GetH2Nybins(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2NYbins");
  if ( ! h2d ) return 0;
  
  return h2d->axis_y().bins();
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymin(G4int id) const
{
// Returns xmin value with applied unit and histogram function

  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymin");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2Ymin");
  return info->fYFcn(h2d->axis_y().lower_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2Ymax(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Ymax");
  if ( ! h2d ) return 0;
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2Ymax");
  return info->fYFcn(h2d->axis_y().upper_edge()*info->fYUnit);
}  

//_____________________________________________________________________________
G4double G4H2ToolsManager::GetH2YWidth(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YWidth", true, false);
  if ( ! h2d ) return 0;
  
  G4int nbins = h2d->axis_y().bins();
  if ( ! nbins ) {
    G4ExceptionDescription description;
    description << "    nbins = 0 (for h2 id = " << id << ").";
    G4Exception("G4H2ToolsManager::GetH2Width",
                "Analysis_W014", JustWarning, description);
    return 0;
  }              
  
  G4HnInformation* info = fHnManager->GetHnInformation(id, "GetH2YWidth");
  return ( info->fYFcn(h2d->axis_y().upper_edge()) 
           - info->fYFcn(h2d->axis_y().lower_edge()))*info->fYUnit/nbins;
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2Title(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2Title");
  if ( ! h2d ) return false;
  
  return h2d->set_title(title);
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2XAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2XAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2YAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2YAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_x_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::SetH2ZAxisTitle(G4int id, const G4String& title)
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "SetH2ZAxisTitle");
  if ( ! h2d ) return false;
  
  h2d->add_annotation(tools::histo::key_axis_z_title(), title);
  return true;  
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2Title(G4int id) const
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2Title");
  if ( ! h2d ) return "";
  
  return h2d->title();
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2XAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2XAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_x_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get x_axis title for h2 id = " << id << ").";
    G4Exception("G4H2ToolsManager::GetH2XAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
} 

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2YAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2YAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_y_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get y_axis title for h2 id = " << id << ").";
    G4Exception("G4H2ToolsManager::GetH2YAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4String G4H2ToolsManager::GetH2ZAxisTitle(G4int id) const 
{
  tools::histo::h2d* h2d = GetH2InFunction(id, "GetH2ZAxisTitle");
  if ( ! h2d ) return "";
  
  G4String title;
  G4bool result = h2d->annotation(tools::histo::key_axis_z_title(), title);
  if ( ! result ) {
    G4ExceptionDescription description;
    description << "    Failed to get z_axis title for h2 id = " << id << ").";
    G4Exception("G4H2ToolsManager::GetH2ZAxisTitle",
                "Analysis_W014", JustWarning, description);
    return "";
  }
  
  return title;              
}  

//_____________________________________________________________________________
G4bool G4H2ToolsManager::WriteOnAscii(std::ofstream& /*output*/)
{
// Write selected objects on ASCII file
// According to the implementation by Michel Maire, originally in
// extended examples.
// Not yet available for H2

  return false;
} 

//
// public methods
// 

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

