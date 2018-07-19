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
// $Id: G4HnManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4HnManager.hh"
#include "G4AnalysisUtilities.hh"

using namespace G4Analysis;

//_____________________________________________________________________________
G4HnManager::G4HnManager(const G4String& hnType,
                         const G4AnalysisManagerState& state)
  : G4BaseAnalysisManager(state),
    fHnType(hnType),
    fNofActiveObjects(0),
    fNofAsciiObjects(0),
    fNofPlottingObjects(0),
    fHnVector()
{
}

//_____________________________________________________________________________
G4HnManager::~G4HnManager()
{
  for ( auto hnInformation : fHnVector ) {
    delete hnInformation;
  }
}

// 
// public methods
//

//_____________________________________________________________________________
G4HnInformation*  G4HnManager::AddHnInformation(const G4String& name, G4int nofDimensions)
{
  auto hnInformation = new G4HnInformation(name, nofDimensions);
  fHnVector.push_back(hnInformation);
  ++fNofActiveObjects;

  return hnInformation;
}  

//_____________________________________________________________________________
G4HnInformation* G4HnManager::GetHnInformation(G4int id, 
                                 G4String functionName, G4bool warn) const
{
  G4int index = id - fFirstId;
  if ( index < 0 || index >= G4int(fHnVector.size()) ) {
    if ( warn ) {
      G4String inFunction = "G4HnManager::";
      if ( functionName.size() )
        inFunction += functionName;
      else
        inFunction += "GetHnInformation"; 
      G4ExceptionDescription description;
      description << "      " << fHnType << " histogram " << id 
                  << " does not exist.";
      G4Exception(inFunction, "Analysis_W011", JustWarning, description);
    }  
    return nullptr;         
  }
  return fHnVector[index];
}    

//_____________________________________________________________________________
G4HnDimensionInformation* G4HnManager::GetHnDimensionInformation(G4int id, 
                                G4int dimension,
                                G4String functionName, G4bool warn) const
{
  auto hnInformation = GetHnInformation(id, functionName, warn);
  if ( ! hnInformation ) return nullptr; 

  return hnInformation->GetHnDimensionInformation(dimension);
}    

//_____________________________________________________________________________
G4bool G4HnManager::IsActive() const
{
  return ( fNofActiveObjects > 0 );
}  

//_____________________________________________________________________________
G4bool G4HnManager::IsAscii() const
{
  return ( fNofAsciiObjects > 0 );
}  

//_____________________________________________________________________________
G4bool G4HnManager::IsPlotting() const
{
  return ( fNofPlottingObjects > 0 );
}  

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4int id, G4bool activation)
{
// Set activation to a given object

  auto info = GetHnInformation(id, "SetActivation");

  if ( ! info ) return;

  // Do nothing if activation does not change
  if ( info->GetActivation() == activation ) return;
  
  // Change activation and account it in fNofActiveObjects
  info->SetActivation(activation);
  if ( activation ) 
    fNofActiveObjects++;
  else
    fNofActiveObjects--;   
}    

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4bool activation)
{
// Set activation to all objects of the given type

  //std::vector<G4HnInformation*>::iterator it;
  //for ( it = fHnVector.begin(); it != fHnVector.end(); it++ ) {
  //   G4HnInformation* info = *it;

  for ( auto info : fHnVector )  {

    // Do nothing if activation does not change
    if ( info->GetActivation() == activation ) continue;
  
    // Change activation and account it in fNofActiveObjects
    info->SetActivation(activation);
    if ( activation ) 
      fNofActiveObjects++;
    else
      fNofActiveObjects--; 
  }     
}    

//_____________________________________________________________________________
void  G4HnManager::SetAscii(G4int id, G4bool ascii)
{
  auto info = GetHnInformation(id, "SetAscii");

  if ( ! info ) return;

  // Do nothing if ascii does not change
  if ( info->GetAscii() == ascii ) return;
  
  // Change ascii and account it in fNofAsciiObjects
  info->SetAscii(ascii);
  if ( ascii ) 
    fNofAsciiObjects++;
  else
    fNofAsciiObjects--;   
}    

//_____________________________________________________________________________
void  G4HnManager::SetPlotting(G4int id, G4bool plotting)
{
  auto info = GetHnInformation(id, "SetPlotting");

  if ( ! info ) return;

  // Do nothing if ascii does not change
  if ( info->GetPlotting() == plotting ) return;
  
  // Change Plotting and account it in fNofPlottingObjects
  info->SetPlotting(plotting);
  if ( plotting ) 
    fNofPlottingObjects++;
  else
    fNofPlottingObjects--;   
}    

//_____________________________________________________________________________
void  G4HnManager::SetPlotting(G4bool plotting)
{
// Set plotting to all objects of the given type

  for ( auto info : fHnVector )  {

    // Do nothing if plotting does not change
    if ( info->GetPlotting() == plotting ) continue;
  
    // Change plotting and account it in fNofActiveObjects
    info->SetPlotting(plotting);
    if ( plotting ) 
      fNofPlottingObjects++;
    else
      fNofPlottingObjects--; 
  }     
}    

//_____________________________________________________________________________
G4String G4HnManager::GetName(G4int id) const
{
  auto info = GetHnInformation(id, "GetName");

  if ( ! info ) return "";
    
  return info->GetName();
}    

//_____________________________________________________________________________
G4double G4HnManager::GetXUnit(G4int id) const
{
  auto info = GetHnDimensionInformation(id, kX, "GetXUnit");

  if ( ! info ) return 1.0;
  
  return info->fUnit;
}    

//_____________________________________________________________________________
G4double G4HnManager::GetYUnit(G4int id) const
{
  auto info = GetHnDimensionInformation(id, kY, "GetYUnit");

  if ( ! info ) return 1.0;
  
  return info->fUnit;
}    

//_____________________________________________________________________________
G4double G4HnManager::GetZUnit(G4int id) const
{
  auto info = GetHnDimensionInformation(id, kZ, "GetZUnit");

  if ( ! info ) return 1.0;
  
  return info->fUnit;
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetActivation(G4int id) const
{
  auto info = GetHnInformation(id, "GetActivation");

  if ( ! info ) return true;
  
  return info->GetActivation();
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetAscii(G4int id) const
{
  auto info = GetHnInformation(id, "GetAscii");

  if ( ! info ) return false;
  
  return info->GetAscii();
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetPlotting(G4int id) const
{
  auto info = GetHnInformation(id, "GetPlotting");

  if ( ! info ) return false;
  
  return info->GetPlotting();
}    
