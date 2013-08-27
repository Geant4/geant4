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

//_____________________________________________________________________________
G4HnManager::G4HnManager(const G4String& hnType,
                         const G4AnalysisManagerState& state)
  : G4BaseAnalysisManager(state),
    fHnType(hnType),
    fNofActiveObjects(0),
    fNofAsciiObjects(0),
    fHnVector()
{
}

//_____________________________________________________________________________
G4HnManager::~G4HnManager()
{
  std::vector<G4HnInformation*>::iterator it;
  for (it = fHnVector.begin(); it != fHnVector.end(); it++ ) {
    delete (*it);
  }  
}

// 
// public methods
//

//_____________________________________________________________________________
void G4HnManager::AddH1Information(const G4String& name,
                                   const G4String& unitName,
                                   const G4String& fcnName,
                                   G4double unit, G4Fcn fcn,
                                   G4BinScheme binScheme)
{
  fHnVector.push_back(
    new G4HnInformation(name, unitName, unitName, fcnName, fcnName,
                        unit, unit, fcn, fcn, binScheme, binScheme));
  ++fNofActiveObjects;
}  

//_____________________________________________________________________________
void  G4HnManager::AddH2Information(const G4String& name,
                                    const G4String& xunitName, 
                                    const G4String& yunitName,
                                    const G4String& xfcnName,
                                    const G4String& yfcnName,
                                    G4double xunit, G4double yunit,
                                    G4Fcn xfcn, G4Fcn yfcn, 
                                    G4BinScheme xBinScheme, G4BinScheme yBinScheme)
{
  fHnVector.push_back(
    new G4HnInformation(name, xunitName, yunitName, xfcnName, yfcnName,
                        xunit, yunit, xfcn, yfcn, xBinScheme, yBinScheme));
  ++fNofActiveObjects;
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
      G4Exception(inFunction, "Analysis_W007", JustWarning, description);
    }  
    return 0;         
  }
  return fHnVector[index];
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
void  G4HnManager::SetActivation(G4int id, G4bool activation)
{
// Set activation to a given object

  G4HnInformation* info = GetHnInformation(id, "SetActivation");

  if ( ! info ) return;

  // Do nothing if activation does not change
  if ( info->fActivation == activation ) return;
  
  // Change activation and account it in fNofActiveObjects
  info->fActivation = activation;
  if ( activation ) 
    fNofActiveObjects++;
  else
    fNofActiveObjects--;   
}    

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4bool activation)
{
// Set activation to all objects of the given type

  std::vector<G4HnInformation*>::iterator it;
  for ( it = fHnVector.begin(); it != fHnVector.end(); it++ ) {
    G4HnInformation* info = *it;

    // Do nothing if activation does not change
    if ( info->fActivation == activation ) continue;
  
    // Change activation and account it in fNofActiveObjects
    info->fActivation = activation;
    if ( activation ) 
      fNofActiveObjects++;
    else
      fNofActiveObjects--; 
  }     
}    

//_____________________________________________________________________________
void  G4HnManager::SetAscii(G4int id, G4bool ascii)
{
  G4HnInformation* info = GetHnInformation(id, "SetAscii");

  if ( ! info ) return;

  // Do nothing if ascii does not change
  if ( info->fAscii == ascii ) return;
  
  // Change ascii and account it in fNofAsciiObjects
  info->fAscii = ascii;
  if ( ascii ) 
    fNofAsciiObjects++;
  else
    fNofAsciiObjects--;   
}    

//_____________________________________________________________________________
G4String G4HnManager::GetName(G4int id) const
{
  G4HnInformation* info = GetHnInformation(id, "GetName");

  if ( ! info ) return "";
    
  return info->fName;
}    

//_____________________________________________________________________________
G4double G4HnManager::GetXUnit(G4int id) const
{
  G4HnInformation* info = GetHnInformation(id, "GetXUnit");

  if ( ! info ) return 1.0;
  
  return info->fXUnit;
}    

//_____________________________________________________________________________
G4double G4HnManager::GetYUnit(G4int id) const
{
  G4HnInformation* info = GetHnInformation(id, "GetYUnit");

  if ( ! info ) return 1.0;
  
  return info->fYUnit;
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetActivation(G4int id) const
{
  G4HnInformation* info = GetHnInformation(id, "GetActivation");

  if ( ! info ) return true;
  
  return info->fActivation;
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetAscii(G4int id) const
{
  G4HnInformation* info = GetHnInformation(id, "GetAscii");

  if ( ! info ) return false;
  
  return info->fAscii;
}    
