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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4HnManager.hh"
#include "G4VFileManager.hh"
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
    fNofFileNameObjects(0),
    fHnVector(),
    fFileManager(nullptr)
{
}

//_____________________________________________________________________________
G4HnManager::~G4HnManager()
{
  for ( auto info : fHnVector ) {
    delete info;
  }
}

//
// private methods
//

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4HnInformation* info, G4bool activation)
{
// Set activation to a given object

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
void  G4HnManager::SetPlotting(G4HnInformation* info, G4bool plotting)
{
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
void  G4HnManager::SetFileName(G4HnInformation* info, const G4String& fileName)
{
  // Do nothing if file name does not change
  if ( info->GetFileName() == fileName ) return;

  // Save the info and account a new file name if file manager
  info->SetFileName(fileName);
  if (fFileManager) {
    fFileManager->AddFileName(fileName);
  } else {
    G4ExceptionDescription description;
    description
      << "Failed to set fileName " << fileName << " for object " << info->GetName() << G4endl
      << "File manager is not set.";
    G4Exception("G4HnManager::SetFileName",
                "Analysis_W012", JustWarning, description);
    return;
  }

  if ( fileName != "" ) {
    fNofFileNameObjects++;
  } else {
    fNofFileNameObjects--;
  }
}

// 
// public methods
//

//_____________________________________________________________________________
G4HnInformation*  G4HnManager::AddHnInformation(const G4String& name, G4int nofDimensions)
{
  auto info = new G4HnInformation(name, nofDimensions);
  fHnVector.push_back(info);
  ++fNofActiveObjects;

  return info;
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
  auto info = GetHnInformation(id, functionName, warn);
  if ( ! info ) return nullptr;

  return info->GetHnDimensionInformation(dimension);
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
G4bool G4HnManager::IsFileName() const
{
  return ( fNofFileNameObjects > 0 );
}

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4int id, G4bool activation)
{
// Set activation to a given object

  auto info = GetHnInformation(id, "SetActivation");

  if ( ! info ) return;

  SetActivation(info, activation);
}    

//_____________________________________________________________________________
void  G4HnManager::SetActivation(G4bool activation)
{
// Set activation to all objects of the given type

  //std::vector<G4HnInformation*>::iterator it;
  //for ( it = fHnVector.begin(); it != fHnVector.end(); it++ ) {
  //   G4HnInformation* info = *it;

  for ( auto info : fHnVector )  {
    SetActivation(info, activation);
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

  SetPlotting(info, plotting);
}

//_____________________________________________________________________________
void  G4HnManager::SetPlotting(G4bool plotting)
{
// Set plotting to all objects of the given type

  for ( auto info : fHnVector )  {
    SetPlotting(info, plotting);
  }
}

//_____________________________________________________________________________
void  G4HnManager::SetFileName(G4int id, const G4String& fileName)
{
  auto info = GetHnInformation(id, "SetFileName");

  if ( ! info ) return;

  SetFileName(info, fileName);
}    

//_____________________________________________________________________________
void  G4HnManager::SetFileName(const G4String& fileName)
{
// Set plotting to all objects of the given type

  for ( auto info : fHnVector )  {
    SetFileName(info, fileName);
  }
}    

//_____________________________________________________________________________
G4bool G4HnManager::SetXAxisIsLog(G4int id, G4bool isLog)
{
  auto info = GetHnInformation(id, "SetXAxisIsLog");

  if ( ! info ) return false;

  info->SetIsLogAxis(kX, isLog);
  return true;
}

//_____________________________________________________________________________
G4bool  G4HnManager::SetYAxisIsLog(G4int id, G4bool isLog)
{
  auto info = GetHnInformation(id, "SetYAxisIsLog");

  if ( ! info ) return false;

  info->SetIsLogAxis(kY, isLog);
  return true;
}

//_____________________________________________________________________________
G4bool  G4HnManager::SetZAxisIsLog(G4int id, G4bool isLog)
{
  auto info = GetHnInformation(id, "SetZAxisIsLog");

  if ( ! info ) return false;

  info->SetIsLogAxis(kZ, isLog);
  return true;
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
G4bool G4HnManager::GetXAxisIsLog(G4int id) const
{
  auto info = GetHnInformation(id, "GetXAxisIsLog");

  if ( ! info ) return false;
  
  return info->GetIsLogAxis(kX);
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetYAxisIsLog(G4int id) const
{
  auto info = GetHnInformation(id, "GetYAxisIsLog");

  if ( ! info ) return 1.0;
  
  return info->GetIsLogAxis(kY);
}    

//_____________________________________________________________________________
G4bool G4HnManager::GetZAxisIsLog(G4int id) const
{
  auto info = GetHnInformation(id, "GetZAxisIsLog");

  if ( ! info ) return 1.0;
  
  return info->GetIsLogAxis(kZ);
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

//_____________________________________________________________________________
G4String G4HnManager::GetFileName(G4int id) const
{
  auto info = GetHnInformation(id, "GetFileName");

  if ( ! info ) return "";
    
  return info->GetFileName();
}
