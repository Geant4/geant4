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
// $Id: G4BaseAnalysisManager.cc 70604 2013-06-03 11:27:06Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4BaseAnalysisManager.hh"
#include "G4AnalysisManagerState.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4BaseAnalysisManager::G4BaseAnalysisManager(
                         const G4AnalysisManagerState& state)
  : fState(state),
    fFirstId(0),
    fLockFirstId(false)
{
}

//_____________________________________________________________________________
G4BaseAnalysisManager::~G4BaseAnalysisManager()
{
}

// 
// protected methods
//

//_____________________________________________________________________________
G4double G4BaseAnalysisManager::GetUnitValue(const G4String& unit) const
{
   G4double value = 1.;
   if ( unit != "none" ) {
     value = G4UnitDefinition::GetValueOf(unit);
     if ( value == 0. ) value = 1.; 
   }  
   return value;
}   

//_____________________________________________________________________________
G4Fcn G4BaseAnalysisManager::GetFunction(const G4String& fcnName) const
{
  G4Fcn fcn = G4FcnIdentity;
   if ( fcnName != "none" ) {
    if      ( fcnName == "log" )  fcn = std::log;
    else if ( fcnName == "log10") fcn = std::log10;
    else if ( fcnName == "exp" )  fcn = std::exp;
    else {
      G4ExceptionDescription description;
      description 
        << "    \"" << fcnName << "\" function is not supported." << G4endl
        << "    " << "No function will be applied to h1 values.";
      G4Exception("G4AnalysisMessenger::GetFunction",
                "Analysis_W013", JustWarning, description);
    }              
  }
  return fcn;            
}
    
//_____________________________________________________________________________
void G4BaseAnalysisManager::UpdateTitle(G4String& title, 
                                        const G4String& unitName, 
                                        const G4String& fcnName) const
{
  if ( fcnName != "none" )  { title += " "; title += fcnName; title += "("; }
  if ( unitName != "none" ) { title += " ["; title += unitName; title += "]";}
  if ( fcnName != "none" )  { title += ")"; }
}  
                                                          
//_____________________________________________________________________________
void  G4BaseAnalysisManager::ExceptionForHistograms(
                                         const G4String& functionName) const
{
  G4String inFunction = "G4";
  inFunction += fState.GetType();
  inFunction += "AnalysisManager::";
  inFunction += functionName;

  G4ExceptionDescription description;
  description << "      " 
              << "Histograms are not supported." ;

  G4Exception(inFunction, "Analysis_W005", JustWarning, description);
}  

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4BaseAnalysisManager::SetFirstId(G4int firstId) 
{
  if ( fLockFirstId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstHistoId as its value was already used.";
    G4Exception("G4VH1Manager::SetFirstHistoId()",
                "Analysis_W009", JustWarning, description);
    return false;
  }              

  fFirstId = firstId;
  return true;
}  
