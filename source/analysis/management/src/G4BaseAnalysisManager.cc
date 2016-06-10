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

#include <iostream>

//_____________________________________________________________________________
G4BaseAnalysisManager::G4BaseAnalysisManager(
                         const G4AnalysisManagerState& state)
  : fState(state),
    fFirstId(0),
    fLockFirstId(false)
{}

//_____________________________________________________________________________
G4BaseAnalysisManager::~G4BaseAnalysisManager()
{}

// 
// public methods
//

//_____________________________________________________________________________
G4bool G4BaseAnalysisManager::SetFirstId(G4int firstId) 
{
  if ( fLockFirstId ) {
    G4ExceptionDescription description;
    description 
      << "Cannot set FirstId as its value was already used.";
    G4Exception("G4BaseAnalysisManager::SetFirstId()",
                "Analysis_W013", JustWarning, description);
    return false;
  }              

  fFirstId = firstId;
  return true;
}  
