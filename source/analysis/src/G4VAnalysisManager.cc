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

// Author: Ivana Hrivnacova, 15/06/2011  (ivana@ipno.in2p3.fr)

#include "G4VAnalysisManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4VAnalysisManager::G4VAnalysisManager(const G4String& type)
  : fVerboseLevel(0),
    fFirstHistoId(0),
    fFirstNtupleId(0),
    //fHistoDirectoryName("histograms"), 
    fHistoDirectoryName("histo"), 
    fNtupleDirectoryName("ntuple"),
    fVerboseL1(type,1),
    fVerboseL2(type,2),
    fpVerboseL1(0),
    fpVerboseL2(0)
{}

//_____________________________________________________________________________
G4VAnalysisManager::~G4VAnalysisManager()
{}

//_____________________________________________________________________________
void G4VAnalysisManager::SetVerboseLevel(G4int verboseLevel) 
{
  if ( verboseLevel == fVerboseLevel || verboseLevel < 0 ) return;
  
  fVerboseLevel = verboseLevel;
  
  if ( verboseLevel == 0 ) {
    fpVerboseL1 = 0;
    fpVerboseL2 = 0;
  }
  else if ( verboseLevel == 1 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = 0;
  }
  else {
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
  }

  G4cout << "fpVerboseL1: " << fpVerboseL1 << G4endl;
  G4cout << "fpVerboseL2: " << fpVerboseL2 << G4endl;

}  


