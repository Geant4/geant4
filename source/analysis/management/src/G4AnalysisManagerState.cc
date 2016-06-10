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
// $Id: G4AnalysisManagerState.cc 66310 2012-12-17 11:56:35Z ihrivnac $

// Author: Ivana Hrivnacova, 09/07/2013  (ivana@ipno.in2p3.fr)

#include "G4AnalysisManagerState.hh"
#include "G4UnitsTable.hh"

#include <iostream>

//_____________________________________________________________________________
G4AnalysisManagerState::G4AnalysisManagerState(
                                         const G4String& type, G4bool isMaster)
 : fType(type),
   fIsMaster(isMaster),
   fIsActivation(false),
   fVerboseLevel(0),
   fCompressionLevel(1),
   fVerboseL1(type,1),
   fVerboseL2(type,2),
   fVerboseL3(type,3),
   fVerboseL4(type,4),
   fpVerboseL1(0),
   fpVerboseL2(0),
   fpVerboseL3(0),
   fpVerboseL4(0)
{
}

//
// private methods
//

//_____________________________________________________________________________
void G4AnalysisManagerState::SetVerboseLevel(G4int verboseLevel) 
{
  if ( verboseLevel == fVerboseLevel || verboseLevel < 0 ) return;
  
  fVerboseLevel = verboseLevel;
  
  if ( verboseLevel == 0 ) {
    fpVerboseL1 = 0;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 1 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = 0;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 2 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = 0;
    fpVerboseL4 = 0;
  }
  else if ( verboseLevel == 3 ) {  
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = &fVerboseL3;
    fpVerboseL4 = 0;
  }
  else {
    fpVerboseL1 = &fVerboseL1;
    fpVerboseL2 = &fVerboseL2;
    fpVerboseL3 = &fVerboseL3;
    fpVerboseL4 = &fVerboseL4;
  }
}  

