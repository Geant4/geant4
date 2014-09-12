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
// 
/// \file ApplicationParameters.hh
/// \brief Definition of the ApplicationParameters namespace

#ifndef ApplicationParameters_h
#define ApplicationParameters_h 1

#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace ApplicationParameters
{
  // Functions
  void SetTestAll(G4bool testAll);
  void SetTestHn(G4bool testAllHn);
  void SetTestPn(G4bool testAllPn);

  // Parameters
  extern G4bool TestH1;  
  extern G4bool TestH2;  
  extern G4bool TestH3;  
  extern G4bool TestP1;  
  extern G4bool TestP2;  
  extern G4bool TestNtuple; 
  extern G4bool TestRead;  
  extern G4bool TestWrite;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

