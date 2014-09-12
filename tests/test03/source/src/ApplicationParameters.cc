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
/// \file AplicationParameters.cc
/// \brief Implementation of the AplicationParameters namespace

#include "ApplicationParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace ApplicationParameters 
{

  G4bool TestH1 = true;  
  G4bool TestH2 = true;  
  G4bool TestH3 = true;  
  G4bool TestP1 = true;  
  G4bool TestP2 = true;  
  G4bool TestNtuple = true; 
  G4bool TestRead = false; 
  G4bool TestWrite = true; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetTestAll(G4bool option) 
{
  TestH1 = option;  
  TestH2 = option;  
  TestH3 = option;  
  TestP1 = option;  
  TestP2 = option;  
  TestNtuple = option;
}   

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SetTestHn(G4bool option)
{
  TestH1 = option;  
  TestH2 = option;  
  TestH3 = option;  
}

void SetTestPn(G4bool option)
{
  TestP1 = option;  
  TestP2 = option;  
}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
