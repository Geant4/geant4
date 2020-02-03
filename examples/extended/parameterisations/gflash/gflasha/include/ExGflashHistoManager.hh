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
//
/// \file ExGflashHistoManager.hh
/// \brief Definition of the ExGflasHistoManager class

#ifndef ExGflashHistoManager_h
#define ExGflashHistoManager_h 1

#include "globals.hh"
#include "g4root.hh"

class ExGflashDetectorConstruction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class ExGflashHistoManager
{
  public:
    ExGflashHistoManager(ExGflashDetectorConstruction* myDet);
   ~ExGflashHistoManager();

  void InitializePerEvent();
  void FillPerEvent();

  inline void FillPerTrack(G4double,G4double);
  inline void FillPerStep (G4double,G4int,G4int);

private:
  void Book();
  G4String fFileName;
  
  ExGflashDetectorConstruction*   fDet;
  
  G4int    fVerbose;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
