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
// $Id: Tst14ActionInitializer.cc 66241 2012-12-13 18:34:42Z gunter $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Tst14ActionInitializer.hh"
#include "Tst14PrimaryGeneratorAction.hh"
#include "Tst14RunAction.hh"
#include "Tst14SteppingAction.hh"
#include "Tst14TrackingAction.hh"
#include "Tst14DetectorConstruction.hh"

#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Tst14ActionInitializer::Tst14ActionInitializer(Tst14DetectorConstruction* det) 
  : detector(det)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst14ActionInitializer::Build() const 
{
  /*
  const Tst14DetectorConstruction* detector = 
    static_cast<const Tst14DetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
  */
  // primary generator
  Tst14PrimaryGeneratorAction* primary = new Tst14PrimaryGeneratorAction(detector);
  SetUserAction(primary);
    
  //Thread-local run action
  //SetUserAction(new Tst14RunAction());
  //SetUserAction(new Tst14SteppingAction());
  //SetUserAction(new Tst14TrackingAction());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Tst14ActionInitializer::BuildForMaster() const
{
  SetUserAction(new Tst14RunAction());
}

