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
/// \file visualization/standalone/src/StandaloneVisAction.cc
/// \brief Implementation of the StandaloneVisAction class
//
//

#include "ICRP110StandaloneVisAction.hh"

#include "ICRP110PhantomConstruction.hh"
#include "ICRP110PhantomMessenger.hh"
#include "G4VVisManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ICRP110StandaloneVisAction::ICRP110StandaloneVisAction()
{
  fICRP110PhantomConstruction = new ICRP110PhantomConstruction;
  new ICRP110PhantomMessenger(fICRP110PhantomConstruction);
  // Don't instantiate the detector here - give time for /phantom/ commands.
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ICRP110StandaloneVisAction::Draw() {
  // Instantiate the detector "just in time" *after* possible /phantom/ commands.
  // Instantiatie as "static" so that it is instantiated only once.
  static auto theWorld = fICRP110PhantomConstruction->Construct();
  G4VVisManager* pVisManager = G4VVisManager::GetConcreteInstance();
  if (pVisManager) {
    pVisManager->DrawGeometry(theWorld);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
