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
// $Id: exGPSActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file exGPSActionInitialization.cc
/// \brief Implementation of the exGPSActionInitialization class

#include "exGPSActionInitialization.hh"
#include "exGPSHistoManager.hh"
#include "exGPSPrimaryGeneratorAction.hh"
#include "exGPSRunAction.hh"
#include "exGPSEventAction.hh"
#include "exGPSGeometryConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSActionInitialization::exGPSActionInitialization()
 : G4VUserActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

exGPSActionInitialization::~exGPSActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSActionInitialization::BuildForMaster() const
{
  // Histo manager
  exGPSHistoManager*  histo = new exGPSHistoManager();
  
  // Actions
  SetUserAction(new exGPSRunAction(histo));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void exGPSActionInitialization::Build() const
{
  // Histo manager
  exGPSHistoManager*  histo = new exGPSHistoManager();
  
  // Actions
  //
  SetUserAction(new exGPSPrimaryGeneratorAction());
  
  exGPSRunAction* runAction = new exGPSRunAction(histo);
  SetUserAction(runAction);
  
  exGPSEventAction* eventAction = new exGPSEventAction(histo);
  SetUserAction(eventAction);


}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
