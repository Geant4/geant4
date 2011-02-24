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

#include "Test2SDConstruction.hh"

#include "G4LogicalVolume.hh"

Test2SDConstruction::Test2SDConstruction(const G4String& name,G4int Segment[3])
  :fname(name)
{
  nSegment[0] = Segment[0];
  nSegment[1] = Segment[1];
  nSegment[2] = Segment[2];
}

Test2SDConstruction::~Test2SDConstruction()
{;}


#include "G4SDManager.hh"
#include "Test2PhantomSD.hh"
void Test2SDConstruction::SetupSensitivity(G4LogicalVolume* logVol) {
  //
  // sensitive detectors
  //
  G4SDManager * sdManager = G4SDManager::GetSDMpointer();
  sdManager->SetVerboseLevel(1);

  G4String phantomSDName = fname;
  Test2PhantomSD * phantomSD = new Test2PhantomSD(phantomSDName,nSegment);
  sdManager->AddNewDetector(phantomSD);
  logVol->SetSensitiveDetector(phantomSD);

  sdManager->SetVerboseLevel(0);

}
