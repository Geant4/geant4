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
// $Id: Tst52TrackerSD.hh,v 1.1.2.1 2007-12-10 16:34:03 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//
// History:
// -----------
// 17 May     2003 SG      first implemntation
// -------------------------------------------------------------------

#ifndef Tst52TrackerSD_h
#define Tst52TrackerSD_h 1

#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "Tst52TrackerHit.hh"

class G4Step;
class G4HCofThisEvent;

class Tst52DetectorConstruction;

class Tst52TrackerSD : public G4VSensitiveDetector
{
public:
  Tst52TrackerSD(G4String,
                 Tst52DetectorConstruction*);
  ~Tst52TrackerSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);
  void ChangeVoxelNumber(G4int newNumber);
  void SetSDParameters(G4double ZSize, G4int voxelZ);

private:
  Tst52DetectorConstruction* Detector;
  G4int numberOfVoxelZ;
  G4int i;
  G4double targetZ; // Size along the Z axis 
  G4double totalEnergyDeposit; 
  Tst52TrackerHitsCollection* trackerCollection;
};
#endif

