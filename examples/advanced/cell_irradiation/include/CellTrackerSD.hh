//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//    **************************************
//    *                                    *
//    *           CellTrackerSD.hh         *
//    *                                    *
//    **************************************
//
//
// Author: Susanna Guatelli (guatelli@ge.infn.it)
//	   Barbara Mascialino (Barbara.Mascialino@ge.infn.it)
//
// History:
// -----------
// 20 September 2006 S. Guatelli, B. Mascialino      first implemntation
// -------------------------------------------------------------------

#ifndef CellTrackerSD_h
#define CellTrackerSD_h 1

#include "globals.hh"
#include "G4VSensitiveDetector.hh"
#include "CellTrackerHit.hh"

class G4Step;
class G4HCofThisEvent;

class CellDetectorConstruction;

class CellTrackerSD : public G4VSensitiveDetector
{
public:
  CellTrackerSD(G4String,
		 G4double, G4double, G4double, // X, Y, Z Size 
		 G4int, G4int, G4int,// Number of voxels along X, Y, Z
                 CellDetectorConstruction*);
  ~CellTrackerSD();

  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step*, G4TouchableHistory*);
  void EndOfEvent(G4HCofThisEvent*);

private:
  CellTrackerHitsCollection* trackerCollection;
  G4int*                   hitID;
  CellDetectorConstruction* Detector;
  G4int numberOfVoxelX;
  G4int numberOfVoxelY; 
  G4int numberOfVoxelZ;
       
  G4double targetZ; // Size along the Z axis
  G4double targetX; 
  G4double targetY;     
  G4double totalEnergyDeposit;
  
};
#endif

