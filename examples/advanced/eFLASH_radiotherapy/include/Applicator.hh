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
//

#ifndef Applicator_H
#define Applicator_H 1

#include "FlashDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "globals.hh"

class G4VPhysicalVolume;

class Applicator {
public:
  Applicator(G4VPhysicalVolume *);
  ~Applicator();

  G4double fFinalApplicatorXPositionFlash;
  G4double fHightFinalApplicatorFlash;

  void SetOuterRadius(G4double radius);
   void SetApplicatorLength(G4double length);

private:
void ConstructCollimator(G4VPhysicalVolume *);

  void FlashBeamLineVacuumSource();

  void FlashBeamLineTitaniumWindows();
  void FlashVWAlcover();
    void FlashAlCover2();
      void FlashExitBit();
        void FlashToroid();
          void OverCover();
            void OverCover2();
              void MonitorChamber();
                void Flash_connector();
                                void Bigconnector();
                                                void Bigconnector2();
                                                                void Bigconnector3();
                                                                void FlashBeamLineApplicator();


  G4double fInitial_pos;
  
  G4double fInnerRadiusFirstApplicatorFlash;
  G4VPhysicalVolume *fMotherPhys;

  G4Material* Fe;
  G4Material* PVDF;
  G4Material* FILM;
  G4Material* aluminumNist;
  G4Material* PMMA;

  G4double fOutRadiusVSFlash;
  G4double fHightVSFlash;
  G4double fXPositionVSFlash;
  G4double fHightFTFlash;
  
    G4double fOutRadiusFTFlash;
    G4double fOutRadius;
    G4double fToroid_outRadius;
    G4double fToroid_hight;
    G4double fToroid_XPosition;
    G4double fBigcover_hight;
    G4double fBigcover_XPosition;
    G4double fChamberpos;
  
  
  void SetDefaultDimensions();

  void ConstructApplicator();

  G4VisAttributes *blue;
  G4VisAttributes *gray;
  G4VisAttributes *white;
  G4VisAttributes *red;
  G4VisAttributes *yellow;
  G4VisAttributes *green;
  G4VisAttributes *darkGreen;
  G4VisAttributes *darkOrange3;
  G4VisAttributes *skyBlue;
  G4VisAttributes *magenta;

  G4double fOuterRadiusFirstApplicatorFlash;
  G4Tubs *fSolidFirstApplicatorFlash;
  G4VPhysicalVolume *fPhysiFirstApplicatorFlash;
  G4Material *fFirstApplicatorMaterialFlash;

  
  //  Titanium Window
  G4Tubs *solidFTFlash;
  G4VPhysicalVolume *physiFTFlash;
  G4Material *FTFlashMaterialFlash;

  //  Vacuum Source
  G4Tubs *solidVSFlash;
  G4VPhysicalVolume *physiVSFlash;
  G4Material *VSFlashMaterialFlash;
};
#endif
