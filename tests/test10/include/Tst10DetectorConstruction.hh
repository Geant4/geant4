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
// $Id: Tst10DetectorConstruction.hh,v 1.6 2006-06-29 21:38:12 gunter Exp $
// ------------------------------------------------------------
//  GEANT 4 class header file 
//
//      This class is a class derived from G4VUserDetectorConstruction
//      for constructing all particles and processes.
//
//  History
//        first version              09 Sept. 1998 by S.Magni
// ------------------------------------------------------------

#ifndef Tst10DetectorConstruction_h
#define Tst10DetectorConstruction_h 1

#include "Tst10DetectorMessenger.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Tst10DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Tst10DetectorConstruction();
    ~Tst10DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     G4VPhysicalVolume* SelectDetector (const G4String& val);
     void               SwitchDetector();
     G4bool             CleanGeometry();
     void               SetMaterial();
     G4double           GetHallSize(){return fHallSize;}

  private:
     Tst10DetectorMessenger* detectorMessenger;
     G4VSolid* aVolume;
     G4VPhysicalVolume* PhysicalVolume;
     G4Material* Water;
     G4Material* Water1;
     G4OpticalSurface* aSurface;
     G4double fHallSize;
};

#endif
