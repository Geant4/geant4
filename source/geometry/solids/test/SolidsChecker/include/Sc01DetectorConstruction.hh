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
// $Id: Sc01DetectorConstruction.hh,v 1.2 2004-11-10 07:43:14 grichine Exp $
// ------------------------------------------------------------
//  GEANT 4 class header file 
//
//      This class is a class derived from G4VUserDetectorConstruction
//      for constructing all particles and processes.
//
//  History
//        first version              09 Sept. 1998 by S.Magni
//        modified for geometry test  11.02.04 V. Grichine 
// ------------------------------------------------------------

#ifndef Sc01DetectorConstruction_h
#define Sc01DetectorConstruction_h 1

#include "Sc01DetectorMessenger.hh"
#include "G4VSolid.hh"
#include "G4Material.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalVolume.hh"

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class Sc01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    Sc01DetectorConstruction();
    ~Sc01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     G4VPhysicalVolume* SelectDetector (const G4String& val);
     void               SwitchDetector();
     void               SetMaterial();
     G4double           GetHallSize(){return fHallSize;};

     G4LogicalVolume* GetConePolycone();

  private:

     Sc01DetectorMessenger* detectorMessenger;
     G4VSolid* aVolume;
     G4VPhysicalVolume* PhysicalVolume;
     G4Material* Water;
     G4Material* Water1;
     G4OpticalSurface* aSurface;
     G4double fHallSize;
};

#endif
