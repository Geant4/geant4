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

#ifndef TstDrawVox01DetectorConstruction_h
#define TstDrawVox01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;
class TstDrawVox01DetectorMessenger;

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class TstDrawVox01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    TstDrawVox01DetectorConstruction();
    ~TstDrawVox01DetectorConstruction();

  public:
     G4VPhysicalVolume* Construct();
     void SwitchDetector();
     void SelectDetector(G4String val);
     void SelectMaterial(G4String val);

  private:
     G4VPhysicalVolume* SelectDetector();
     void ConstructDetectors();
     void SelectMaterialPointer();

     G4LogicalVolume*   simpleBoxLog;
     G4VPhysicalVolume* simpleBoxDetector;
     G4VPhysicalVolume* honeycombDetector;
     G4Material* Air;
     G4Material* Al;
     G4Material* Pb;
     G4Material* selectedMaterial;
     G4int detectorChoice;
     G4String materialChoice;
     TstDrawVox01DetectorMessenger * detectorMessenger;
};

#endif

