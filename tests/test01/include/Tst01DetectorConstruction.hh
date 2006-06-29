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

#ifndef Tst01DetectorConstruction_h
#define Tst01DetectorConstruction_h 1

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4AssemblyVolume;
class G4DisplacedSolid;
class G4VSolid ;
class G4Material;
class Tst01DetectorMessenger;

#include <vector>
#include "globals.hh"

#include "G4VUserDetectorConstruction.hh"

class Tst01DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

     Tst01DetectorConstruction();
    ~Tst01DetectorConstruction();

     G4VPhysicalVolume* Construct();

     void SelectDetector(const G4String& val) ;
     void SwitchDetector() ;

     void SelectMaterial(const G4String& val) ;

  // Select/Switch CSG/Boolean

     void SelectCSG(const G4String& name) ;
     void SwitchCSG() ;

     void SelectBoolean(const G4String& name) ;
     void SwitchBoolean() ;

     G4double GetWorldSize() const { return wSize; }

  private:

  // Private functions

     G4VPhysicalVolume* SelectDetector();
     void ConstructDetectors();
     void SelectMaterialPointer();

  // Class members

     Tst01DetectorMessenger* detectorMessenger ;

     G4double wSize;

     G4LogicalVolume*   simpleBoxLog ;
     G4VPhysicalVolume* simpleBoxDetector ;
     G4VPhysicalVolume* honeycombDetector ;
     G4VPhysicalVolume* fWorldPhysVol ;

  // CSGs /logic/physics volumes

     G4VSolid*          fTestCSG ;
     G4LogicalVolume*   fTestLog ;
     G4VPhysicalVolume* fTestVol ;

  // Bools /logic/physics volumes + daughters

     G4DisplacedSolid*  fDisPb ;
     G4VSolid*          fPb1, *fPb2, *fPb3, *fSphere;

     G4VSolid*          fTestBool ;
     G4LogicalVolume*   fTestBoolLog ;
     G4VPhysicalVolume* fTestBoolVol ;
     G4LogicalVolume*   fTestD1Log ;
     G4VPhysicalVolume* fTestD1Vol ;

     G4Material* Air ;
     G4Material* Al ;
     G4Material* Pb ;
     G4Material* selectedMaterial ;

     G4int       detectorChoice ;
     G4String    materialChoice ;

     G4int       fChoiceCSG ;
     G4int       fChoiceBool ;

     // ------------------------------------------
     // Assembly
     G4LogicalVolume*                AssemblyDetectorLog;
     G4VPhysicalVolume*              AssemblyDetector;
     G4AssemblyVolume*               AssemblyCalo;
     G4LogicalVolume*                AssemblyCellLog;

     // ------------------------------------------
     // Assembly2 (assembly of assemblies)
     G4VPhysicalVolume*              AssemblyDetector2;

     // ------------------------------------------
     // Assembly3 (assembly with reflections )
     G4VPhysicalVolume*              AssemblyDetector3;
};

#endif

