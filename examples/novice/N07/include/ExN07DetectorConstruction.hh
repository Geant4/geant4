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
//
// $Id: ExN07DetectorConstruction.hh,v 1.3 2005/11/22 22:20:55 asaim Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// 

#ifndef ExN07DetectorConstruction_h
#define ExN07DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVReplica;
class G4Material;
class G4Box;
class ExN07DetectorMessenger;

class ExN07DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN07DetectorConstruction();
    virtual ~ExN07DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
     
  public:
    void PrintCalorParameters() const;
    void SetAbsorberMaterial(G4String materialChoice);     
    G4String GetAbsorberMaterial() const;
    void SetGapMaterial(G4String materialChoice);     
    G4String GetGapMaterial() const;
    void SetSerialGeometry(G4bool ser);
    void SetNumberOfLayers(G4int nl);
    inline G4int GetNumberOfLayers() const
    { return numberOfLayers; }
    inline G4bool IsSerial() const
    { return serial; }
     
  private:
    void DefineMaterials();
    void SetupGeometry();
    void SetupDetectors();

  private:
    G4int              numberOfLayers;
    G4double           totalThickness; // total thinkness of one calorimeter
    G4double           layerThickness; // = totalThickness / numberOfLayers

    G4bool             constructed;
    
    G4String           calName[3];

    G4Material*        worldMaterial;
    G4Material*        absorberMaterial;
    G4Material*        gapMaterial;

    G4Box*             layerSolid;
    G4Box*             gapSolid;

    G4LogicalVolume*   worldLogical;
    G4LogicalVolume*   calorLogical[3];
    G4LogicalVolume*   layerLogical[3];
    G4LogicalVolume*   gapLogical[3];

    G4VPhysicalVolume* worldPhysical;
    G4VPhysicalVolume* calorPhysical[3];
    G4PVReplica*       layerPhysical[3];
    G4VPhysicalVolume* gapPhysical[3];

    G4bool             serial;

    ExN07DetectorMessenger* detectorMessenger; 
      
};

#endif

