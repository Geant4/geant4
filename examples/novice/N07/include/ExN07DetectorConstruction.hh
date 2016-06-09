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
// $Id: ExN07DetectorConstruction.hh,v 1.2 2004/11/17 03:07:30 kurasige Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
// 

#ifndef ExN07DetectorConstruction_h
#define ExN07DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVParameterised;
class G4Material;
class G4Box;
class ExN07GapParameterisation;
class ExN07DetectorMessenger;

class ExN07DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN07DetectorConstruction();
    virtual ~ExN07DetectorConstruction();

  public:
    virtual G4VPhysicalVolume* Construct();
     
  public:
    void CreateMaterial(G4String materialChoice);
    void UpdateGeometry();
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

  private:
    G4int              numberOfLayers;

    G4Material*        worldMaterial;
    G4Material*        absorberMaterial;
    G4Material*        gapMaterial;

    G4Box*             gapSolid;

    G4LogicalVolume*   worldLogical;
    G4LogicalVolume*   calorLogical[3];
    G4LogicalVolume*   layerLogical[3];

    G4VPhysicalVolume* worldPhysical;
    G4VPhysicalVolume* calorPhysical[3];
    G4PVParameterised* layerPhysical[3];

    G4bool             serial;

    ExN07GapParameterisation* gapParam;
    ExN07DetectorMessenger* detectorMessenger; 
      
};

#endif

