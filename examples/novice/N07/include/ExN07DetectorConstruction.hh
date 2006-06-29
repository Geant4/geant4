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
// $Id: ExN07DetectorConstruction.hh,v 1.5 2006-06-29 17:54:42 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

    void  AddMaterial();
  
    G4int GetVerboseLevel() const;
    void SetVerboseLevel(G4int val);
 
     
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
    
    G4int              verboseLevel;
      
};

inline  G4int ExN07DetectorConstruction::GetVerboseLevel() const
{
  return  verboseLevel;
}

inline  void ExN07DetectorConstruction::SetVerboseLevel(G4int val)
{
  verboseLevel = val;
}



#endif

