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
/// \file GB03DetectorConstruction.hh
/// \brief Definition of the GB03DetectorConstruction class

#ifndef GB03DetectorConstruction_h
#define GB03DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVReplica;
class G4Material;
class G4Box;
class GB03DetectorMessenger;

class GB03DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    GB03DetectorConstruction();
    virtual ~GB03DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
  
    void PrintCalorParameters() const;
    void SetAbsorberMaterial(G4String materialChoice);     
    G4String GetAbsorberMaterial() const;
    void SetGapMaterial(G4String materialChoice);     
    G4String GetGapMaterial() const;
    void SetNumberOfLayers(G4int nl);
    static G4int GetNumberOfLayers() { return fNumberOfLayers; }

    G4int GetVerboseLevel() const { return  fVerboseLevel; }
    void SetVerboseLevel(G4int val) { fVerboseLevel = val; }
     
  private:
    void DefineMaterials();
    void SetupGeometry();
    void SetupDetectors();
    void SetupBiasing();

    // data members
    static G4int       fNumberOfLayers;

    G4double           fTotalThickness; /// total thinkness of one calorimeter
    G4double           fLayerThickness; /// = fTotalThickness / fNumberOfLayers

    G4bool             fConstructed;
    static G4ThreadLocal G4bool fConstructedSDandField;
  
    G4String           fCalName;

    G4Material*        fWorldMaterial;
    G4Material*        fAbsorberMaterial;
    G4Material*        fGapMaterial;

    G4Box*             fLayerSolid;
    G4Box*             fGapSolid;

    G4LogicalVolume*   fWorldLogical;
    G4LogicalVolume*   fCalorLogical;
    G4LogicalVolume*   fLayerLogical;
    G4LogicalVolume*   fGapLogical;

    G4VPhysicalVolume* fWorldPhysical;
    G4VPhysicalVolume* fCalorPhysical;
    G4PVReplica*       fLayerPhysical;
    G4VPhysicalVolume* fGapPhysical;

    GB03DetectorMessenger* fDetectorMessenger; 
    
    G4int              fVerboseLevel;
      
};

#endif

