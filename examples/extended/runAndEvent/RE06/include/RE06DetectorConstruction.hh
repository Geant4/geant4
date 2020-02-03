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
/// \file RE06/include/RE06DetectorConstruction.hh
/// \brief Definition of the RE06DetectorConstruction class
//
// 

#ifndef RE06DetectorConstruction_h
#define RE06DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4PVReplica;
class G4Material;
class G4Box;
class RE06DetectorMessenger;

class RE06DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    RE06DetectorConstruction();
    virtual ~RE06DetectorConstruction();

    virtual G4VPhysicalVolume* Construct();
  
    void ConstructSDandField();
  
    void PrintCalorParameters() const;
    void SetAbsorberMaterial(G4String materialChoice);     
    G4String GetAbsorberMaterial() const;
    void SetGapMaterial(G4String materialChoice);     
    G4String GetGapMaterial() const;
    void SetSerialGeometry(G4bool ser);
    void SetNumberOfLayers(G4int nl);
    G4int GetNumberOfLayers() const { return fNumberOfLayers; }
    G4bool IsSerial() const { return fSerial; }

    void  AddMaterial();
  
    G4int GetVerboseLevel() const { return  fVerboseLevel; }
    void SetVerboseLevel(G4int val) { fVerboseLevel = val; }
     
  private:
    void DefineMaterials();
    void SetupGeometry();
    void SetupDetectors();

    // data members
    G4int              fNumberOfLayers;

    G4double           fTotalThickness; ///< total thinkness of one calorimeter
    G4double           fLayerThickness; ///< = fTotalThickness / fNumberOfLayers

    G4bool             fConstructed;
    static G4ThreadLocal G4bool fConstructedSDandField;
  
    G4String           fCalName[3];

    G4Material*        fWorldMaterial;
    G4Material*        fAbsorberMaterial;
    G4Material*        fGapMaterial;

    G4Box*             fLayerSolid;
    G4Box*             fGapSolid;

    G4LogicalVolume*   fWorldLogical;
    G4LogicalVolume*   fCalorLogical[3];
    G4LogicalVolume*   fLayerLogical[3];
    G4LogicalVolume*   fGapLogical[3];

    G4VPhysicalVolume* fWorldPhysical;
    G4VPhysicalVolume* fCalorPhysical[3];
    G4PVReplica*       fLayerPhysical[3];
    G4VPhysicalVolume* fGapPhysical[3];

    G4bool             fSerial;

    RE06DetectorMessenger* fDetectorMessenger; 
    
    G4int              fVerboseLevel;
      
};


#endif

