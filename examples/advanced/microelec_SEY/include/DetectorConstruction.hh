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
/// \file electromagnetic/TestEm14/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 84208 2014-10-10 14:44:50Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "MicroElecSdSey.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Sphere.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4UserLimits.hh"
#include "G4VisAttributes.hh"



class G4LogicalVolume;
class G4Material;
class DetectorMessenger;
class G4UserLimits;
class G4Region;
class ElectroMagneticFieldSetup;
const G4int kMaxAbsor = 1;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
	  G4Region* GetTargetRegion() { return fRegion; }

     G4VPhysicalVolume* Construct() override;
	 void ConstructSDandField() override;

     void SetSize     (G4double);              
     void SetWidth	  (G4double);
	 void SetMaterial (G4String);
	 void SetSizeSurface(G4double);
	 void SetMaterialSurface(G4String);
	 void SetSizeLayer1(G4double);
	 void SetMaterialLayer1(G4String);
	 void SetSizeLayer2(G4double);
	 void SetMaterialLayer2(G4String);
	 void SetSizeLayer3(G4double);
	 void SetMaterialLayer3(G4String);
	 void SetSizeLayer4(G4double);
	 void SetMaterialLayer4(G4String);


	 void UpdateGeometry();
	 G4Material* DefaultMaterial;
	 G4Material* My_TiN;
	 G4Material* My_BN;
	 G4Material* My_SiO2;
	 G4Material* My_Kapton;

  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return fPBox;};           
                    
     G4double           GetSize()       {return fBoxSize;};
	 G4double           GetWidth()		{ return fBoxWidth; };
     G4Material*        GetMaterial()   {return fMaterial;};
	 G4double           GetSizeSurface() { return fBoxSizeSurface; };
	 G4Material*        GetMaterialSurface() { return fMaterialSurface; };
	 G4double           GetSizeLayer1() { return fBoxSizeLayer1; };
	 G4Material*        GetMaterialLayer1() { return fMaterialLayer1; };
	 G4double           GetSizeLayer2() { return fBoxSizeLayer2; };
	 G4Material*        GetMaterialLayer2() { return fMaterialLayer2; };
	 G4double           GetSizeLayer3() { return fBoxSizeLayer3; };
	 G4Material*        GetMaterialLayer3() { return fMaterialLayer3; };
	 G4double           GetSizeLayer4() { return fBoxSizeLayer4; };
	 G4Material*        GetMaterialLayer4() { return fMaterialLayer4; };

	
	 //SDTrajectoires* GetSensitiveDetector() { return pSDTrajectoires; };
	 
     void               PrintParameters();
                       
  private:
  
	  G4Region*          fRegion;

	 G4double		fWorldSizeX;
	 G4double		fWorldSizeY;
	 G4double		fWorldSizeZ;
	 G4Material*	World_Material;

     G4Box*				fSolidWorld;
	 G4LogicalVolume*   fLogicWorld;
	 G4VPhysicalVolume* fPhysiWorld;
	 
	 G4double			WorldDim;
	 G4double			WorldRay;
	 
	 // substrate
	 G4Box*					fSBox;
     G4LogicalVolume*      fLBox;
	 G4VPhysicalVolume*    fPBox;
     G4double              fBoxSize;
	 G4double              fBoxWidth;
     G4Material*           fMaterial;
	 // surface
	 G4Box*					fSBoxSurface;
	 G4VPhysicalVolume*    fPBoxSurface;
	 G4LogicalVolume*      fLBoxSurface;
	 G4double              fBoxSizeSurface;
	 G4Material*           fMaterialSurface;
	 // Layer 1 just below the surface
	 G4Box*					fSBoxLayer1;
	 G4VPhysicalVolume*    fPBoxLayer1;
	 G4LogicalVolume*      fLBoxLayer1;
	 G4double              fBoxSizeLayer1;
	 G4Material*           fMaterialLayer1;
	 // Layer 2 just below the layer1
	 G4Box*					fSBoxLayer2;
	 G4VPhysicalVolume*    fPBoxLayer2;
	 G4LogicalVolume*      fLBoxLayer2;
	 G4double              fBoxSizeLayer2;
	 G4Material*           fMaterialLayer2;
	 // Layer 3 just below the layer2
	 G4Box*					fSBoxLayer3;
	 G4VPhysicalVolume*    fPBoxLayer3;
	 G4LogicalVolume*      fLBoxLayer3;
	 G4double              fBoxSizeLayer3;
	 G4Material*           fMaterialLayer3;
	 // Layer 4 just below the layer3 on top of the substrate
	 G4Box*					fSBoxLayer4;
	 G4VPhysicalVolume*    fPBoxLayer4;
	 G4LogicalVolume*      fLBoxLayer4;
	 G4double              fBoxSizeLayer4;
	 G4Material*           fMaterialLayer4;

	 G4Sphere*              Detector_1_box;
	 G4LogicalVolume*		Detector_1_log; 
	 G4VPhysicalVolume*		Detector_1_phys;  
	 G4Material* detectorMat;
	 G4double DetectorRay;
	 

     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();


};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

