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
// $Id: Tst20DetectorConstruction.hh,v 1.6 2007-11-09 18:32:59 pia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 


#ifndef Tst20DetectorConstruction_h
#define Tst20DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Tst20DetectorMessenger;
class Tst20CalorimeterSD;


class Tst20DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Tst20DetectorConstruction();
   ~Tst20DetectorConstruction();

  public:
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);          
      
     void SetAbsorberZpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);

     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void CleanGeometry();
     void UpdateGeometry();
     
  public:
  
     void PrintCalorParameters(); 
                    
     G4Material* GetWorldMaterial() {return worldMaterial;};
     G4double GetWorldSizeZ() {return worldSizeZ;}; 
     G4double GetWorldSizeR() {return worldSizeR;};
     
     G4double GetAbsorberZpos() {return zAbsorber;}; 
     G4double GetZstartAbs() {return zStartAbs;};
     G4double GetZendAbs() {return zEndAbs;};

     G4Material* GetAbsorberMaterial() {return absorberMaterial;};
     G4double GetAbsorberThickness() {return absorberThickness;};      
     G4double GetAbsorberRadius() {return absorberRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
     const G4VPhysicalVolume* GetAbsorber() {return physiAbsorber;};
                 
  private:
     
     G4bool worldChanged;
     G4Material* absorberMaterial;
     G4double absorberThickness;
     G4double absorberRadius;
 
     G4double zAbsorber ;
     G4double zStartAbs;
     G4double zEndAbs ;
     
     G4Material* worldMaterial;
     G4double worldSizeR;
     G4double worldSizeZ;
            
     G4Tubs* solidWorld;    //pointer to the solid World 
     G4LogicalVolume* logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

     G4Tubs* solidAbsorber; //pointer to the solid Absorber
     G4LogicalVolume* logicAbsorber; //pointer to the logical Absorber
     G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
     
     G4UniformMagField* magneticField;      //pointer to the magnetic field
     
     Tst20DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Tst20CalorimeterSD* calorimeterSD;  //pointer to the sensitive detector
      
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};



#endif

