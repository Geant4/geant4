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
// $Id: Tst50DetectorConstruction.hh,v 1.3 2003-01-07 15:29:39 guatelli Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst50DetectorConstruction_h
#define Tst50DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class  Tst50TrackerSD;
class Tst50DetectorMessenger;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst50DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Tst50DetectorConstruction();
   ~Tst50DetectorConstruction();

  public:
     
     void SetTargetMaterial (G4String);     
     void SetTargetThickness(G4double);     
     void SetTargetX(G4double); 
     void SetTargetY(G4double); 
  G4double GetTargetMaterialDensity() {return density;};
          
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintParameters(); 
                    
        
     G4Material* GetTargetMaterial()  {return TargetMaterial;};
     G4double    GetTargetThickness() {return TargetThickness;};      
     
        
   
                 
  private:
     
     G4Material*        TargetMaterial;
     G4double           TargetThickness;
     G4double targetX;
     G4double targetY; 
  G4double density;

     G4Material*        defaultMaterial;
     G4double           WorldSizeYZ;
     G4double           WorldSizeX;
            
     G4Box*             solidWorld;    //pointer to the solid World 
     G4LogicalVolume*   logicWorld;    //pointer to the logical World
     G4VPhysicalVolume* physiWorld;    //pointer to the physical World

   
         
     G4Box*             solidTarget; //pointer to the solid Target
     G4LogicalVolume*   logicTarget; //pointer to the logical Target
     G4VPhysicalVolume* physiTarget; //pointer to the physical Target
     
        
     Tst50DetectorMessenger* detectorMessenger;  //pointer to the Messenger
     Tst50TrackerSD* pTargetSD;  //pointer to the sensitive detector
      
  private:
    
     void DefineMaterials();
     void ComputeParameters();
     G4VPhysicalVolume* ConstructWorld();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif




