// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2DetectorConstruction.hh,v 1.1 1999-10-11 15:08:36 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em2DetectorConstruction_h
#define Em2DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

class G4Tubs;
class G4LogicalVolume;
class G4UniformMagField;
class Em2DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em2DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em2DetectorConstruction();
   ~Em2DetectorConstruction();

  public:
     
     void SetMaterial(G4String);
     void SetLBining (G4ThreeVector);
     void SetRBining (G4ThreeVector);      
     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     const 
     G4VPhysicalVolume* GetEcal() {return physiEcal;};
     G4Material*    GetMaterial() {return myMaterial;};

     G4int    GetnLtot()      {return nLtot;};
     G4int    GetnRtot()      {return nRtot;};
     G4double GetdLradl()     {return dLradl;}; 
     G4double GetdRradl()     {return dRradl;};
     G4double GetfullLength() {return EcalLength;};
     G4double GetfullRadius() {return EcalRadius;};   
                 
  private:
     
     G4int    nLtot,  nRtot;          // nb of bins: longitudinal and radial
     G4double dLradl, dRradl;         // bin thickness (in radl unit)
     
     G4Material* myMaterial;          //pointer to the material
     G4UniformMagField* magField;     //pointer to the mag field
                                     
     G4double EcalLength;             //full length of the Calorimeter
     G4double EcalRadius;             //radius  of the Calorimeter
       
     G4Tubs*            solidEcal;    //pointer to the solid calorimeter 
     G4LogicalVolume*   logicEcal;    //pointer to the logical calorimeter
     G4VPhysicalVolume* physiEcal;    //pointer to the physical calorimeter
          
     G4Tubs*            solidSlice;   //pointer to the solid  L-slice 
     G4LogicalVolume*   logicSlice;   //pointer to the logical L-slide
     G4VPhysicalVolume* physiSlice;   //pointer to the physical L-slide
     
     G4Tubs*            solidRing;    //pointer to the solid  R-slice 
     G4LogicalVolume*   logicRing;    //pointer to the logical R-slide
     G4VPhysicalVolume* physiRing;    //pointer to the physical R-slide
     
     Em2DetectorMessenger* detectorMessenger;  //pointer to the Messenger object  
      
  private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#endif

