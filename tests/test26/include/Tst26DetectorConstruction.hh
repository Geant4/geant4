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
// $Id: Tst26DetectorConstruction.hh,v 1.1 2003-01-31 18:43:56 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

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
     
     Em2DetectorMessenger* detectorMessenger;  //pointer to the Messenger   
      
  private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

