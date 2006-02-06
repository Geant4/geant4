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
// $Id: Em10DetectorConstruction.hh,v 1.14 2006-02-06 16:21:43 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

#ifndef Em10DetectorConstruction_h
#define Em10DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ios.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4Region;
class G4UniformMagField;
class Em10DetectorMessenger;
class Em10CalorimeterSD;
class G4Region;
class Em10Materials;


class Em10DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em10DetectorConstruction();
   ~Em10DetectorConstruction();

  public:
       
     
     void SetAbsorberMaterial (G4String);     
     void SetAbsorberThickness(G4double);     
     void SetAbsorberRadius(G4double);                
     void SetAbsorberZpos(G4double);

     void SetRadiatorMaterial (G4String);     
     void SetRadiatorThickness(G4double);
 
     void SetGasGapThickness(G4double);     
   
     void SetFoilNumber (G4int    i)  {fFoilNumber = i;  };     

     void SetWorldMaterial(G4String);
     void SetWorldSizeZ(G4double);
     void SetWorldSizeR(G4double);
     void SetDetectorSetUp(G4String s) {fSetUp = s;};


     void SetMagField(G4double);
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
  public:
  
     void PrintGeometryParameters(); 
                    
     G4Material* GetWorldMaterial()    {return fWorldMaterial;};
     G4double GetWorldSizeZ()          {return fWorldSizeZ;}; 
     G4double GetWorldSizeR()          {return fWorldSizeR;};
     
     G4double GetAbsorberZpos()        {return fAbsorberZ;}; 

     G4Material* GetAbsorberMaterial()  {return fAbsorberMaterial;};
     G4double    GetAbsorberThickness() {return fAbsorberThickness;};      
     G4double GetAbsorberRadius()       {return fAbsorberRadius;};
     
     const G4VPhysicalVolume* GetphysiWorld() {return fPhysicsWorld;};           
     const G4VPhysicalVolume* GetAbsorber()   {return fPhysicsAbsorber;};
     G4LogicalVolume* GetLogicalAbsorber()    {return fLogicAbsorber;};

     G4LogicalVolume* GetLogicalRadiator()    {return fLogicRadiator;};
     G4double         GetFoilThick()          {return fRadThickness;};      
     G4double         GetGasThick()           {return fGasGap;};      
     G4int            GetFoilNumber()         {return fFoilNumber;};      
     G4Material* GetFoilMaterial()  {return fFoilMat;};
     G4Material* GetGasMaterial()  {return fGasMat;};
      
private:
    
  G4VPhysicalVolume* ConstructDetectorXTR(); 

    
  G4VPhysicalVolume* SimpleSetUpALICE();     
  G4VPhysicalVolume* SetUpHarris73(); 



  void TestOld();    
                
private:
     
  G4bool             fWorldChanged;
  G4Material*        fAbsorberMaterial;
  G4double           fAbsorberThickness;
  G4double           fAbsorberRadius;

  G4Material*        fWindowMat ;
  G4double           fWindowThick ;

  G4Material*        fElectrodeMat ;
  G4double           fElectrodeThick ;

  G4Material*        fGapMat ;
  G4double           fGapThick ;
 
  G4double           fAbsorberZ ;
  //  G4double           zstartAbs , zendAbs;
  G4String           fSetUp;
     
  G4Material*        fWorldMaterial;
  G4double           fWorldSizeR;
  G4double           fWorldSizeZ;
            
  G4Box*             fSolidWorld;    //pointer to the solid World 
  G4LogicalVolume*   fLogicWorld;    //pointer to the logical World
  G4VPhysicalVolume* fPhysicsWorld;    //pointer to the physical World

  // TR radiator volumes and dimensions
          
  G4Box*             fSolidRadSlice;   // pointer to the solid  z-slice 
  G4LogicalVolume*   fLogicRadSlice;   // pointer to the logical z-slide
  G4VPhysicalVolume* fPhysicRadSlice;  // pointer to the physical z-slide

  G4Box*             fSolidRadRing;    // pointer to the solid  R-slice 
  G4LogicalVolume*   fLogicRadRing;    // pointer to the logical R-slide
  G4VPhysicalVolume* fPhysicRadRing;   // pointer to the physical R-slide

  G4Box*             fSolidRadiator;
  G4LogicalVolume*   fLogicRadiator; 
  G4VPhysicalVolume* fPhysicsRadiator;

  G4Material* fRadiatorMat;        // pointer to the mixed TR radiator material
  G4Material* fFoilMat;            // pointer to the TR foil radiator material
  G4Material* fGasMat;             // pointer to the TR gas radiator material

  G4double fRadThickness ;
  G4double fGasGap       ;
  G4double foilGasRatio  ;

  G4int fFoilNumber ;

  G4double fDetThickness ;
  G4double fDetLength    ;
  G4double fDetGap       ;

  G4double fStartR       ;
  G4double fStartZ       ;

  G4int fModuleNumber ;   // the number of Rad-Det modules

  G4double fRadThick;
  G4double fRadZ;
  G4double fWindowZ;
  G4double fGapZ;
  G4double fElectrodeZ;

  G4Box*             fSolidAbsorber; //pointer to the solid Absorber
  G4LogicalVolume*   fLogicAbsorber; //pointer to the logical Absorber
  G4VPhysicalVolume* fPhysicsAbsorber; //pointer to the physical Absorber
     
  G4UniformMagField* fMagField;      //pointer to the magnetic field

  // G4double fElectronCut, fGammaCut, fPositronCut;
       
  Em10DetectorMessenger* fDetectorMessenger;  //pointer to the Messenger
  Em10CalorimeterSD*     fCalorimeterSD;  //pointer to the sensitive detector
  G4Region*             fRegGasDet;
  G4Region*             fRadRegion;
  Em10Materials*        fMat;  

};

#endif







