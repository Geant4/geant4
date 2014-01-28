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
// $Id$
//
// 

#ifndef Tst14DetectorConstruction_h
#define Tst14DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "globals.hh"
#include "tls.hh"

class G4Box;
class G4Tubs;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class Tst14DetectorMessenger;
class Tst14CalorimeterSD;

class Tst14DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  Tst14DetectorConstruction();
  virtual ~Tst14DetectorConstruction();

public:
  
  void SetAbsorberMaterial (G4String);     
  void SetAbsorberThickness(G4double);     
  void SetAbsorberRadius(G4double);          
  
  void SetAbsorberZpos(G4double);
  
  void SetWorldMaterial(G4String);
  void SetWorldSizeZ(G4double);
  void SetWorldSizeR(G4double);
    
  G4VPhysicalVolume* Construct();

  void ConstructSDandField();

  void CleanGeometry();
  void UpdateGeometry();
  
public:
  
  void PrintCalorParameters(); 
  
  G4Material* GetWorldMaterial()    {return WorldMaterial;};
  G4double GetWorldSizeZ()          {return WorldSizeZ;}; 
  G4double GetWorldSizeR()          {return WorldSizeR;};
  
  G4double GetAbsorberZpos()        {return zAbsorber;}; 
  G4double GetzstartAbs() const       {return zstartAbs;};
  G4double GetzendAbs()             {return zendAbs;};
  
  G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};
  G4double    GetAbsorberThickness() const {return AbsorberThickness;};      
  G4double GetAbsorberRadius() const      {return AbsorberRadius;};
  
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
  const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
  
private:
  
  //G4bool             worldchanged;
  G4Material*        AbsorberMaterial;
  G4double           AbsorberThickness;
  G4double           AbsorberRadius;
  
  G4double           zAbsorber ;
  G4double           zstartAbs , zendAbs ;
  
  G4Material*        WorldMaterial;
  G4double           WorldSizeR;
  G4double           WorldSizeZ;
  
  G4Tubs*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Tubs*             solidAbsorber; //pointer to the solid Absorber
  G4LogicalVolume*   logicAbsorber; //pointer to the logical Absorber
  G4VPhysicalVolume* physiAbsorber; //pointer to the physical Absorber
  
  Tst14DetectorMessenger* detectorMessenger;  //pointer to the Messenger

  G4Cache<Tst14CalorimeterSD*> calorimeterSD;  //pointer to the sensitive detector
  //Automatic management of the UI to the mag field
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
  
private:
  
  void DefineMaterials();
  void ComputeCalorParameters();
};



#endif

