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
// $Id: DetectorConstruction.hh,v 1.13 2006-06-29 16:51:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4UniformMagField;
class DetectorMessenger;

     const G4int MaxAbsor = 10;			// 0 + 9  
     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  
  DetectorConstruction();
  ~DetectorConstruction();

public:
  
  void SetNbOfAbsor     (G4int);      
  void SetAbsorMaterial (G4int,const G4String&);     
  void SetAbsorThickness(G4int,G4double);
          
  void SetWorldMaterial (const G4String&);
  void SetCalorSizeYZ   (G4double);          
  void SetNbOfLayers    (G4int);   

  void SetMagField   (G4double);
     
  G4VPhysicalVolume* Construct();

  void UpdateGeometry();
     
public:
  
  void PrintCalorParameters(); 
                    
  G4double GetWorldSizeX()           {return WorldSizeX;}; 
  G4double GetWorldSizeYZ()          {return WorldSizeYZ;};
     
  G4double GetCalorThickness()       {return CalorThickness;}; 
  G4double GetCalorSizeYZ()          {return CalorSizeYZ;};
      
  G4int GetNbOfLayers()              {return NbOfLayers;}; 
     
  G4int       GetNbOfAbsor()             {return NbOfAbsor;}; 
  G4Material* GetAbsorMaterial(G4int i)  {return AbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return AbsorThickness[i];};      

  const G4VPhysicalVolume* GetphysiWorld()        {return physiWorld;};
  const G4Material*        GetWorldMaterial()     {return defaultMaterial;};
  const G4VPhysicalVolume* GetAbsorber(G4int i)   {return physiAbsor[i];};

private:

  G4int              NbOfAbsor;
  G4Material*        AbsorMaterial [MaxAbsor];
  G4double           AbsorThickness[MaxAbsor];

  G4int              NbOfLayers;
  G4double           LayerThickness;

  G4double           CalorSizeYZ;
  G4double           CalorThickness;

  G4Material*        defaultMaterial;
  G4double           WorldSizeYZ;
  G4double           WorldSizeX;

  G4Box*             solidWorld;
  G4LogicalVolume*   logicWorld;
  G4VPhysicalVolume* physiWorld;

  G4Box*             solidCalor;
  G4LogicalVolume*   logicCalor;
  G4VPhysicalVolume* physiCalor;

  G4Box*             solidLayer;
  G4LogicalVolume*   logicLayer;
  G4VPhysicalVolume* physiLayer;

  G4Box*             solidAbsor[MaxAbsor];
  G4LogicalVolume*   logicAbsor[MaxAbsor];
  G4VPhysicalVolume* physiAbsor[MaxAbsor];

  G4UniformMagField* magField;

  DetectorMessenger* detectorMessenger;

private:

  void DefineMaterials();
  void ComputeCalorParameters();
  G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

