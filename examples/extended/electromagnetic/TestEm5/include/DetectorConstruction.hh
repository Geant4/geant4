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
// $Id: DetectorConstruction.hh,v 1.5 2004/06/21 10:57:10 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "globals.hh"

class G4Box;
class G4VPhysicalVolume;
class G4Material;
class G4MaterialCutsCouple;
class G4UniformMagField;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetAbsorberMaterial (G4String);
     void SetAbsorberThickness(G4double);
     void SetAbsorberSizeYZ   (G4double);

     void SetAbsorberXpos(G4double);

     void SetWorldMaterial(G4String);
     void SetWorldSizeX   (G4double);
     void SetWorldSizeYZ  (G4double);

     void SetMagField(G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

  public:

     void PrintCalorParameters();

     G4Material* GetAbsorberMaterial()  {return AbsorberMaterial;};	     
     G4double    GetAbsorberThickness() {return AbsorberThickness;};
     G4double    GetAbsorberSizeYZ()    {return AbsorberSizeYZ;};

     G4double    GetAbsorberXpos()      {return XposAbs;};
     G4double    GetxstartAbs()         {return xstartAbs;};
     G4double    GetxendAbs()           {return xendAbs;};

     G4Material* GetWorldMaterial()     {return WorldMaterial;};
     G4double    GetWorldSizeX()        {return WorldSizeX;};
     G4double    GetWorldSizeYZ()       {return WorldSizeYZ;};

     const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};
     const G4VPhysicalVolume* GetAbsorber()   {return physiAbsorber;};
     const G4MaterialCutsCouple* GetAbsorbMaterialCut()  const
                             {return logicAbsorber->GetMaterialCutsCouple();};
  private:

     G4Material*        AbsorberMaterial;
     G4double           AbsorberThickness;
     G4double           AbsorberSizeYZ;

     G4double           XposAbs;
     G4double           xstartAbs, xendAbs;

     G4Material*        WorldMaterial;
     G4double           WorldSizeX;
     G4double           WorldSizeYZ;

     G4bool             defaultWorld;

     G4Box*             solidWorld;
     G4LogicalVolume*   logicWorld;
     G4VPhysicalVolume* physiWorld;

     G4Box*             solidAbsorber;
     G4LogicalVolume*   logicAbsorber;
     G4VPhysicalVolume* physiAbsorber;

     G4UniformMagField* magField;
     
     DetectorMessenger* detectorMessenger;
      
  private:
    
     void DefineMaterials();
     void ComputeCalorParameters();
     G4VPhysicalVolume* ConstructCalorimeter();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

