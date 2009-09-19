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
// $Id: DetectorConstruction.hh,v 1.1 2009-09-19 16:09:44 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"

#include "globals.hh"

class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

    DetectorConstruction();
   ~DetectorConstruction();

  public:

     void SetMaterialScatter (G4String);
     void SetThicknessScatter(G4double);

     G4VPhysicalVolume* Construct();

     void UpdateGeometry();

  public:

     void PrintGeometry();

     G4double    GetThicknessWorld()         {return thickness_World;};
     G4double    GetRadiusWorld()            {return radius_World;};     
     G4double    GetThicknessFrame()         {return thickness_Frame;};     
  
     G4Material* GetMaterialScatter()        {return material_ScatterFoil;};
     G4double    GetThicknessScatter()       {return thickness_ScatterFoil;};
          
     const G4VPhysicalVolume* GetpvolWorld() {return pvol_World;};
     const G4VPhysicalVolume* GetpvolFrame() {return pvol_Frame;};

  private:
  
     G4Material*        material_World;
     G4double           radius_World;     
     G4double           thickness_World;

     G4Material*        material_Frame;
     G4double           thickness_Frame;
     G4double           zfront_Frame;
     
     G4Material*        material_ExitWindow;
     G4double           thickness_ExitWindow;
     G4double           zfront_ExitWindow;

     G4Material*        material_ScatterFoil;
     G4double           thickness_ScatterFoil;
     G4double           zfront_ScatterFoil;

     G4Material*        material_MonitorChbr;
     G4double           thickness_MonitorChbr;
     G4double           zfront_MonitorChbr;

     G4Material*        material_Bag;
     G4double           thickness_Bag;
     G4double           zfront_Bag;

     G4Material*        material_Gas;
     G4double           thickness_Gas;

     G4Material*        material_Ring;
     G4double           thickness_Ring;
     G4double           inneradius_Ring;
     

     G4VPhysicalVolume* pvol_World;
     G4VPhysicalVolume* pvol_Frame;
     
     DetectorMessenger* detectorMessenger;
      
  private:
    
     void DefineMaterials();
     void GeometryParameters();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

