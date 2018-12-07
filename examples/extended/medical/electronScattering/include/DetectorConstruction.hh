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
/// \file medical/electronScattering/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
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

     virtual G4VPhysicalVolume* Construct();

     void UpdateGeometry();

  public:

     void PrintGeometry();

     G4double    GetThicknessWorld()         {return fThickness_World;};
     G4double    GetRadiusWorld()            {return fRadius_World;};     
     G4double    GetThicknessFrame()         {return fThickness_Frame;};     
  
     G4Material* GetMaterialScatter()        {return fMaterial_ScatterFoil;};
     G4double    GetThicknessScatter()       {return fThickness_ScatterFoil;};

     G4double    GetZdist_foil_detector() 
                                   {return fThickness_Frame-fZfront_ScatterFoil
                                         -0.5*fThickness_ScatterFoil;};
          
     const G4VPhysicalVolume* GetpvolWorld() {return fPvol_World;};
     const G4VPhysicalVolume* GetpvolFrame() {return fPvol_Frame;};

  private:
  
     G4Material*        fMaterial_World;
     G4double           fRadius_World;     
     G4double           fThickness_World;

     G4Material*        fMaterial_Frame;
     G4double           fThickness_Frame;
     G4double           fZfront_Frame;
     
     G4Material*        fMaterial_ExitWindow;
     G4double           fThickness_ExitWindow;
     G4double           fZfront_ExitWindow;

     G4Material*        fMaterial_ScatterFoil;
     G4double           fThickness_ScatterFoil;
     G4double           fZfront_ScatterFoil;

     G4Material*        fMaterial_MonitorChbr;
     G4double           fThickness_MonitorChbr;
     G4double           fZfront_MonitorChbr;

     G4Material*        fMaterial_Bag;
     G4double           fThickness_Bag;
     G4double           fZfront_Bag;

     G4Material*        fMaterial_Gas;
     G4double           fThickness_Gas;

     G4Material*        fMaterial_Ring;
     G4double           fThickness_Ring;
     G4double           fInnerRadius_Ring;
     

     G4VPhysicalVolume* fPvol_World;
     G4VPhysicalVolume* fPvol_Frame;
     
     DetectorMessenger* fDetectorMessenger;
      
  private:
    
     void DefineMaterials();
     void GeometryParameters();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

