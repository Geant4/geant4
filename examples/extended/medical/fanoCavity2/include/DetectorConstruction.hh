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
/// \file medical/fanoCavity2/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 103507 2017-04-11 14:15:33Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4NistManager.hh"
#include "globals.hh"

class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     void SetWallMaterial    (G4String);
       
     void SetWallThickness   (G4double);       
     void SetCavityThickness (G4double);            
      
     void SetWorldRadius     (G4double);
                 
     virtual G4VPhysicalVolume* Construct();
     void               UpdateGeometry();
     
  public:
    
     G4double     GetWallThickness()   {return fWallThickness;};
     G4Material*  GetWallMaterial()    {return fWallMaterial;};
     G4VPhysicalVolume* GetWall()      {return fWall;};
                     
     G4double     GetCavityThickness() {return fCavityThickness;};
     G4Material*  GetCavityMaterial()  {return fCavityMaterial;}; 
     G4VPhysicalVolume* GetCavity()    {return fCavity;};
          
     G4double     GetWorldThickness()  {return fWorldThickness;};
     G4double     GetWorldRadius()     {return fWorldRadius;};        
     
     void         PrintParameters();
                       
  private:
   
     G4double            fWallThickness;
     G4Material*         fWallMaterial;
     G4VPhysicalVolume*  fWall;
     
     G4double            fCavityThickness;
     G4Material*         fCavityMaterial;
     G4VPhysicalVolume*  fCavity;
     
     G4double            fWorldThickness;
     G4double            fWorldRadius;

     DetectorMessenger* fDetectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

