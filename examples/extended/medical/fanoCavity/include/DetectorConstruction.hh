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
/// \file medical/fanoCavity/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
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
  
     void SetWallThickness   (G4double);
     void SetWallMaterial    (const G4String&);
       
     void SetCavityThickness (G4double);
     void SetCavityRadius    (G4double);           
     void SetCavityMaterial  (const G4String&);
          
     virtual G4VPhysicalVolume* Construct();
     
     G4double     GetWallThickness()   {return fWallThickness;};
     G4double     GetWallRadius()      {return fWallRadius;};           
     G4Material*  GetWallMaterial()    {return fWallMaterial;};
     G4VPhysicalVolume* GetWall()      {return fWall;};
                     
     G4double     GetCavityThickness() {return fCavityThickness;};
     G4double     GetCavityRadius()    {return fCavityRadius;};           
     G4Material*  GetCavityMaterial()  {return fCavityMaterial;}; 
     G4VPhysicalVolume* GetCavity()    {return fCavity;};
          
     G4double     GetTotalThickness()  {return fTotalThickness;};     
     
     void         PrintParameters();
                       
  private:

     void         DefineMaterials();
   
     G4double            fWallThickness;
     G4double            fWallRadius;     
     G4Material*         fWallMaterial;
     G4VPhysicalVolume*  fWall;
     
     G4double            fCavityThickness;
     G4double            fCavityRadius;     
     G4Material*         fCavityMaterial;
     G4VPhysicalVolume*  fCavity;
     
     G4double            fTotalThickness;

     DetectorMessenger*  fDetectorMessenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

