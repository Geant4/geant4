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
// $Id: DetectorConstruction.hh,v 1.1 2006-05-09 14:03:02 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     void SetSize     (G4double);              
     void SetMaterial (G4String);            

     void UpdateGeometry();
     
  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return pWorld;};
     const
     G4VPhysicalVolume* GetBox()        {return pBox;};
                     
     G4double           GetWorldSize()  {return worldSize;};                          
     G4double           GetBoxSize()    {return boxSize;};
     G4Material*        GetMaterial()   {return aMaterial;};
     
     void               PrintParameters();
                       
  private:
  
     G4VPhysicalVolume*    pWorld;  
     G4VPhysicalVolume*    pBox;
     
     G4double              boxSize;
     G4double              worldSize;
     G4Material*           aMaterial;
     G4Material*           wMaterial;     
     
     DetectorMessenger* detectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

