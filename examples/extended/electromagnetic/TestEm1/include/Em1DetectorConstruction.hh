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
// $Id: Em1DetectorConstruction.hh,v 1.3 2001-07-11 09:57:12 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#ifndef Em1DetectorConstruction_h
#define Em1DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

class G4LogicalVolume;
class G4Material;
class G4UniformMagField;
class Em1DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class Em1DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Em1DetectorConstruction();
   ~Em1DetectorConstruction();

  public:
  
     G4VPhysicalVolume* Construct();
     
     void SetSize     (G4double);              
     void SetMaterial (G4String);            
     void SetMagField (G4double);

     void UpdateGeometry();
     
  public:
  
     const
     G4VPhysicalVolume* GetWorld()      {return pBox;};           
                    
     G4double           GeSize()        {return BoxSize;};      
     G4Material*        GetMaterial()   {return aMaterial;};
     
     void               PrintParameters();
                       
  private:
  
     G4VPhysicalVolume*    pBox;
     G4LogicalVolume*      lBox;
     
     G4double              BoxSize;
     G4Material*           aMaterial;     
     G4UniformMagField*    magField;
     
     Em1DetectorMessenger* detectorMessenger;

  private:
    
     void               DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....


#endif

