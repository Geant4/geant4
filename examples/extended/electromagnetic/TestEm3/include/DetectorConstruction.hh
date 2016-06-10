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
/// \file electromagnetic/TestEm3/include/DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// $Id: DetectorConstruction.hh 78655 2014-01-14 11:13:41Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4Cache.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class DetectorMessenger;

class G4GlobalMagFieldMessenger;


     const G4int MaxAbsor = 10;                        // 0 + 9  
     
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
  
  virtual G4VPhysicalVolume* Construct();
  virtual void ConstructSDandField();
     
public:
  
  void PrintCalorParameters(); 
                    
  G4double GetWorldSizeX()           {return fWorldSizeX;}; 
  G4double GetWorldSizeYZ()          {return fWorldSizeYZ;};
     
  G4double GetCalorThickness()       {return fCalorThickness;}; 
  G4double GetCalorSizeYZ()          {return fCalorSizeYZ;};
      
  G4int GetNbOfLayers()              {return fNbOfLayers;}; 
     
  G4int       GetNbOfAbsor()             {return fNbOfAbsor;}; 
  G4Material* GetAbsorMaterial(G4int i)  {return fAbsorMaterial[i];};
  G4double    GetAbsorThickness(G4int i) {return fAbsorThickness[i];};      

  const G4VPhysicalVolume* GetphysiWorld()        {return fPhysiWorld;};
  const G4Material*        GetWorldMaterial()     {return fDefaultMaterial;};
  const G4VPhysicalVolume* GetAbsorber(G4int i)   {return fPhysiAbsor[i];};

private:

  G4int              fNbOfAbsor;
  G4Material*        fAbsorMaterial [MaxAbsor];
  G4double           fAbsorThickness[MaxAbsor];

  G4int              fNbOfLayers;
  G4double           fLayerThickness;

  G4double           fCalorSizeYZ;
  G4double           fCalorThickness;

  G4Material*        fDefaultMaterial;
  G4double           fWorldSizeYZ;
  G4double           fWorldSizeX;

  G4Box*             fSolidWorld;
  G4LogicalVolume*   fLogicWorld;
  G4VPhysicalVolume* fPhysiWorld;

  G4Box*             fSolidCalor;
  G4LogicalVolume*   fLogicCalor;
  G4VPhysicalVolume* fPhysiCalor;

  G4Box*             fSolidLayer;
  G4LogicalVolume*   fLogicLayer;
  G4VPhysicalVolume* fPhysiLayer;

  G4Box*             fSolidAbsor[MaxAbsor];
  G4LogicalVolume*   fLogicAbsor[MaxAbsor];
  G4VPhysicalVolume* fPhysiAbsor[MaxAbsor];

  DetectorMessenger* fDetectorMessenger;
  G4Cache<G4GlobalMagFieldMessenger*> fFieldMessenger;
  
private:

  void DefineMaterials();
  void ComputeCalorParameters();
  G4VPhysicalVolume* ConstructCalorimeter();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

