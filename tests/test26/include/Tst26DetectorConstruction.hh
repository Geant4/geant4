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
// $Id: Tst26DetectorConstruction.hh,v 1.5 2006-06-29 21:53:05 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
/////////////////////////////////////////////////////////////////////////
//
// test26: Cut per region physics
//
// Created: 31.01.03 V.Ivanchenko
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef Tst26DetectorConstruction_h
#define Tst26DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4VPhysicalVolume.hh"

class G4Tubs;
class G4LogicalVolume;
class G4UniformMagField;
class Tst26DetectorMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Tst26DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    Tst26DetectorConstruction();
   ~Tst26DetectorConstruction();

  public:
     
     void SetEcalMaterial(const G4String&);
     void SetAbsMaterial(const G4String&);
     void SetEcalLength (G4double val)   {ecalLength = val;};
     void SetEcalWidth  (G4double val)   {ecalWidth = val;};      
     void SetVertexLength (G4double val) {vertexLength = val;};
     void SetPadLength  (G4double val)   {padLength = val;};      
     void SetPadWidth  (G4double val)    {padWidth = val;};      
     void SetAbsLength(G4double val)     {absLength = val;};
     
     G4VPhysicalVolume* Construct();

     void UpdateGeometry();
     
     G4double GetWorldSizeZ() {return worldZ;}
     
  private:
     
     G4double ecalLength;    
     G4double ecalWidth;
     G4double vertexLength;
     G4double padLength;
     G4double padWidth;
     G4double absLength;
     G4double worldZ;
  
     G4Material* calMaterial;          
     G4Material* vertMaterial;          
     G4Material* absMaterial;          
     G4Material* worldMaterial;
     G4Material* yorkMaterial;

     G4LogicalVolume* logicC;
     G4LogicalVolume* logicA1;
     G4LogicalVolume* logicA2;
     G4LogicalVolume* logicA3;
     G4LogicalVolume* logicA4;
        
     Tst26DetectorMessenger* detectorMessenger;  //pointer to the Messenger   
      
  private:
    
     void DefineMaterials();
     G4VPhysicalVolume* ConstructVolumes();     
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

#endif

