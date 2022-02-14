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
//    ****************************************
//    *                                      *
//    *    BrachyDetectorConstruction.hh     *
//    *                                      *
//    ****************************************
//
// S. Guatelli, A. Le, University of Wollongong

#ifndef BrachyDetectorConstruction_H
#define BrachyDetectorConstruction_H 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"

class BrachyDetectorMessenger;
class BrachyFactory;
class G4LogicalVolume;
class G4Material;
class G4Box;
class G4Colour;
class G4VPhysicalVolume;
class G4VPhysicalVolume;

class BrachyDetectorConstruction : public G4VUserDetectorConstruction
{
public:

  explicit BrachyDetectorConstruction();
  ~BrachyDetectorConstruction();

  G4VPhysicalVolume*   Construct() override;  
  void SwitchBrachytherapicSeed(); //Change brachy source through GUI
  void SelectBrachytherapicSeed(G4String val);
  void ConstructPhantom(); 
  void PrintDetectorParameters(); 
  void SetPhantomMaterial(G4String); 

private:
  BrachyFactory* fFactory;
  BrachyDetectorMessenger* fDetectorMessenger;   
  
  G4Box*             fWorld;        //pointer to the solid World 
  G4LogicalVolume*   fWorldLog;     //pointer to the logical World
  G4VPhysicalVolume* fWorldPhys;    //pointer to the physical World

  G4Box*              fPhantom;  //pointer to solid phantom
  G4LogicalVolume*    fPhantomLog; //pointer to logic phantom
  G4VPhysicalVolume*  fPhantomPhys; //pointer to physical phantom
 
  G4double fPhantomSizeX; //Phantom X Size
  G4double fPhantomSizeY; //Phantom Y Size
  G4double fPhantomSizeZ; //Phantom Z Size  
  
  G4double fWorldSizeX; //World X Size
  G4double fWorldSizeY; //World Y Size
  G4double fWorldSizeZ; //World X Size
  
  G4int fDetectorChoice; //Select brachytherapic seed
};
#endif
