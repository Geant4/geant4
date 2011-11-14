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
// This is the *BASIC* version of IORT, a Geant4-based application
//
// Main Authors: G.Russo(a,b), C.Casarino*(c), G.C. Candiano(c), G.A.P. Cirrone(d), F.Romano(d)
// Contributor Authors: S.Guatelli(e)
// Past Authors: G.Arnetta(c), S.E.Mazzaglia(d)
//    
//   (a) Fondazione Istituto San Raffaele G.Giglio, Cefalù, Italy
//   (b) IBFM-CNR , Segrate (Milano), Italy
//   (c) LATO (Laboratorio di Tecnologie Oncologiche), Cefalù, Italy
//   (d) Laboratori Nazionali del Sud of the INFN, Catania, Italy
//   (e) University of Wallongong, Australia
//
//   *Corresponding author, email to carlo.casarino@polooncologicocefalu.it
//////////////////////////////////////////////////////////////////////////////////////////////

#ifndef Collimator70BeamLine_H
#define Collimator70BeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

#include "G4Cons.hh"  
#include "G4SubtractionSolid.hh"  

class G4VPhysicalVolume;
class IORTDetectorConstruction;
class Collimator70BeamLineMessenger;

class Collimator70BeamLine : public G4VUserDetectorConstruction
{
public:

 
   Collimator70BeamLine();
  ~Collimator70BeamLine();

  G4VPhysicalVolume* Construct();  


  void IortBeamLineVacuumSource();
  void IortBeamLineTitaniumWindows();
  void IortBeamLineMonitorChambers();
  void IortBeamLineBlocks() ;
  void IortBeamLineJunctions(); 
  void IortBeamLineFinalCollimator();
 
  
  void SetInnerRadiusFinalCollimatorIORT(G4double);
  // This method allows to change the size of the inner radius of the 
  // final collimator

  void SetOuterRadiusFinalCollimatorIORT(G4double);
  // This method allows to change the size of the outer radius of the 
  // final collimator
  
  //  void SetFinalCollimatorIORTMaterial(G4String);
  // This method allows to change the material 
  // of the final collimator

  
 
private:
//Collimator70 line dimensions
  void SetDefaultDimensions(); 
  void ConstructCollimator70BeamLine();

 
  Collimator70BeamLineMessenger* collimatorMessenger;  
  G4VPhysicalVolume* physicalTreatmentRoom;
  IORTDetectorConstruction* iortDetectorConstruction; 

  G4Material* kapton;

  G4VisAttributes* blue;
  G4VisAttributes* gray;
  G4VisAttributes* white;
  G4VisAttributes* red;
  G4VisAttributes* yellow;
  G4VisAttributes* green;
  G4VisAttributes* darkGreen;
  G4VisAttributes* darkOrange3;
  G4VisAttributes* skyBlue;

  // Final Collimator IORT
  G4double innerRadiusFinalCollimatorIORT;
  G4double OuterRadiusFinalCollimatorIORT;
  G4Tubs* solidFinalCollimatorIORT; 
  G4LogicalVolume* logicFinalCollimatorIORT;
  G4VPhysicalVolume* physiFinalCollimatorIORT;
  G4Material* finalCollimatorMaterialIORT;

  // Junction 1 FINAL COLLIMATOR IORT
  G4double innerRadiusGiunz1FinalCollIORT;
  G4double OuterRadiusGiunz1FinalCollIORT;
  G4Tubs* solidGiunz1FinalCollIORT; 
  G4LogicalVolume* logicGiunz1FinalCollIORT;  
  G4VPhysicalVolume* physiGiunz1FinalCollIORT;
  G4Material* Giunz1FinalCollMaterialIORT;

  // Junction 2 FINAL COLLIMATOR IORT
  G4double innerRadiusGiunz2FinalCollIORT;
  G4double OuterRadiusGiunz2FinalCollIORT;
  G4Tubs* solidGiunz2FinalCollIORT; 
  G4LogicalVolume* logicGiunz2FinalCollIORT;
  G4VPhysicalVolume* physiGiunz2FinalCollIORT;
  G4Material* Giunz2FinalCollMaterialIORT;

  // Junction 3 FINAL COLLIMATOR IORT
  G4double innerRadiusGiunz3FinalCollIORT;
  G4double OuterRadiusGiunz3FinalCollIORT;
  G4Tubs* solidGiunz3FinalCollIORT; 
  G4LogicalVolume* logicGiunz3FinalCollIORT;  
  G4VPhysicalVolume* physiGiunz3FinalCollIORT;
  G4Material* Giunz3FinalCollMaterialIORT;   

  // Junction 3 FINAL COLLIMATOR intIORT
  G4Cons* solidGiunz3FinalCollIntIORT; 
  G4LogicalVolume* logicGiunz3FinalCollIntIORT;
  G4Material* Giunz3FinalCollMaterialIntIORT;
  G4SubtractionSolid* solidsottrazione;
  G4LogicalVolume* logicsottrazione;
  G4VPhysicalVolume* physiGiunz3FinalCollIntIORT;

  // Junction 4 FINAL COLLIMATOR IORT
  G4double innerRadiusGiunz4FinalCollIORT;
  G4double OuterRadiusGiunz4FinalCollIORT;
  G4Tubs* solidGiunz4FinalCollIORT; 
  G4LogicalVolume* logicGiunz4FinalCollIORT;  
  G4VPhysicalVolume* physiGiunz4FinalCollIORT;
  G4Material* Giunz4FinalCollMaterialIORT;

  // Junction 5 FINAL COLLIMATOR IORT
  G4double innerRadiusGiunz5FinalCollIORT;
  G4double OuterRadiusGiunz5FinalCollIORT;
  G4Tubs* solidGiunz5FinalCollIORT; 
  G4LogicalVolume* logicGiunz5FinalCollIORT;
  G4VPhysicalVolume* physiGiunz5FinalCollIORT;
  G4Material* Giunz5FinalCollMaterialIORT;

  // Block 1 Diameter 30 mm 
  G4double innerRadiusBlocco1IORT;
  G4double OuterRadiusBlocco1IORT;
  G4Tubs* solidBlocco1IORT; 
  G4LogicalVolume* logicBlocco1IORT;
  G4VPhysicalVolume* physiBlocco1IORT;
  G4Material* Blocco1IORTMaterialIORT;

  // Block 2 Diameter 30 mm 
  G4double innerRadiusBlocco2IORT;
  G4double OuterRadiusBlocco2IORT;
  G4Tubs* solidBlocco2IORT; 
  G4LogicalVolume* logicBlocco2IORT;
  G4VPhysicalVolume* physiBlocco2IORT;
  G4Material* Blocco2IORTMaterialIORT;

  // Block 3 Diameter 30 mm 
  G4double innerRadiusBlocco3IORT;
  G4double OuterRadiusBlocco3IORT;
  G4Tubs* solidBlocco3IORT; 
  G4LogicalVolume* logicBlocco3IORT;  
  G4VPhysicalVolume* physiBlocco3IORT;
  G4Material* Blocco3IORTMaterialIORT;

  // Block 4 Internal Diameter 20 mm 
  G4double innerRadiusBlocco20mmIORT;
  G4double OuterRadiusBlocco20mmIORT;
  G4Tubs* solidBlocco20mmIORT; 
  G4LogicalVolume* logicBlocco20mmIORT;  
  G4VPhysicalVolume* physiBlocco20mmIORT;
  G4Material* Blocco20mmIORTMaterialIORT;

  // First Chamber Monitor  1 of 2 
  G4double innerRadiusCM1_1_2IORT;
  G4double OuterRadiusCM1_1_2IORT;
  G4Tubs* solidCM1_1_2IORT; 
  G4LogicalVolume* logicCM1_1_2IORT;  
  G4VPhysicalVolume* physiCM1_1_2IORT;
  G4Material* CM1_1_2IORTMaterialIORT;

  // First Chamber Monitor  2 of 2 
  G4double innerRadiusCM1_2_2IORT;
  G4double OuterRadiusCM1_2_2IORT;
  G4Tubs* solidCM1_2_2IORT; 
  G4LogicalVolume* logicCM1_2_2IORT;  
  G4VPhysicalVolume* physiCM1_2_2IORT;
  G4Material* CM1_2_2IORTMaterialIORT;

  // Second Chamber Monitor  1 of 2 
  G4double innerRadiusCM2_1_2IORT;
  G4double OuterRadiusCM2_1_2IORT;
  G4Tubs* solidCM2_1_2IORT; 
  G4LogicalVolume* logicCM2_1_2IORT; 
  G4VPhysicalVolume* physiCM2_1_2IORT;
  G4Material* CM2_1_2IORTMaterialIORT;

  // Second Chamber Monitor  2 of 2 
  G4double innerRadiusCM2_2_2IORT;
  G4double OuterRadiusCM2_2_2IORT;
  G4Tubs* solidCM2_2_2IORT; 
  G4LogicalVolume* logicCM2_2_2IORT;
  G4VPhysicalVolume* physiCM2_2_2IORT;
  G4Material* CM2_2_2IORTMaterialIORT;

  // Monitor Chambers Cylinder 
  G4double innerRadiusCCMIORT;
  G4double OuterRadiusCCMIORT;
  G4Tubs* solidCCMIORT; 
  G4LogicalVolume* logicCCMIORT;
  G4VPhysicalVolume* physiCCMIORT;
  G4Material* CCMIORTMaterialIORT;

  //  Superior Final Monitor Chamber Part 1
  G4double innerRadiusPFS1IORT;
  G4double OuterRadiusPFS1IORT;
  G4Tubs* solidPFS1IORT; 
  G4LogicalVolume* logicPFS1IORT;
  G4VPhysicalVolume* physiPFS1IORT;
  G4Material* PFS1IORTMaterialIORT;
  
  //  Superior Final Monitor Chamber Part 2
  G4double innerRadiusPFS2IORT;
  G4double OuterRadiusPFS2IORT;
  G4Tubs* solidPFS2IORT; 
  G4LogicalVolume* logicPFS2IORT;
  G4VPhysicalVolume* physiPFS2IORT;
  G4Material* PFS2IORTMaterialIORT;

  //  Superior Final Monitor Chamber Part 3
  G4double innerRadiusPFS3IORT;
  G4double OuterRadiusPFS3IORT;
  G4Tubs* solidPFS3IORT; 
  G4LogicalVolume* logicPFS3IORT;
  G4VPhysicalVolume* physiPFS3IORT;
  G4Material* PFS3IORTMaterialIORT;
  
  //  Titanium Window
  G4double innerRadiusFTIORT;
  G4double OuterRadiusFTIORT;
  G4Tubs* solidFTIORT; 
  G4LogicalVolume* logicFTIORT;  
  G4VPhysicalVolume* physiFTIORT;
  G4Material* FTIORTMaterialIORT;
  
  //  Vacuum Source
  G4double innerRadiusVSIORT;
  G4double OuterRadiusVSFTIORT;
  G4Tubs* solidVSIORT; 
  G4LogicalVolume* logicVSIORT;
  G4VPhysicalVolume* physiVSIORT;
  G4Material* VSIORTMaterialIORT;

 
};
#endif



