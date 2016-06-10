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
// Hadrontherapy advanced example for Geant4
// See more at: https://twiki.cern.ch/twiki/bin/view/Geant4/AdvancedExamplesHadrontherapy

#ifndef LaserDrivenBeamLine_H
#define LaserDrivenBeamLine_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"

class G4VPhysicalVolume;
class HadrontherapyDetectorConstruction;
class LaserDrivenBeamLineMessenger;
class G4Mag_UsualEqRhs;
class G4MagIntegratorStepper;
class G4ChordFinder;
class G4UniformMagField;
class G4MagInt_Driver; 
class G4EqMagElectricField;
//class G4TransportationManager;
class G4FieldManager;
class G4MagneticField; 
class HadrontherapyMagneticField3D;
class HadrontherapyDetectorROGeometry;

class HadrontherapyElectricTabulatedField3D;


class LaserDrivenBeamLine : public G4VUserDetectorConstruction
{
public:

  LaserDrivenBeamLine();
  ~LaserDrivenBeamLine();

  G4VPhysicalVolume* Construct();
  void ConstructSDandField();

  void RemoveESS();	
  void SetFirstCollimatorRadius(G4double value);
  void SetFirstCollimatorThickness(G4double value);
  void SetFirstCollimatorPositionZ(G4double value);
  void SetSecondCollimatorRadius(G4double value);
  void SetSecondCollimatorThickness(G4double value);
  void SetSecondCollimatorPositionZ(G4double value);
  void SetThicknessSlit(G4double value);
  void SetSlitHoleDimensionY(G4double value);
  void SetSlitHoleDimensionZ(G4double value);
  void SetSlitHolePositionZ(G4double value);
  void RemoveQuads();

private:
  void SetDefaultDimensions(); 
  void ConstructLaserDrivenBeamLine();
  void EnergySelectorChamber();
  void Collimator();
  void Magnet_1();
  void Magnet_2();
  void Magnet_3();
  void Magnet_4();
  void Slit();
  void FinalCollimator();
  void ExitPipe();
  void ExitWindow();
  void Exithole();
  void Entrancehole();
  void EntrancePipe();
  void Quadrupole();
  void FaradayCup();
	
  LaserDrivenBeamLineMessenger *laserDrivenMessenger;  
  HadrontherapyDetectorConstruction* hadrontherapydetectorconstruction;

  HadrontherapyDetectorROGeometry* RO;

  // Variables definition for the World
  // (called treatment room)
  G4Box *solidTreatmentRoom;
  G4LogicalVolume *logicTreatmentRoom;
  G4VPhysicalVolume *physicTreatmentRoom;

  // Variables definition for the triplet of quadrupoles
  G4Material *QuadMaterial;

  G4double QuadChamberWallPosX;
  G4double QuadChamberWallPosY;
  G4double QuadChamberWallPosZ;
  G4Box *SQuadChamberWall, *SQuadChamber;
  G4LogicalVolume *LQuadChamberWall, *LQuadChamber;
  G4VPhysicalVolume *PQuadChamberWall, *PQuadChamber; 

  G4double InnerRadiusQuad;
  G4double InnerRadiusTriplet;
  G4double ExternalRadiusQuad;
  G4double FirstQuadThickness;
  G4double SecondQuadThickness;
  G4double ThirdQuadThickness;
  G4double FourthQuadThickness;
  G4double startAngleQuad;
  G4double spanningAngleQuad;
  G4double FirstQuadXPosition;
  G4double FirstQXPosition;
  G4double SecondQuadXPosition;
  G4double SecondQXPosition;
  G4double ThirdQuadXPosition;
  G4double ThirdQXPosition;
  G4double FourthQuadXPosition;
  G4double FourthQXPosition;
  G4double QuadYPosition;
  G4double QYPosition;
  G4double QuadZPosition;
  G4double QZPosition;
  
  G4Tubs *SFirstTriplet, *SSecondTriplet, *SThirdTriplet, *SFourthTriplet;
  G4LogicalVolume *LFirstTriplet, *LSecondTriplet, *LThirdTriplet, *LFourthTriplet;
  G4VPhysicalVolume *PFirstTriplet, *PSecondTriplet, *PThirdTriplet, *PFourthTriplet;

  G4Tubs *solidFirstQuad, *solidSecondQuad, *solidThirdQuad, *solidFourthQuad;
  G4LogicalVolume *logicFirstQuad, *logicSecondQuad, *logicThirdQuad, *logicFourthQuad;
  G4VPhysicalVolume *physicFirstQuad, *physicSecondQuad, *physicThirdQuad, *physicFourthQuad;  

  // Variables definition for the vacuum chamber containing
  // the spectrometer
  G4Material *externalChamberMaterial;
  G4Material *internalChamberMaterial;
	
  G4Box *solidExternalChamber;
  G4LogicalVolume *logicExternalChamber;
  G4VPhysicalVolume *physicExternalChamber; 
	
  G4Box *solidInternalChamber;
  G4LogicalVolume *logicInternalChamber;
  G4VPhysicalVolume *physicInternalChamber; 

  G4double VaccumChamberWallThickness;

  /////////////////// Declaration of Magnetic Variables ///////////////
  G4FieldManager   *pFieldMgr, *pFieldMgrQuadFourth, *pFieldMgrQuadThird, *pFieldMgrQuadSecond, *pFieldMgrQuadFirst;
  G4MagneticField* PurgMagField, *PurgMagFieldQuadFourth, *PurgMagFieldQuadThird, *PurgMagFieldQuadSecond, *PurgMagFieldQuadFirst;
  G4ChordFinder *pChordFinder, *pChordFinderQuadFourth, *pChordFinderQuadThird, *pChordFinderQuadSecond, *pChordFinderQuadFirst;
  G4Mag_UsualEqRhs* fEquation, *fEquationQuadFourth, *fEquationQuadThird, *fEquationQuadSecond, *fEquationQuadFirst;
  G4MagInt_Driver* pIntgrDriver, *pIntgrDriverQuadFourth, *pIntgrDriverQuadThird, *pIntgrDriverQuadSecond, *pIntgrDriverQuadFirst;
  G4MagIntegratorStepper* fstepper, *fstepperQuadFourth, *fstepperQuadThird, *fstepperQuadSecond, *fstepperQuadFirst;
 
  /////////////////// Declaration of Exit Window Variables ///////////////
  G4double InnerRadiusExitWindow;
  G4double ExternalRadiusExitWindow;
  G4double ExitWindowThickness;
  
  G4double ExitWindowXPosition;
  G4double ExitWindowYPosition;
  G4double ExitWindowZPosition;
  
  G4double startAngleExitWindow;
  G4double spanningAngleExitWindow;
  /////////////////// Declaration of Exit Pipe Variables ///////////////
  G4double ExitPipeheight; 
  G4double InnerRadiusExitPipe;
  G4double ExternalRadiusExitPipe;
  
  G4double ExitPipeXPosition;
  G4double ExitPipeYPosition;
  G4double ExitPipeZPosition;
  
  G4double startAngleExitPipe;
  G4double spanningAngleExitPipe;
   /////////////////// Declaration of Entrance Pipe Variables ///////////////
  G4double EntrancePipeheight; 
  G4double InnerRadiusEntrancePipe;
  G4double ExternalRadiusEntrancePipe;
  
  G4double EntrancePipeXPosition;
  G4double EntrancePipeYPosition;
  G4double EntrancePipeZPosition;
  
  G4double startAngleEntrancePipe;
  G4double spanningAngleEntrancePipe;
  /////////////////////// Declaration of Exit hole in vessel Variables ///////////////////
  G4double InnerRadiusExithole;
  G4double ExternalRadiusExithole;
  G4double ExitholeThickness;
  G4double ExitholeXPosition;
  G4double ExitholeYPosition;
  G4double ExitholeZPosition;
  
  G4double startAngleExithole;
  G4double spanningAngleExithole;
  /////////////////////// Declaration of Entrance hole in vessel Variables ///////////////////
  G4double InnerRadiusEntrancehole;
  G4double ExternalRadiusEntrancehole;
  G4double EntranceholeThickness;
  G4double EntranceholeXPosition;
  G4double EntranceholeYPosition;
  G4double EntranceholeZPosition;
  G4double EntranceholeQuadXPosition;

  G4double startAngleEntrancehole;
  G4double spanningAngleEntrancehole;
  ////////////////////////////////////////////////////////////////////////////////////////////////
  
  ////////////////////////////////////////// Relative Distances Declaration  /////////////////////// 
  G4double ExitholeToFinalCollimator;
  G4double FinalCollimatorToMagnet4;
  G4double Magnet4ToMagnet3;
  G4double Magnet3ToMagnet2;
  G4double Magnet2ToMagnet1;
  G4double Magnet1ToFirstCollimator;
	  
  ////////////////////////////////////////// Chamber dimentions  /////////////////////// 
  G4double externalChamberXSize;
  G4double externalChamberYSize;
  G4double externalChamberZSize;
  G4double internalChamberXSize;
  G4double internalChamberYSize;
  G4double internalChamberZSize;
	
  G4double externalChamberXPosition;
  G4double externalChamberYPosition;
  G4double externalChamberZPosition;
	
  // Variables definition of the collimator
  G4double defaultInnerRadiusCollimator;
  G4double innerRadiusCollimator;
	
  G4double defaultThicknessCollimator;
  G4double thicknessCollimator;
	
  G4double defaultOuterRadiusCollimator;
  G4double outerRadiusCollimator; 
	
  G4double defaultStartAngleCollimator;
  G4double startAngleCollimator;
	
  G4double defaultSpanningAngleCollimator;
  G4double spanningAngleCollimator;
	
  G4double defaultCollimatorXPosition;
  G4double collimatorXPosition;
	
  G4double defaultCollimatorYPosition;
  G4double collimatorYPosition;
	
  G4double defaultCollimatorZPosition;
  G4double collimatorZPosition;
	
  G4double collimatorBoxYSize;
  G4double collimatorBoxZSize;

  G4double collimatorBox_XPosition;
  G4double collimatorBox_YPosition;
  G4double collimatorBox_ZPosition;

  G4Box *solidCollimator;
  G4LogicalVolume *logicCollimator;
  G4VPhysicalVolume *physicCollimator;

  G4Tubs *solidCollimatorHole;
  G4LogicalVolume *logicCollimatorHole;
  G4VPhysicalVolume *physicCollimatorHole;
  
  G4Material *collimatorHoleMaterial;
  G4Material *collimatorMaterial;
  // Variables definition of the final collimator
  
  G4double defaultInnerRadiusFinalCollimator;
  G4double innerRadiusFinalCollimator;
	
  G4double defaultFinalCollimatorThickness;
  G4double FinalCollimatorThickness;
	
  G4double defaultOuterRadiusFinalCollimator;
  G4double outerRadiusFinalCollimator; 
	
  G4double defaultStartAngleFinalCollimator;
  G4double startAngleFinalCollimator;
	
  G4double defaultSpanningAngleFinalCollimator;
  G4double spanningAngleFinalCollimator;
	
  G4double defaultFinalCollimatorXPosition;
  G4double FinalcollimatorXPosition;
	
  G4double defaultFinalCollimatorYPosition;
  G4double FinalcollimatorYPosition;
	
  G4double defaultFinalCollimatorZPosition;
  G4double FinalcollimatorZPosition;
   
  G4double collimatorFinalBoxXSize;
  G4double collimatorFinalBoxYSize;
  G4double collimatorFinalBoxZSize;

  G4double collimatorFinalBox_XPosition;
  G4double collimatorFinalBox_YPosition;
  G4double collimatorFinalBox_ZPosition;

  G4Box *solidFinalCollimator;
  G4LogicalVolume *logicFinalCollimator;
  G4VPhysicalVolume *physicFinalCollimator;

  G4Tubs *solidFinalCollimatorHole;
  G4LogicalVolume *logicFinalCollimatorHole;
  G4VPhysicalVolume *physicFinalCollimatorHole;

  G4Material *FinalcollimatorMaterial;
  G4Material *FinalcollimatorHoleMaterial;

  G4Material *WindowMaterial;
  
  G4Material *PipeMaterial;
 

  // Variables definition of the magnetic component
  G4Material *MotherMaterial;
  
  G4Material *externalMagnet_1Material,  *externalMagnet_2Material, *externalMagnet_3Material, *externalMagnet_4Material, *externalSlitMaterial, *internalSlitMaterial;

  G4Box *solidExternalMagnet_1;
  G4LogicalVolume *logicExternalMagnet_1;
  G4VPhysicalVolume *physicExternalMagnet_1; 
  G4VPhysicalVolume *physicExternalMagnet_1Down;
 ///////000000ooooooo0000000//////////
  G4Box *solidMagnet_1;
  G4LogicalVolume *logicMagnet_1;
  G4VPhysicalVolume *physicMagnet_1Right; 
  G4VPhysicalVolume *physicMagnet_1Left;  
	
  G4Box *solidExternalMagnet_2;
  G4LogicalVolume *logicExternalMagnet_2;
  G4VPhysicalVolume *physicExternalMagnet_2; 
  G4VPhysicalVolume *physicExternalMagnet_2Down;

  G4Box *solidMagnet_2;
  G4LogicalVolume *logicMagnet_2;
  G4VPhysicalVolume *physicMagnet_2Right; 
  G4VPhysicalVolume *physicMagnet_2Left; 

  G4Box *solidExternalMagnet_3;
  G4LogicalVolume *logicExternalMagnet_3;
  G4VPhysicalVolume *physicExternalMagnet_3; 
  G4VPhysicalVolume *physicExternalMagnet_3Down; 

  G4Box *solidMagnet_3;
  G4LogicalVolume *logicMagnet_3;
  G4VPhysicalVolume *physicMagnet_3Right; 
  G4VPhysicalVolume *physicMagnet_3Left; 

  G4Box *solidExternalMagnet_4;
  G4LogicalVolume *logicExternalMagnet_4;
  G4VPhysicalVolume *physicExternalMagnet_4; 
  G4VPhysicalVolume *physicExternalMagnet_4Down;
	
  G4Box *solidMagnet_4;
  G4LogicalVolume *logicMagnet_4;
  G4VPhysicalVolume *physicMagnet_4Right; 
  G4VPhysicalVolume *physicMagnet_4Left; 

  G4Box *solidExternalSlit;
  G4LogicalVolume *logicExternalSlit;
  G4VPhysicalVolume *physicExternalSlit; 
	
  G4Box *solidInternalSlit;
  G4LogicalVolume *logicInternalSlit;
  G4VPhysicalVolume *physicInternalSlit;


  G4double externalMagnet_1XSize;
  G4double externalMagnet_1YSize;
  G4double externalMagnet_1ZSize;
	
  G4double externalMagnet_2XSize;
  G4double externalMagnet_2YSize;
  G4double externalMagnet_2ZSize;
	
  G4double externalMagnet_3XSize;
  G4double externalMagnet_3YSize;
  G4double externalMagnet_3ZSize;
	
  G4double externalMagnet_4XSize;
  G4double externalMagnet_4YSize;
  G4double externalMagnet_4ZSize;
	
  G4double externalMagnet_1XPosition;
  G4double externalMagnet_1YPosition;
  G4double externalMagnet_1ZPosition;
	
  G4double externalMagnet_2XPosition;
  G4double externalMagnet_2YPosition;
  G4double externalMagnet_2ZPosition;
	
  G4double externalMagnet_3XPosition;
  G4double externalMagnet_3YPosition;
  G4double externalMagnet_3ZPosition;
	
  G4double externalMagnet_4XPosition;
  G4double externalMagnet_4YPosition;
  G4double externalMagnet_4ZPosition;

  G4double externalSlitXPosition;
  G4double externalSlitYPosition;
  G4double externalSlitZPosition;
	
  G4double externalSlitXSize;
  G4double externalSlitYSize;
  G4double externalSlitZSize; 
  
  G4Tubs *solidExitPipe; 
  G4LogicalVolume *logicExitPipe; 
  G4VPhysicalVolume *physicExitPipe; 
  
  G4Tubs *solidExitWindow; 
  G4LogicalVolume *logicExitWindow; 
  G4VPhysicalVolume *physicExitWindow; 
  
  G4Tubs *solidExithole; 
  G4LogicalVolume *logicExithole; 
  G4VPhysicalVolume *physicExithole; 
  
  G4Tubs *solidEntrancePipe; 
  G4LogicalVolume *logicEntrancePipe; 
  G4VPhysicalVolume *physicEntrancePipe; 
    
  G4Tubs *solidEntrancehole; 
  G4LogicalVolume *logicEntrancehole; 
  G4VPhysicalVolume *physicEntrancehole; 
  G4VPhysicalVolume *physicEntranceholeESSChamber; 

  G4double Magnet_1XPosition;
  G4double Magnet_1YPosition;
  G4double Magnet_1ZPosition;
	
  G4double Magnet_1XSize;
  G4double Magnet_1YSize;
  G4double Magnet_1ZSize;

  G4double Magnet_2XPosition;
  G4double Magnet_2YPosition;
  G4double Magnet_2ZPosition;
	
  G4double Magnet_2XSize;
  G4double Magnet_2YSize;
  G4double Magnet_2ZSize;

  G4double Magnet_3XPosition;
  G4double Magnet_3YPosition;
  G4double Magnet_3ZPosition;
	
  G4double Magnet_3XSize;
  G4double Magnet_3YSize;
  G4double Magnet_3ZSize;

  G4double Magnet_4XPosition;
  G4double Magnet_4YPosition;
  G4double Magnet_4ZPosition;
	
  G4double Magnet_4XSize;
  G4double Magnet_4YSize;
  G4double Magnet_4ZSize;

  G4double internalSlitXPosition;
  G4double internalSlitYPosition;
  G4double internalSlitZPosition;
	
  G4double internalSlitXSize;
  G4double internalSlitYSize;
  G4double internalSlitZSize;
    
////////////////////////////////////////// Faraday Cup /////////////////////////////////////////////


  G4Material *KaptonEntranceWindowMaterial;
  G4Material *MassRingMaterial;
  G4Material *GuardRingMaterial;
  G4Material *FaradayCupBottomMaterial;
  G4Material *CupMaterial;
    
    G4Box *virtualMag;
    G4LogicalVolume *logicVirtualMag;
    G4VPhysicalVolume *physicVirtualMag;	

    G4Box*Box;
    G4Tubs*Cylinder;
    G4LogicalVolume* logicBeveledCylinder;
    G4VPhysicalVolume* physicBeveledCylinder; 

    G4Tubs *KaptonEntranceWindow; 
    G4LogicalVolume *logicKaptonEntranceWindow;
    G4VPhysicalVolume *physicKaptonEntranceWindow;

    G4Tubs *MassRing;
    G4LogicalVolume *logicMassRing;
    G4VPhysicalVolume *physicMassRing;

    G4Tubs *VirtualWindow;
    G4LogicalVolume *logicVirtualWindow;
    G4VPhysicalVolume *physicVirtualWindow;  

    G4Tubs *GuardRing;
    G4LogicalVolume *logicGuardRing;
    G4VPhysicalVolume *physicGuardRing;
  
    G4Tubs *VirtualMiddle;
    G4LogicalVolume *logicVirtualMiddle;
    G4VPhysicalVolume *physicVirtualMiddle;
 
    G4Tubs *FaradayCupBottom; 
    G4LogicalVolume *logicFaradayCupBottom; 
    G4VPhysicalVolume *physicFaradayCupBottom;
 
    G4Tubs *VirtualBottom;
    G4LogicalVolume *logicVirtualBottom;
    G4VPhysicalVolume *physicVirtualBottom;
 
    G4Tubs *Cup; 
    G4LogicalVolume *logicCup; 
    G4VPhysicalVolume *physicCup;

    G4Tubs *VirtualOverBottom;
    G4LogicalVolume *logicVirtualOverBottom;
    G4VPhysicalVolume *physicVirtualOverBottom;

    G4Tubs *VirtualLateral;
    G4LogicalVolume *logicVirtualLateral;
    G4VPhysicalVolume *physicVirtualLateral;

  // Colours
  G4VisAttributes* blue;
  G4VisAttributes* gray;
  G4VisAttributes* white;
  G4VisAttributes* red;
  G4VisAttributes* yellow;
  G4VisAttributes* green;
  G4VisAttributes* darkGreen;
  G4VisAttributes* darkOrange3;
  G4VisAttributes* skyBlue;	
  G4VisAttributes* black;
};
#endif
