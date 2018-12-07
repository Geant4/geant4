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

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4NistManager.hh"
#include "G4NistElementBuilder.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "HadrontherapyDetectorConstruction.hh"
#include "LaserDrivenBeamLine.hh"
#include "LaserDrivenBeamLineMessenger.hh"
//
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4Material.hh"
//
#include "G4FieldManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
//#include "G4TransportationManager.hh"
#include "G4EqMagElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4PropagatorInField.hh"
#include "G4VisCommandsViewer.hh"
#include "G4UImanager.hh"
#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4SubtractionSolid.hh"

//
#include "G4UniformElectricField.hh"
#include "G4ElectricField.hh"
#include "HadrontherapyElectricTabulatedField3D.hh"

#include "HadrontherapyMagneticField3D.hh"
//
//G4bool LaserDrivenBeamLine::doCalculation = false;
/////////////////////////////////////////////////////////////////////////////
LaserDrivenBeamLine::LaserDrivenBeamLine():
hadrontherapydetectorconstruction(0), physicTreatmentRoom(0),
PFirstTriplet(0),PSecondTriplet(0),PThirdTriplet(0),PFourthTriplet(0), physicFirstQuad(0),physicSecondQuad(0),physicThirdQuad(0),physicFourthQuad(0),
solidExternalChamber(0),logicExternalChamber(0),physicExternalChamber(0),
solidInternalChamber(0),logicInternalChamber(0),physicInternalChamber(0),
solidCollimator(0),logicCollimator(0),physicCollimator(0),
solidCollimatorHole(0),logicCollimatorHole(0),physicCollimatorHole(0),
solidFinalCollimator(0), logicFinalCollimator(0),physicFinalCollimator(0),
solidFinalCollimatorHole(0),logicFinalCollimatorHole(0),physicFinalCollimatorHole(0),
solidExternalMagnet_1(0),logicExternalMagnet_1(0),physicExternalMagnet_1(0), physicExternalMagnet_1Down(0),
solidMagnet_1(0),logicMagnet_1(0),physicMagnet_1Right(0),physicMagnet_1Left(0), solidExternalMagnet_2(0),logicExternalMagnet_2(0),
physicExternalMagnet_2(0),physicExternalMagnet_2Down(0),solidMagnet_2(0),logicMagnet_2(0),physicMagnet_2Right(0),physicMagnet_2Left(0), solidExternalMagnet_3(0),logicExternalMagnet_3(0),physicExternalMagnet_3(0),physicExternalMagnet_3Down(0),
solidMagnet_3(0),logicMagnet_3(0),physicMagnet_3Right(0),physicMagnet_3Left(0),
solidExternalMagnet_4(0),logicExternalMagnet_4(0),physicExternalMagnet_4(0),physicExternalMagnet_4Down(0),
solidMagnet_4(0),logicMagnet_4(0),physicMagnet_4Right(0),physicMagnet_4Left(0),
solidExternalSlit(0), logicExternalSlit(0), physicExternalSlit(0),
solidInternalSlit(0),logicInternalSlit(0),physicInternalSlit(0),
physicExitPipe(0),physicExitWindow(0),physicExithole(0),physicEntrancePipe(0),physicEntrancehole(0)
{
    laserDrivenMessenger = new LaserDrivenBeamLineMessenger(this);
    
    //***************************** PW ***************************************
    
    static G4String ROGeometryName = "DetectorROGeometry";
    RO = new HadrontherapyDetectorROGeometry(ROGeometryName);
    
    G4cout << "Going to register Parallel world...";
    RegisterParallelWorld(RO);
    G4cout << "... done" << G4endl;
    //***************************** PW ***************************************
}

/////////////////////////////////////////////////////////////////////////////
LaserDrivenBeamLine::~LaserDrivenBeamLine()
{
    //delete laserDrivenMessenger;
    delete hadrontherapydetectorconstruction;
}

/////////////////////////////////////////////////////////////////////////////
G4VPhysicalVolume* LaserDrivenBeamLine::Construct()
{
    // Sets default geometry and materials
    SetDefaultDimensions();
    
    // Construct the energyselector (magnetic part and slit) and detector plane
    ConstructLaserDrivenBeamLine();
    
    //***************************** PW ***************************************
    if (!hadrontherapydetectorconstruction)
        
        //***************************** PW ***************************************
        
        // HadrontherapyDetectorConstruction builds ONLY the phantom and the detector with its associated ROGeometry
        hadrontherapydetectorconstruction = new HadrontherapyDetectorConstruction(physicTreatmentRoom);
    G4cout<<"HadrontherapyDetectorConstruction"<<G4endl;
    //***************************** PW ***************************************
    
    hadrontherapydetectorconstruction->InitializeDetectorROGeometry(RO,hadrontherapydetectorconstruction->GetDetectorToWorldPosition());
    
    //***************************** PW ***************************************
    return physicTreatmentRoom;
}
/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::SetDefaultDimensions()
{
    ///////////////////////////////////////////////////////////////////////
    // Definition of the colour sets
    white = new G4VisAttributes( G4Colour(1.,1.,1., 0.2));
    white -> SetVisibility(true);
    white -> SetForceSolid(true);
    white -> SetForceWireframe(true);
    
    blue = new G4VisAttributes(G4Colour(0. ,0. ,1.));
    blue -> SetVisibility(true);
    //blue -> SetForceSolid(true);
    
    gray = new G4VisAttributes( G4Colour(0.5, 0.5, 0.5, 0.5 ));
    gray-> SetVisibility(true);
    gray-> SetForceSolid(true);
    
    red = new G4VisAttributes(G4Colour(1. ,0. ,0., 0.2));
    red-> SetVisibility(true);
    red-> SetForceSolid(true);
    //red -> SetForceWireframe(true);
    
    yellow = new G4VisAttributes(G4Colour(1., 1., 0., 0.2));
    yellow-> SetVisibility(true);
    yellow-> SetForceSolid(true);
    
    green = new G4VisAttributes( G4Colour(25/255. , 255/255. ,  25/255., 0.4));
    green -> SetVisibility(true);
    green -> SetForceWireframe(true);
    green -> SetForceSolid(true);
    
    black = new G4VisAttributes( G4Colour(255/255. , 255/255.,  255/255.));
    black -> SetVisibility(true);
    black -> SetForceSolid(true);
    
    darkGreen = new G4VisAttributes( G4Colour(0/255. , 100/255. ,  0/255.));
    darkGreen -> SetVisibility(true);
    darkGreen -> SetForceSolid(true);
    
    darkOrange3 = new G4VisAttributes( G4Colour(205/255. , 102/255. ,  000/255., 0.7));
    darkOrange3 -> SetVisibility(true);
    darkOrange3 -> SetForceSolid(true);
    
    skyBlue = new G4VisAttributes( G4Colour(135/255. , 206/255. ,  235/255., 0.1));
    skyBlue -> SetVisibility(true);
    skyBlue -> SetForceSolid(true);
    
    // DEFAULT DIMENSIONS AND POSITIONS ARE PROVIDED HERE.
    /////////////////////// Exit Window ///////////////////////////////////////////////
    G4double defaultInnerRadiusExitWindow=0. *mm;
    InnerRadiusExitWindow=defaultInnerRadiusExitWindow;
    
    G4double defaultExternalRadiusExitWindow=55*mm;
    ExternalRadiusExitWindow=defaultExternalRadiusExitWindow;
    
    G4double defaultExitWindowThickness=25 *um;
    ExitWindowThickness=defaultExitWindowThickness;
    
    G4double defaultExitWindowXPosition=-ExitWindowThickness/2.;
    ExitWindowXPosition=defaultExitWindowXPosition;
    
    G4double defaultExitWindowYPosition=0.;
    ExitWindowYPosition=defaultExitWindowYPosition;
    
    G4double defaultExitWindowZPosition=0.0*mm;
    ExitWindowZPosition=defaultExitWindowZPosition;
    
    G4double defaultStartAngleExitWindow = 0.0 *deg;
    startAngleExitWindow = defaultStartAngleExitWindow;
    
    G4double defaultSpanningAngleExitWindow = 360.*deg;
    spanningAngleExitWindow = defaultSpanningAngleExitWindow;
    ////////////////////////////// Exit pipe ////////////////////////////////
    G4double defaultExitPipeheight=105. *mm;
    ExitPipeheight=defaultExitPipeheight;
    
    G4double defaultInnerRadiusExitPipe=50. *mm;
    InnerRadiusExitPipe=defaultInnerRadiusExitPipe;
    
    G4double defaultExternalRadiusExitPipe=55 *mm;
    ExternalRadiusExitPipe=defaultExternalRadiusExitPipe;
    
    G4double defaultExitPipeXPosition=-ExitPipeheight/2-ExitWindowThickness;
    ExitPipeXPosition=defaultExitPipeXPosition;
    
    G4double defaultExitPipeYPosition=0;
    ExitPipeYPosition=defaultExitPipeYPosition;
    
    G4double defaultExitPipeZPosition=0.0*mm;
    ExitPipeZPosition=defaultExitPipeZPosition;
    
    G4double defaultStartAngleExitPipe = 0.0 *deg;
    startAngleExitPipe = defaultStartAngleExitPipe;
    
    G4double defaultSpanningAngleExitPipe = 360.*deg;
    spanningAngleExitPipe = defaultSpanningAngleExitPipe;
    //////////////////////////////////////////////// Vacuum chamber //////////////////////////////
    G4double defaultExternalChamberXSize = 79.6*cm;
    externalChamberXSize = defaultExternalChamberXSize;
    
    G4double defaultExternalChamberYSize = 50. *cm;
    externalChamberYSize = defaultExternalChamberYSize;
    
    G4double defaultExternalChamberZSize = 50. *cm;
    externalChamberZSize = defaultExternalChamberZSize;
    
    G4double defaultExternalChamberXPosition = -(externalChamberXSize/2.+ExitPipeheight/2.)+ ExitPipeXPosition;
    externalChamberXPosition = defaultExternalChamberXPosition;
    
    G4double defaultExternalChamberYPosition = 0.0 *mm;
    externalChamberYPosition = defaultExternalChamberYPosition;
    
    G4double defaultExternalChamberZPosition = 0.0 *mm;
    externalChamberZPosition = defaultExternalChamberZPosition;
    
    // Defaults of the internal chamber dimensions
    // The position of its center is in the center
    // of the internal chamber while the dimension are
    // authomatically calculated respect to the external chamber ones
    G4double defaultVaccumChamberWallThickness=5 *mm;
    VaccumChamberWallThickness=defaultVaccumChamberWallThickness;
    
    G4double defaultInternalChamberXSize =externalChamberXSize - 2*VaccumChamberWallThickness;
    internalChamberXSize = defaultInternalChamberXSize;
    
    G4double defaultInternalChamberYSize =externalChamberYSize - 2*VaccumChamberWallThickness;
    internalChamberYSize = defaultInternalChamberYSize;
    
    G4double defaultInternalChamberZSize = externalChamberZSize - 2*VaccumChamberWallThickness;
    internalChamberZSize = defaultInternalChamberZSize;
    /////////////////////// Exit hole in vessel ///////////////////////////////////////////////
    G4double defaultInnerRadiusExithole=0.*mm;
    InnerRadiusExithole=defaultInnerRadiusExithole;
    
    G4double defaultExternalRadiusExithole=50.*mm;
    ExternalRadiusExithole=defaultExternalRadiusExithole;
    
    G4double defaultExitholeThickness=VaccumChamberWallThickness;
    ExitholeThickness=defaultExitholeThickness;
    
    G4double defaultExitholeXPosition=(externalChamberXSize/2.-ExitholeThickness/2.);
    ExitholeXPosition=defaultExitholeXPosition;
    
    G4double defaultExitholeYPosition=0.;
    ExitholeYPosition=defaultExitholeYPosition;
    
    G4double defaultExitholeZPosition=0.*mm;
    ExitholeZPosition=defaultExitholeZPosition;
    
    G4double defaultStartAngleExithole = 0.0 *deg;
    startAngleExithole= defaultStartAngleExithole;
    
    G4double defaultSpanningAngleExithole = 360.*deg;
    spanningAngleExithole = defaultSpanningAngleExithole;
    /////////////////////////////////Final collimator //////////////////////////////
    // The Final Collimator is located after  the 4th magnet
    G4double defaultExitholeToFinalCollimator=70 *mm;
    ExitholeToFinalCollimator=defaultExitholeToFinalCollimator;
    
    defaultInnerRadiusFinalCollimator = 0.0 *mm;
    innerRadiusFinalCollimator = defaultInnerRadiusFinalCollimator;
    
    defaultOuterRadiusFinalCollimator = 2.50 *mm;
    outerRadiusFinalCollimator = defaultOuterRadiusFinalCollimator;
    
    defaultFinalCollimatorThickness = 3.0 *mm;
    FinalCollimatorThickness = defaultFinalCollimatorThickness;
    
    defaultStartAngleFinalCollimator = 0.0 *deg;
    startAngleFinalCollimator = defaultStartAngleFinalCollimator;
    
    defaultSpanningAngleFinalCollimator = 360.*deg;
    spanningAngleFinalCollimator = defaultSpanningAngleFinalCollimator;
    
    defaultFinalCollimatorXPosition = internalChamberXSize/2.-ExitholeToFinalCollimator-FinalCollimatorThickness/2.;
    collimatorFinalBox_XPosition=defaultFinalCollimatorXPosition;
    FinalcollimatorXPosition = 0.0*mm; //HOLE IN THE FINAL COLLIMATOR
    
    defaultFinalCollimatorYPosition = 0.0*mm;
    collimatorFinalBox_YPosition=defaultFinalCollimatorYPosition;
    FinalcollimatorYPosition = defaultFinalCollimatorYPosition;
    
    defaultFinalCollimatorZPosition = 0.0*mm;
    collimatorFinalBox_ZPosition=0.0*mm;
    FinalcollimatorZPosition =defaultFinalCollimatorZPosition;
    
    defaultThicknessCollimator =3.0 *mm;
    collimatorFinalBoxXSize=defaultFinalCollimatorThickness;
    collimatorFinalBoxYSize=82.0*mm;
    collimatorFinalBoxZSize=210.0*mm;
    //////////////////ooooooooooOOOOOOOO000000000000OOOOOOOOOOOOooooooooooo/////////////////
    
    //Magnet characteristics
    G4double defaultExternalMagnet_XSize = 88.0*mm;
    G4double defaultExternalMagnet_YSizeTotal=87.*mm;
    G4double defaultInternalMagnet_YSize = 10. *mm;
    G4double defaultExternalMagnet_YSize =(defaultExternalMagnet_YSizeTotal-defaultInternalMagnet_YSize)/2.;
    G4double defaultExternalMagnet_ZSize = 104 *mm;
    
    G4double defaultExternalMagnet_YPosition =defaultInternalMagnet_YSize/2.+defaultExternalMagnet_YSize/2.;
    G4double defaultExternalMagnet_ZPosition = 0.0 *mm;
    
    G4double defaultMagnet_XSize=defaultExternalMagnet_XSize;
    G4double defaultMagnet_YSize=defaultExternalMagnet_YSizeTotal;
    G4double defaultMagnet_ZSize=19*mm;
    
    // Defaults of the external part of the magnet 4:
    G4double defaultFinalCollimatorToMagnet4=25.*mm;
    FinalCollimatorToMagnet4=defaultFinalCollimatorToMagnet4;
    
    externalMagnet_4XSize = defaultExternalMagnet_XSize;
    externalMagnet_4YSize = defaultExternalMagnet_YSize;
    externalMagnet_4ZSize = defaultExternalMagnet_ZSize;
    
    Magnet_4XSize=defaultMagnet_XSize;
    Magnet_4YSize=defaultMagnet_YSize;
    Magnet_4ZSize=defaultMagnet_ZSize;
    
    G4double defaultExternalMagnet_4XPosition = -(FinalCollimatorThickness/2.+FinalCollimatorToMagnet4+defaultExternalMagnet_XSize/2.)+ collimatorFinalBox_XPosition;
    externalMagnet_4XPosition = defaultExternalMagnet_4XPosition;
    
    externalMagnet_4YPosition = defaultExternalMagnet_YPosition;
    externalMagnet_4ZPosition = defaultExternalMagnet_ZPosition;
    
    Magnet_4XPosition=externalMagnet_4XPosition;
    Magnet_4YPosition=0.0*mm;
    Magnet_4ZPosition=(defaultExternalMagnet_ZSize+defaultMagnet_ZSize)/2.;
    //////////////////ooooooooooOOOOOOOO000000000000OOOOOOOOOOOOooooooooooo/////////////////
    // Defaults of the external part of the magnet 3:
    externalMagnet_3XSize = defaultExternalMagnet_XSize;
    externalMagnet_3YSize = defaultExternalMagnet_YSize;
    externalMagnet_3ZSize = defaultExternalMagnet_ZSize;
    
    Magnet_3XSize=defaultMagnet_XSize;
    Magnet_3YSize=defaultMagnet_YSize;
    Magnet_3ZSize=defaultMagnet_ZSize;
    
    G4double defaultMagnet4ToMagnet3=65.*mm; //85.*mm ANTONELLA
    Magnet4ToMagnet3=defaultMagnet4ToMagnet3;
    
    G4double defaultExternalMagnet_3XPosition =-(Magnet4ToMagnet3+defaultExternalMagnet_XSize/2.+defaultExternalMagnet_XSize/2.)+externalMagnet_4XPosition;
    externalMagnet_3XPosition = defaultExternalMagnet_3XPosition;
    
    externalMagnet_3YPosition =defaultExternalMagnet_YPosition;
    externalMagnet_3ZPosition = defaultExternalMagnet_ZPosition;
    
    
    Magnet_3XPosition=externalMagnet_3XPosition;
    Magnet_3YPosition=0.0*mm;
    Magnet_3ZPosition=(defaultExternalMagnet_ZSize+defaultMagnet_ZSize)/2.;
    //////////////////ooooooooooOOOOOOOO000000000000OOOOOOOOOOOOooooooooooo/////////////////
    // Defaults of the external part of the magnet 2:
    externalMagnet_2XSize = defaultExternalMagnet_XSize;
    externalMagnet_2YSize = defaultExternalMagnet_YSize;
    externalMagnet_2ZSize = defaultExternalMagnet_ZSize;
    
    Magnet_2XSize=defaultMagnet_XSize;
    Magnet_2YSize=defaultMagnet_YSize;
    Magnet_2ZSize=defaultMagnet_ZSize;
    
    G4double defaultMagnet3ToMagnet2=10 *mm;
    Magnet3ToMagnet2=defaultMagnet3ToMagnet2;
    
    G4double defaultExternalMagnet_2XPosition =-(Magnet3ToMagnet2+defaultExternalMagnet_XSize/2.+defaultExternalMagnet_XSize/2.)+externalMagnet_3XPosition;
    externalMagnet_2XPosition = defaultExternalMagnet_2XPosition;
    
    externalMagnet_2YPosition = defaultExternalMagnet_YPosition;
    externalMagnet_2ZPosition = defaultExternalMagnet_ZPosition;
    
    Magnet_2XPosition=externalMagnet_2XPosition;
    Magnet_2YPosition=0.0*mm;
    Magnet_2ZPosition=(defaultExternalMagnet_ZSize+defaultMagnet_ZSize)/2.;
    //////////////////ooooooooooOOOOOOOO000000000000OOOOOOOOOOOOooooooooooo/////////////////
    // Defaults of the external part of the magnet 1:
    externalMagnet_1XSize=defaultExternalMagnet_XSize;
    externalMagnet_1YSize = defaultExternalMagnet_YSize;
    externalMagnet_1ZSize = defaultExternalMagnet_ZSize;
    
    Magnet_1XSize=defaultMagnet_XSize;
    Magnet_1YSize=defaultMagnet_YSize;
    Magnet_1ZSize=defaultMagnet_ZSize;
    
    G4double defaultMagnet2ToMagnet1=85 *mm;
    Magnet2ToMagnet1=defaultMagnet2ToMagnet1;
    
    G4double defaultExternalMagnet_1XPosition = -(Magnet2ToMagnet1+defaultExternalMagnet_XSize/2.+defaultExternalMagnet_XSize/2.)+externalMagnet_2XPosition;
    externalMagnet_1XPosition = defaultExternalMagnet_1XPosition;
    
    externalMagnet_1YPosition = defaultExternalMagnet_YPosition;
    externalMagnet_1ZPosition = defaultExternalMagnet_ZPosition;
    
    Magnet_1XPosition=defaultExternalMagnet_1XPosition;
    Magnet_1YPosition=0.0*mm;
    Magnet_1ZPosition=(defaultExternalMagnet_ZSize+defaultMagnet_ZSize)/2.;
    
    // Defaults of the external part of the Slit
    G4double defaultExternalSlitXSize = 8.0 *mm;
    externalSlitXSize = defaultExternalSlitXSize;
    
    G4double defaultExternalSlitYSize = 82. *mm;
    externalSlitYSize = defaultExternalSlitYSize;
    
    G4double defaultExternalSlitZSize = 210. *mm;
    externalSlitZSize = defaultExternalSlitZSize;
    
    G4double defaultExternalSlitXPosition = -(Magnet3ToMagnet2/2.+defaultExternalMagnet_XSize/2.)+externalMagnet_3XPosition;
    externalSlitXPosition = defaultExternalSlitXPosition;
    
    G4double defaultExternalSlitYPosition = 0.0 *mm;
    externalSlitYPosition = defaultExternalSlitYPosition;
    
    G4double defaultExternalSlitZPosition = 0.0 *mm;
    externalSlitZPosition = defaultExternalSlitZPosition;
    
    // Defaults of the internal part of the Slit:
    internalSlitXSize = defaultExternalSlitXSize;
    
    G4double defaultInternalSlitYSize = 3 *mm;
    internalSlitYSize = defaultInternalSlitYSize;
    
    G4double defaultInternalSlitZSize = 3 *mm;
    internalSlitZSize = defaultInternalSlitZSize;
    
    G4double defaultInternalSlitXPosition = 0.0 *mm;
    internalSlitXPosition = defaultInternalSlitXPosition;
    
    G4double defaultInternalSlitYPosition = 0.0 *mm;
    internalSlitYPosition = defaultInternalSlitYPosition;
    
    G4double defaultInternalSlitZPosition = 40.0 *mm;
    internalSlitZPosition = defaultInternalSlitZPosition;
    
    // Defaults of the particle collimator (First collimator).
    // The Collimator should be located before the 1st magnet
    //
    defaultInnerRadiusCollimator = 0.0 *mm;
    innerRadiusCollimator = defaultInnerRadiusCollimator;
    
    defaultOuterRadiusCollimator = 2.5 *mm;
    outerRadiusCollimator = defaultOuterRadiusCollimator;
    
    thicknessCollimator = defaultThicknessCollimator;
    
    defaultStartAngleCollimator = 0.0 *deg;
    startAngleCollimator = defaultStartAngleCollimator;
    
    defaultSpanningAngleCollimator = 360.*deg;
    spanningAngleCollimator = defaultSpanningAngleCollimator;
    
    G4double defultMagnet1ToFirstCollimator=25.*mm;
    Magnet1ToFirstCollimator=defultMagnet1ToFirstCollimator;
    
    defaultCollimatorXPosition = -(thicknessCollimator/2.+Magnet1ToFirstCollimator+defaultExternalMagnet_XSize/2.)+externalMagnet_1XPosition;
    collimatorBox_XPosition=defaultCollimatorXPosition;
    collimatorXPosition = 0.0*mm;
    
    defaultCollimatorYPosition = 0.0*mm;
    collimatorBox_YPosition=defaultCollimatorYPosition;
    collimatorYPosition = 0.0*mm;
    
    defaultCollimatorZPosition = 0.0*mm;
    collimatorBox_ZPosition=defaultCollimatorZPosition;
    collimatorZPosition = 0.*mm;
    
    collimatorBoxYSize=82.0* mm;
    collimatorBoxZSize=210.0* mm;
    
    //////////////////// Entrance Hole //////////////////////////////////
    G4double defaultInnerRadiusEntrancehole=0. *mm;
    InnerRadiusEntrancehole=defaultInnerRadiusEntrancehole;
    
    G4double defaultExternalRadiusEntrancehole=50.*mm;
    ExternalRadiusEntrancehole=defaultExternalRadiusEntrancehole;
    
    G4double defaultEntranceholeThickness=VaccumChamberWallThickness;
    EntranceholeThickness=defaultEntranceholeThickness;
    
    G4double defaultEntranceholeXPosition=-(externalChamberXSize/2.-EntranceholeThickness/2.);
    EntranceholeXPosition=defaultEntranceholeXPosition;
    
    G4double defaultEntranceholeQuadXPosition=+(externalChamberXSize/2.-EntranceholeThickness/2.);
    EntranceholeQuadXPosition=defaultEntranceholeQuadXPosition;
    
    G4double defaultEntranceholeYPosition=0.;
    EntranceholeYPosition=defaultEntranceholeYPosition;
    
    G4double defaultEntranceholeZPosition=0.0*mm;
    EntranceholeZPosition=defaultEntranceholeZPosition;
    
    G4double defaultStartAngleEntrancehole= 0.0 *deg;
    startAngleEntrancehole= defaultStartAngleEntrancehole;
    
    G4double defaultSpanningAngleEntrancehole= 360.*deg;
    spanningAngleEntrancehole=defaultSpanningAngleEntrancehole;
    
    ///////////////// Entrance Pipe/////////////////////////////////////////////
    
    G4double defaultEntrancePipeheight=105. *mm;
    EntrancePipeheight=defaultEntrancePipeheight;
    
    G4double defaultInnerRadiusEntrancePipe=50. *mm;
    InnerRadiusEntrancePipe=defaultInnerRadiusEntrancePipe;
    
    G4double defaultExternalRadiusEntrancePipe=55 *mm;
    ExternalRadiusEntrancePipe=defaultExternalRadiusEntrancePipe;
    
    G4double defaultEntrancePipeXPosition=-EntrancePipeheight/2-externalChamberXSize/2+externalChamberXPosition;
    EntrancePipeXPosition=defaultEntrancePipeXPosition;
    
    G4double defaultEntrancePipeYPosition=0;
    EntrancePipeYPosition=defaultEntrancePipeYPosition;
    
    G4double defaultEntrancePipeZPosition=0.0*mm;
    EntrancePipeZPosition=defaultEntrancePipeZPosition;
    
    G4double defaultStartAngleEntrancePipe= 0.0 *deg;
    startAngleEntrancePipe= defaultStartAngleEntrancePipe;
    
    G4double defaultSpanningAngleEntrancePipe= 360.*deg;
    spanningAngleEntrancePipe=defaultSpanningAngleEntrancePipe;
    
    /////////////////////////////////////Quadrupole//////////////////////////////////
    G4double defaultQuadChamberWallPosX=-(externalChamberXSize/2.)-EntrancePipeheight/2.+EntrancePipeXPosition;
    QuadChamberWallPosX=defaultQuadChamberWallPosX;
    G4double defaultQuadChamberWallPosY=0.0*cm;
    QuadChamberWallPosY=defaultQuadChamberWallPosY;
    G4double defaultQuadChamberWallPosZ=0.0*cm;
    QuadChamberWallPosZ=defaultQuadChamberWallPosZ;
    
    G4double defaultInnerRadiusQuad=10.0*mm;
    InnerRadiusQuad=defaultInnerRadiusQuad;
    
    G4double defaultInnerRadiusTriplet=0.0*mm;
    InnerRadiusTriplet=defaultInnerRadiusTriplet;
    
    G4double defaultExternalRadiusQuad=30.0*mm;
    ExternalRadiusQuad=defaultExternalRadiusQuad;
    
    G4double defaultFirstQuadThickness=80.0*mm;
    FirstQuadThickness=defaultFirstQuadThickness;
    G4double defaultSecondQuadThickness=40.0*mm;
    SecondQuadThickness=defaultSecondQuadThickness;
    G4double defaultThirdQuadThickness=40.0*mm;
    ThirdQuadThickness=defaultThirdQuadThickness;
    G4double defaultFourthQuadThickness=80.0*mm;
    FourthQuadThickness=defaultFourthQuadThickness;
    
    G4double defaultStartAngleQuad = 0.0 *deg;
    startAngleQuad = defaultStartAngleQuad;
    
    G4double defaultSpanningAngleQuad = 360.*deg;
    spanningAngleQuad = defaultSpanningAngleQuad;
    
    G4double distancefromQuadChamber=100.0*mm;
    G4double defaultFourthQuadXPosition= internalChamberXSize/2.-distancefromQuadChamber-FourthQuadThickness/2.;
    FourthQuadXPosition=defaultFourthQuadXPosition;
    FourthQXPosition=0.0*mm;
    
    G4double distanceFQuadTQuad=100.0*mm;
    G4double defaultThirdQuadXPosition=-ThirdQuadThickness/2.-distanceFQuadTQuad-FourthQuadThickness/2.+FourthQuadXPosition;
    ThirdQuadXPosition=defaultThirdQuadXPosition;
    ThirdQXPosition=0.0*mm;
    
    G4double distanceTQuadSQuad=100.0*mm;
    G4double defaultSecondQuadXPosition=-SecondQuadThickness/2.-distanceTQuadSQuad-ThirdQuadThickness/2.+ThirdQuadXPosition;
    SecondQuadXPosition=defaultSecondQuadXPosition;
    SecondQXPosition=0.0*mm;
    
    G4double distanceSQuadFQuad=100.0*mm;
    G4double defaultFirstQuadXPosition=-FirstQuadThickness/2.-distanceSQuadFQuad-SecondQuadThickness/2.+SecondQuadXPosition;
    FirstQuadXPosition=defaultFirstQuadXPosition;
    FirstQXPosition=0.0*mm;
    
    G4double defaultQuadYPosition=0.0*mm;
    QuadYPosition=defaultQuadYPosition;
    QYPosition=defaultQuadYPosition;
    
    G4double defaultQuadTZPosition= 0.*mm;
    QuadZPosition=defaultQuadTZPosition;
    G4double defaultQuadZPosition=0.0*mm;
    QZPosition=defaultQuadZPosition;
    
    // DEFAULT DEFINITION OF THE MATERIALS
    // All elements and compound definition follows the NIST database
    
    //ELEMENTS
    G4bool isotopes = false;
    G4Element* zincNist = G4NistManager::Instance()->FindOrBuildElement("Zn");
    G4Element* copperNist = G4NistManager::Instance()->FindOrBuildElement("Cu");
    
    //COMPOUNDS
    G4Material* ironNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Fe", isotopes);
    G4Material* aluminiumNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_Al");
    G4Material* kaptonNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_KAPTON", isotopes);
    //G4Material* waterNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_WATER", isotopes);
    G4Material* stainless_steelNist = G4NistManager::Instance()->FindOrBuildMaterial("G4_STAINLESS-STEEL", isotopes);
    
    // Elements and compunds not pre-defined in Geant4
    G4double d; // Density
    G4int nComponents;// Number of components
    G4double fractionmass; // Fraction in mass of an element in a material
    d = 8.40*g/cm3;
    nComponents = 2;
    G4Material* brass = new G4Material("Brass", d, nComponents);
    brass -> AddElement(zincNist, fractionmass = 30 *perCent);
    brass -> AddElement(copperNist, fractionmass = 70 *perCent);
    
    G4double atomicNumber = 1.;
    G4double massOfMole = 1.008*g/mole;
    d = 1.e-25*g/cm3;
    G4double temperature = 2.73*kelvin;
    G4double pressure = 3.e-18*pascal;
    G4Material* vacuum = new G4Material("interGalactic", atomicNumber,massOfMole, d, kStateGas,temperature, pressure);
    
    //***************************** PW ***************************************
    
    // DetectorROGeometry Material
    new G4Material("dummyMat", 1., 1.*g/mole, 1.*g/cm3);
    
    //***************************** PW ***************************************
    
    // MATERIAL ASSIGNMENT
    MotherMaterial=vacuum;
    QuadMaterial=ironNist;
    externalChamberMaterial = stainless_steelNist;
    internalChamberMaterial = vacuum;
    collimatorMaterial = aluminiumNist;
    collimatorHoleMaterial=vacuum;
    FinalcollimatorMaterial=aluminiumNist;
    FinalcollimatorHoleMaterial=vacuum;
    WindowMaterial=kaptonNist;
    PipeMaterial=stainless_steelNist;
    
    externalMagnet_1Material = ironNist;
    externalMagnet_2Material = ironNist;
    externalMagnet_3Material = ironNist;
    externalMagnet_4Material = ironNist;
    
    externalSlitMaterial = brass;
    internalSlitMaterial =vacuum;
    
    //FC Material
    
    KaptonEntranceWindowMaterial=kaptonNist;
    GuardRingMaterial=stainless_steelNist;
    FaradayCupBottomMaterial=aluminiumNist;
    CupMaterial=FaradayCupBottomMaterial;
    MassRingMaterial=GuardRingMaterial;
    
}

/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::ConstructLaserDrivenBeamLine()
{
    // -----------------------------
    // Treatment room - World volume
    //------------------------------
    
    const G4double worldX = 800.0 *cm;
    const G4double worldY = 400.0 *cm;
    const G4double worldZ = 400.0 *cm;
    
    solidTreatmentRoom = new G4Box("TreatmentRoom",
                                   worldX,
                                   worldY,
                                   worldZ);
    
    logicTreatmentRoom = new G4LogicalVolume(solidTreatmentRoom,
                                             MotherMaterial,
                                             "logicTreatmentRoom",
                                             0,
                                             0,
                                             0);
    
    physicTreatmentRoom = new G4PVPlacement(0,
                                            G4ThreeVector(),
                                            "physicalTreatmentRoom",
                                            logicTreatmentRoom,
                                            0,
                                            false,
                                            0);
    
    
    // The treatment room is invisible in the Visualisation
    logicTreatmentRoom -> SetVisAttributes (G4VisAttributes::GetInvisible());
    
    // The various components of the energyselector are constructed calling
    // the following methods
    
    // This method constructs the chamber where the energyselector is located
    EnergySelectorChamber();
    // This method construct the exit window
    ExitWindow();
    // This method construct the exit pipe
    ExitPipe();
    // This method construct the exit hole
    Exithole();
    
    // This method constructs a circular collimator of variable thickness and
    // aperture. It is placed befor the magnet to collimate particles caming from the
    // plasma;
    Collimator();
    
    // This method constructs the magnet 1 and its associated magnetic field
    Magnet_1();
    
    // This method constructs the magnet 2 and its associated magnetic field
    Magnet_2();
    
    // This method constructs the magnet 3 and its associated magnetic field
    Magnet_3();
    
    // This method constructs the magnet 4 and its associated magnetic field
    Magnet_4();
    
    // The selection slit is a square hole moveable inside a metallic plate
    Slit();
    
    FinalCollimator();
    
    // This method construct the quadrupoles
    Quadrupole();
    // This method construct the entrance hole
    Entrancehole();
    // This method construct the entrance pipe
    EntrancePipe();
    
    FaradayCup();
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void LaserDrivenBeamLine::ConstructSDandField()
{
    G4double minEps=1.0e-5;  //   Minimum & value for smallest steps
    G4double maxEps=1.0e-4;
    G4bool allLocal = true;
    // G4int nvar = 8;  For pure magnetic field, the number of integration variables is the default!
    
    //....oooOO0OOooo..........ENERGY SELECTOR SYSTEM FIELD..........oooOO0OOooo....
    if(logicInternalChamber){G4double xOffset =(internalChamberXSize/2.0)+externalSlitXPosition;
        PurgMagField = new HadrontherapyMagneticField3D("field/ESSMagneticField.TABLE", xOffset);
        pFieldMgr =new G4FieldManager();
        pFieldMgr -> SetDetectorField(PurgMagField);
        G4cout << "DeltaStep "<< pFieldMgr -> GetDeltaOneStep()/mm <<"mm" <<endl;
        pFieldMgr -> CreateChordFinder(PurgMagField);
        fEquation = new G4Mag_UsualEqRhs(PurgMagField);
        fstepper = new G4ClassicalRK4(fEquation);
        //////fstepper = new G4HelixImplicitEuler(fEquation);
        pIntgrDriver = new G4MagInt_Driver(1*mm,fstepper,fstepper-> GetNumberOfVariables());
        //the first parameter is the minimum step
        pChordFinder = new G4ChordFinder(pIntgrDriver);
        pFieldMgr->SetChordFinder(pChordFinder);
        pFieldMgr->SetMinimumEpsilonStep(minEps);
        pFieldMgr->SetMaximumEpsilonStep(maxEps);
        pFieldMgr->SetDeltaOneStep(0.5e-3*mm);//default value of DeltaChord is 0.25 mm
        logicInternalChamber -> SetFieldManager(pFieldMgr, allLocal);}
    //....oooOO0OOooo..........QUADS FIELDS..........oooOO0OOooo....
    //....oooOO0OOooo..........FOURTH QUAD FIELD..........oooOO0OOooo....
    if(LFourthTriplet){G4double xOffsetFQ =-(QuadChamberWallPosX+FourthQuadXPosition);
        PurgMagFieldQuadFourth = new HadrontherapyMagneticField3D("field/Quad80MagneticField.TABLE", xOffsetFQ);
        pFieldMgrQuadFourth =  new G4FieldManager();
        pFieldMgrQuadFourth -> SetDetectorField(PurgMagFieldQuadFourth);
        
        pFieldMgrQuadFourth -> CreateChordFinder(PurgMagFieldQuadFourth);
        fEquationQuadFourth = new G4Mag_UsualEqRhs(PurgMagFieldQuadFourth);
        fstepperQuadFourth = new G4ClassicalRK4(fEquationQuadFourth);
        pIntgrDriverQuadFourth = new G4MagInt_Driver(1*mm,fstepperQuadFourth,fstepperQuadFourth-> GetNumberOfVariables());
        //the first parameter is the minimum step
        pChordFinderQuadFourth = new G4ChordFinder(pIntgrDriverQuadFourth);
        pFieldMgrQuadFourth->SetChordFinder(pChordFinderQuadFourth);
        pFieldMgrQuadFourth->SetMinimumEpsilonStep(minEps);
        pFieldMgrQuadFourth->SetMaximumEpsilonStep(maxEps);
        pFieldMgrQuadFourth->SetDeltaOneStep(0.5e-3*mm);//default value of DeltaChord is 0.25 mm
        LFourthTriplet -> SetFieldManager(pFieldMgrQuadFourth, allLocal);}
    //....oooOO0OOooo..........THIRD QUAD FIELD..........oooOO0OOooo....
    if(LThirdTriplet){ G4double xOffsetTQ =-(QuadChamberWallPosX+ThirdQuadXPosition);
        PurgMagFieldQuadThird = new HadrontherapyMagneticField3D("field/Quad40MagneticField.TABLE", xOffsetTQ);
        pFieldMgrQuadThird =  new G4FieldManager();
        pFieldMgrQuadThird -> SetDetectorField(PurgMagFieldQuadThird);
        pFieldMgrQuadThird -> CreateChordFinder(PurgMagFieldQuadThird);
        fEquationQuadThird = new G4Mag_UsualEqRhs(PurgMagFieldQuadThird);
        fstepperQuadThird = new G4ClassicalRK4(fEquationQuadThird);
        pIntgrDriverQuadThird = new G4MagInt_Driver(1*mm,fstepperQuadThird,fstepperQuadThird-> GetNumberOfVariables());
        //the first parameter is the minimum step
        pChordFinderQuadThird = new G4ChordFinder(pIntgrDriverQuadThird);
        pFieldMgrQuadThird->SetChordFinder(pChordFinderQuadThird);
        pFieldMgrQuadThird->SetMinimumEpsilonStep(minEps);
        pFieldMgrQuadThird->SetMaximumEpsilonStep(maxEps);
        pFieldMgrQuadThird->SetDeltaOneStep(0.5e-3*mm);//default value of DeltaChord is 0.25 mm
        LThirdTriplet -> SetFieldManager(pFieldMgrQuadThird, allLocal);}
    //....oooOO0OOooo..........SECOND QUAD FIELD..........oooOO0OOooo....
    if(LSecondTriplet){G4double xOffsetSQ =-(QuadChamberWallPosX+SecondQuadXPosition);
        PurgMagFieldQuadSecond = new HadrontherapyMagneticField3D("field/Quad40MagneticField.TABLE", xOffsetSQ);
        pFieldMgrQuadSecond =  new G4FieldManager();
        pFieldMgrQuadSecond -> SetDetectorField(PurgMagFieldQuadSecond);
        pFieldMgrQuadSecond -> CreateChordFinder(PurgMagFieldQuadSecond);
        fEquationQuadSecond = new G4Mag_UsualEqRhs(PurgMagFieldQuadSecond);
        fstepperQuadSecond = new G4ClassicalRK4(fEquationQuadSecond);
        pIntgrDriverQuadSecond = new G4MagInt_Driver(1*mm,fstepperQuadSecond,fstepperQuadSecond-> GetNumberOfVariables());
        //the first parameter is the minimum step
        pChordFinderQuadSecond = new G4ChordFinder(pIntgrDriverQuadSecond);
        pFieldMgrQuadSecond->SetChordFinder(pChordFinderQuadSecond);
        pFieldMgrQuadSecond->SetMinimumEpsilonStep(minEps);
        pFieldMgrQuadSecond->SetMaximumEpsilonStep(maxEps);
        pFieldMgrQuadSecond->SetDeltaOneStep(0.5e-3*mm);//default value of DeltaChord is 0.25 mm
        LSecondTriplet -> SetFieldManager(pFieldMgrQuadSecond, allLocal);}
    //....oooOO0OOooo..........FIRST QUAD FIELD..........oooOO0OOooo....
    if(LFirstTriplet) {G4double xOffsetFirstQ =-(QuadChamberWallPosX+FirstQuadXPosition);
        PurgMagFieldQuadFirst = new HadrontherapyMagneticField3D("field/Quad80MagneticField.TABLE", xOffsetFirstQ);
        pFieldMgrQuadFirst =  new G4FieldManager();
        pFieldMgrQuadFirst -> SetDetectorField(PurgMagFieldQuadFirst);
        pFieldMgrQuadFirst -> CreateChordFinder(PurgMagFieldQuadFirst);
        fEquationQuadFirst = new G4Mag_UsualEqRhs(PurgMagFieldQuadFirst);
        fstepperQuadFirst = new G4ClassicalRK4(fEquationQuadFirst);
        pIntgrDriverQuadFirst = new G4MagInt_Driver(1*mm,fstepperQuadFirst,fstepperQuadFirst-> GetNumberOfVariables());
        //the first parameter is the minimum step
        pChordFinderQuadFirst = new G4ChordFinder(pIntgrDriverQuadFirst);
        pFieldMgrQuadFirst->SetChordFinder(pChordFinderQuadFirst);
        pFieldMgrQuadFirst->SetMinimumEpsilonStep(minEps);
        pFieldMgrQuadFirst->SetMaximumEpsilonStep(maxEps);
        pFieldMgrQuadFirst->SetDeltaOneStep(0.5e-3*mm);//default value of DeltaChord is 0.25 mm
        LFirstTriplet -> SetFieldManager(pFieldMgrQuadFirst, allLocal);}
    //....oooOO0OOooo..........FARADAY CUP FIELD..........oooOO0OOooo....
    if(logicVirtualMag) {G4double exOffset= -20*cm;
        G4double eyOffset= 0*cm;
        G4double ezOffset= 0*cm;
        G4FieldManager *pEFieldmanager = new G4FieldManager();
        G4ElectricField *ElectricField = new HadrontherapyElectricTabulatedField3D("field/ElectricFieldFC-600V.TABLE", exOffset, eyOffset, ezOffset);
        // UNIFORM FIELD
        // G4ElectroMagneticField* ElectricField = new G4UniformElectricField(G4ThreeVector(0.0, 10.0*volt/m, 0.0)); //G4UniformElectricField
        // The following is only for global field in the whole geometry
        //pEFieldmanager = G4TransportationManager::GetTransportationManager() -> GetFieldManager();
        
        const G4int nvarElectric=8;  // The Equation of motion for Electric (or combined Electric/Magnetic)
        // field requires 8 integration variables
        
        G4EqMagElectricField *fLocalEquation = new G4EqMagElectricField(ElectricField);
        G4MagIntegratorStepper* fLocalStepper = new G4ClassicalRK4(fLocalEquation, nvarElectric);
        G4MagInt_Driver  *pIntgrDriver_E = new G4MagInt_Driver(0.02*mm, fLocalStepper, fLocalStepper -> GetNumberOfVariables() );
        G4ChordFinder *fLocalChordFinder = new G4ChordFinder(pIntgrDriver_E);
        pEFieldmanager -> SetDetectorField(ElectricField);
        pEFieldmanager -> SetChordFinder(fLocalChordFinder);
        //G4double deltainter=0.0001*mm;
        //G4double missdist=0.1*mm;
        //pEFieldmanager->SetDeltaIntersection(deltainter);
        //fLocalChordFinder->SetDeltaChord(missdist);
        pEFieldmanager->SetMinimumEpsilonStep(minEps);
        pEFieldmanager->SetMaximumEpsilonStep(maxEps);
        pEFieldmanager->SetDeltaOneStep( 0.5e-3 * mm );
        //pEFieldmanager -> SetFieldChangesEnergy(true);
        logicVirtualMag -> SetFieldManager(pEFieldmanager, allLocal);}
    //....oooOO0OOooo....................oooOO0OOooo....
    G4cout<<" //....oooOO0OOooo.......... FIELDS HAVE BEEN IMPLEMENTED..........oooOO0OOooo...."<<G4endl;
    return;
}

/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::FaradayCup()
{
    /// FC sizes ///
    
    G4double InnerRadiusFC=25*mm;
    G4double OuterRadiusFC=45*mm;
    G4double MassRingThickness=5*mm;
    G4double GuardRingThickness=180*mm;
    G4double FaradayCupBottomThickness=120*mm;
    G4double CupThickness=10*cm;
    G4double KaptonEntranceWindowThickness=25*um;
    
    /// Virtual Volumes ///
    
    G4double VirtualWindowThickness=1.*um ;
    G4double VirtualMiddleThickness= 1.*um ;
    G4double VirtualBottomThickness= 1. *um ;
    G4double VirtualOverBottomThickness=1. *um ;
    G4double VirtualLateralLength=FaradayCupBottomThickness+CupThickness+VirtualBottomThickness;
    
    
    //// Position ////
    
    G4double virtualMagPosX=31*cm;
    G4double FC_XOffset=20*cm;
    G4double KaptonEntranceWindowPosX=-virtualMagPosX+KaptonEntranceWindowThickness/2+FC_XOffset;
    G4double MassRingPosX=KaptonEntranceWindowPosX+KaptonEntranceWindowThickness/2+MassRingThickness/2;
    G4double VirtualWindowPosX=MassRingPosX+MassRingThickness/2+VirtualWindowThickness/2;
    G4double GuardRingPosX=MassRingPosX+MassRingThickness/2+GuardRingThickness/2+2*mm;
    G4double VirtualMiddlePosX=GuardRingPosX+GuardRingThickness/2+VirtualMiddleThickness/2;
    G4double FaradayCupBottomPosX=GuardRingPosX+GuardRingThickness/2+FaradayCupBottomThickness/2+1*cm;
    G4double VirtualBottomPosX=FaradayCupBottomPosX+FaradayCupBottomThickness/2+VirtualBottomThickness/2;
    G4double CupPosX=VirtualBottomPosX+VirtualBottomThickness/2+CupThickness/2;
    G4double VirtualOverBottomPosX=CupPosX+CupThickness/2+VirtualOverBottomThickness/2;
    G4double VirtualLateralPosX=GuardRingPosX+GuardRingThickness/2+1*cm+(FaradayCupBottomThickness+CupThickness+VirtualBottomThickness)/2;
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    virtualMag= new G4Box("virtualMag", 31.*cm, 6*cm, 6*cm );
    
    logicVirtualMag= new G4LogicalVolume( virtualMag,
                                         internalChamberMaterial,
                                         "LVirtualMag",
                                         0,0,0);
    physicVirtualMag = new G4PVPlacement(0,
                                         G4ThreeVector(virtualMagPosX, 0.*cm, 0*mm),
                                         "PVirtualMag",
                                         logicVirtualMag,
                                         physicTreatmentRoom,
                                         true, 0);
    
    
    logicVirtualMag -> SetVisAttributes(blue);
    
    //// BeveledCylinder ////
    
    G4RotationMatrix *Rot= new G4RotationMatrix;
    Rot->rotateX(14*deg);
    G4ThreeVector trans(0.,22.5*mm,-15*mm);
    Cylinder= new G4Tubs("cylinder",20*mm,22.5*mm,90*mm,0.,2*pi);
    Box= new G4Box("Box",22.5*mm,22.5*mm,90*mm);
    
    G4SubtractionSolid* BeveledCylinder=new G4SubtractionSolid("Cylinder-Box",
                                                               Cylinder,
                                                               Box,
                                                               Rot,
                                                               trans);
    
    logicBeveledCylinder= new G4LogicalVolume                   (BeveledCylinder,
                                                                 GuardRingMaterial,
                                                                 "LBeveledCylinder",
                                                                 0,0,0);
    
    physicBeveledCylinder =new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(GuardRingPosX,0,0)),
                                             "physicBeveledCylinder",
                                             logicBeveledCylinder,
                                             physicVirtualMag,
                                             true,0);
    
    logicBeveledCylinder->SetVisAttributes(green);
    
    
    ///// KaptonEntranceWindow /////
    
    KaptonEntranceWindow= new G4Tubs("KaptonEntranceWindow",
                                     0,
                                     OuterRadiusFC,
                                     KaptonEntranceWindowThickness/2,
                                     0*deg,360*deg);
    
    logicKaptonEntranceWindow=new G4LogicalVolume(          KaptonEntranceWindow,
                                                  // internalChamberMaterial, for track control
                                                  KaptonEntranceWindowMaterial,
                                                  "LKaptonEntranceWindow",
                                                  0,0,0);
    
    physicKaptonEntranceWindow=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(KaptonEntranceWindowPosX,0,0)),
                                                 "PhysicEntranceWindow",
                                                 logicKaptonEntranceWindow,
                                                 physicVirtualMag,true,0);
    logicKaptonEntranceWindow -> SetVisAttributes(gray);
    
    ////// MassRing /////
    
    MassRing=new G4Tubs ("MassRing",
                         InnerRadiusFC,
                         OuterRadiusFC,
                         MassRingThickness/2,
                         0*deg,360*deg);
    
    logicMassRing=new G4LogicalVolume(                      MassRing,
                                      MassRingMaterial,
                                      "logicMassRing",
                                      0,0,0);
    
    physicMassRing=new G4PVPlacement(                       G4Transform3D(rm,G4ThreeVector(MassRingPosX,0,0)),
                                     
                                     "PhysicMassRing",logicMassRing,
                                     
                                     physicVirtualMag,
                                     true,0);
    logicMassRing -> SetVisAttributes(green);
    
    
    
    
    ///// VirtualWindow /////
    
    
    VirtualWindow=new G4Tubs("VirtualWindow",
                             0,
                             OuterRadiusFC,
                             VirtualWindowThickness/2,
                             0*deg,360*deg);
    
    logicVirtualWindow=new G4LogicalVolume(                 VirtualWindow,
                                           internalChamberMaterial,
                                           "logicVirtualWindow",
                                           0,0,0);
    
    physicVirtualWindow=new G4PVPlacement(                   G4Transform3D(rm,G4ThreeVector(VirtualWindowPosX,0,0)),
                                          
                                          "PhysicVirtualWindow",
                                          logicVirtualWindow,
                                          physicVirtualMag,
                                          true,0);
    logicVirtualWindow->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    ///// GuardRing /////
    
    GuardRing=new G4Tubs ("GuardRing",
                          InnerRadiusFC,
                          OuterRadiusFC,
                          GuardRingThickness/2,
                          0*deg,360*deg);
    
    logicGuardRing=new G4LogicalVolume(                      GuardRing,
                                       GuardRingMaterial,
                                       "logicGuardRing",
                                       0,0,0);
    
    physicGuardRing=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(GuardRingPosX,0,0)),
                                      
                                      "PhysicGuardRing", logicGuardRing,
                                      
                                      physicVirtualMag,
                                      true,0);
    logicGuardRing -> SetVisAttributes(red);
    
    
    /////VirtualMiddle /////
    
    
    VirtualMiddle=new G4Tubs ("VirtualMiddle",
                              0,
                              OuterRadiusFC,
                              VirtualMiddleThickness/2,
                              0*deg,360*deg);
    
    logicVirtualMiddle=new G4LogicalVolume(                      VirtualMiddle,
                                           internalChamberMaterial,
                                           "logicVirtualMiddle",
                                           0,0,0);
    
    physicVirtualMiddle=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(VirtualMiddlePosX,0,0)),
                                          
                                          "PhysicVirtualMiddle", logicVirtualMiddle,
                                          
                                          physicVirtualMag,
                                          true,0);
    
    logicVirtualMiddle->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    ///// FaradayCupBottom /////
    
    FaradayCupBottom=new G4Tubs ("FaradayCupBottom",
                                 InnerRadiusFC,
                                 OuterRadiusFC,
                                 FaradayCupBottomThickness/2,
                                 0*deg,360*deg);
    
    logicFaradayCupBottom=new G4LogicalVolume(             FaradayCupBottom,
                                              FaradayCupBottomMaterial,
                                              "logicFaradayCupBottom",
                                              0,0,0);
    
    physicFaradayCupBottom=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(FaradayCupBottomPosX,0,0)),
                                             "PhysicFaradayCupBottom",logicFaradayCupBottom,
                                             physicVirtualMag,
                                             true,0);
    logicFaradayCupBottom -> SetVisAttributes(yellow);
    
    
    ///// Virtual Bottom //////
    
    VirtualBottom=new G4Tubs ("VirtualBottom",
                              0,
                              OuterRadiusFC,
                              VirtualBottomThickness/2,
                              0*deg,360*deg);
    
    logicVirtualBottom=new G4LogicalVolume(                  VirtualBottom,
                                           internalChamberMaterial,
                                           "logicVirtualBottom",
                                           0,0,0);
    
    physicVirtualBottom=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(VirtualBottomPosX,0,0)),
                                          "PhysicVirtualBottom",
                                          logicVirtualBottom,
                                          physicVirtualMag,
                                          true,0);
    
    logicVirtualBottom->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    ///// Cup /////
    
    Cup=new G4Tubs ("Cup",
                    0,
                    OuterRadiusFC,
                    CupThickness/2,
                    0*deg,360*deg);
    
    logicCup=new G4LogicalVolume(                           Cup,
                                 CupMaterial,
                                 "logicCup",
                                 0,0,0);
    
    physicCup=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(CupPosX,0,0)),
                                "PhysicCup", logicCup,
                                
                                physicVirtualMag,
                                true,0);
    
    logicCup -> SetVisAttributes(darkGreen);
    
    
    ///// Virtual OverBottom /////
    
    VirtualOverBottom=new G4Tubs ("VirtualOverBottom",
                                  0,
                                  OuterRadiusFC,
                                  VirtualOverBottomThickness/2,
                                  0*deg,360*deg);
    
    logicVirtualOverBottom=new G4LogicalVolume(             VirtualOverBottom,
                                               internalChamberMaterial,
                                               "logicVirtualOverBottom",
                                               0,0,0);
    
    physicVirtualOverBottom=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(VirtualOverBottomPosX,0,0)),
                                              "PhysicVirtualOverBottom",logicVirtualOverBottom,
                                              
                                              physicVirtualMag,
                                              true,0);
    logicVirtualOverBottom->SetVisAttributes (G4VisAttributes::GetInvisible());
    
    
    ///// Virtual Lateral /////
    
    
    VirtualLateral=new G4Tubs ("VirtualLateral",
                               OuterRadiusFC,
                               OuterRadiusFC+1*um,// the VirtualLateralThickness is 1*um
                               VirtualLateralLength/2,
                               0*deg,360*deg);
    
    logicVirtualLateral=new G4LogicalVolume(                  VirtualLateral,
                                            internalChamberMaterial,
                                            "logicVirtualLateral",
                                            0,0,0);
    
    physicVirtualLateral=new G4PVPlacement(G4Transform3D(rm,G4ThreeVector(VirtualLateralPosX,0,0)),
                                           "VirtualLateral",logicVirtualLateral,
                                           
                                           physicVirtualMag,
                                           true,0);
    
    
    
    logicVirtualLateral->SetVisAttributes (G4VisAttributes::GetInvisible());
}

/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::Quadrupole()
{
    // To rotate the quadrupoles putting their axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    SQuadChamberWall = new G4Box("solidQuadChamberWall",externalChamberXSize/2., externalChamberYSize/2.,externalChamberZSize/2.);
    
    LQuadChamberWall = new G4LogicalVolume(SQuadChamberWall, externalChamberMaterial,"logicQuadChamberWall");
    
    PQuadChamberWall = new G4PVPlacement(0, G4ThreeVector(QuadChamberWallPosX, QuadChamberWallPosY, QuadChamberWallPosZ),
                                         "physQuadChamberWall", LQuadChamberWall,physicTreatmentRoom, false, 0);
    
    
    SQuadChamber = new G4Box("solidQuadChamber", internalChamberXSize/2., internalChamberYSize/2.,internalChamberZSize/2.);
    
    LQuadChamber = new G4LogicalVolume(SQuadChamber, internalChamberMaterial,"logicQuadChamber");
    
    PQuadChamber = new G4PVPlacement(0, G4ThreeVector(0.0, 0.0, 0.0),
                                     "physQuadChamber", LQuadChamber,PQuadChamberWall, false, 0);
    
    LQuadChamberWall -> SetVisAttributes(red);
    LQuadChamber -> SetVisAttributes(white);
    ///////////----------------------------Fourth Quadrupole----------------------------/////////
    SFourthTriplet = new G4Tubs("SolidTQuad", InnerRadiusTriplet, ExternalRadiusQuad,((FourthQuadThickness/2.)+1*mm),
                                startAngleQuad, spanningAngleQuad);
    
    LFourthTriplet = new G4LogicalVolume(SFourthTriplet, internalChamberMaterial,"LogicTQuad", 0, 0, 0);
    
    PFourthTriplet = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(FourthQuadXPosition, QuadYPosition, QuadZPosition)),
                                       "PhysFourthTQuad", LFourthTriplet, PQuadChamber, false, 0);
    
    solidFourthQuad = new G4Tubs("SolidQuad", InnerRadiusQuad, ExternalRadiusQuad, FourthQuadThickness/2.,
                                 startAngleQuad, spanningAngleQuad);
    
    logicFourthQuad = new G4LogicalVolume(solidFourthQuad, QuadMaterial, "LogicQuad", 0, 0, 0);
    
    physicFourthQuad = new G4PVPlacement(0, G4ThreeVector(FourthQXPosition, QYPosition, QZPosition),
                                         "PhysFourthQuad",logicFourthQuad, PFourthTriplet, false, 0);
    
    LFourthTriplet -> SetVisAttributes(yellow);
    logicFourthQuad -> SetVisAttributes(green);
    ///////////----------------------------Third Quadrupole----------------------------/////////
    SThirdTriplet = new G4Tubs("SolidTTQuad", InnerRadiusTriplet, ExternalRadiusQuad, (ThirdQuadThickness/2.)+1*mm,
                               startAngleQuad, spanningAngleQuad);
    
    LThirdTriplet = new G4LogicalVolume(SThirdTriplet, internalChamberMaterial,"LogicTTQuad", 0, 0, 0);
    
    PThirdTriplet = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(ThirdQuadXPosition, QuadYPosition, QuadZPosition)),
                                      "PhysThirdTQuad",LThirdTriplet,PQuadChamber, false, 0);
    
    solidThirdQuad = new G4Tubs("SolidTQuad", InnerRadiusQuad, ExternalRadiusQuad, ThirdQuadThickness/2.,
                                startAngleQuad, spanningAngleQuad);
    
    logicThirdQuad = new G4LogicalVolume(solidThirdQuad, QuadMaterial, "LogicTQuad", 0, 0, 0);
    
    physicThirdQuad = new G4PVPlacement(0, G4ThreeVector(ThirdQXPosition, QYPosition, QZPosition),
                                        "PhysThirdQuad",logicThirdQuad, PThirdTriplet, false, 0);
    
    LThirdTriplet -> SetVisAttributes(yellow);
    logicThirdQuad -> SetVisAttributes(green);
    ///////////----------------------------Second Quadrupole----------------------------/////////
    SSecondTriplet = new G4Tubs("SolidTSQuad", InnerRadiusTriplet, ExternalRadiusQuad, (SecondQuadThickness/2.)+1*mm,
                                startAngleQuad, spanningAngleQuad);
    
    LSecondTriplet = new G4LogicalVolume(SSecondTriplet, internalChamberMaterial,"LogicTSQuad", 0, 0, 0);
    
    PSecondTriplet = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(SecondQuadXPosition, QuadYPosition, QuadZPosition)),
                                       "PhysSecondTQuad", LSecondTriplet, PQuadChamber, false, 0);
    
    solidSecondQuad = new G4Tubs("SolidSQuad", InnerRadiusQuad, ExternalRadiusQuad, SecondQuadThickness/2.,
                                 startAngleQuad, spanningAngleQuad);
    
    logicSecondQuad = new G4LogicalVolume(solidSecondQuad, QuadMaterial, "LogicSQuad", 0, 0, 0);
    
    physicSecondQuad = new G4PVPlacement(0, G4ThreeVector(SecondQXPosition, QYPosition, QZPosition),
                                         "PhysSecondQuad", logicSecondQuad, PSecondTriplet, false, 0);
    
    LSecondTriplet -> SetVisAttributes(yellow);
    logicSecondQuad -> SetVisAttributes(green);
    ///////////----------------------------First Quadrupole----------------------------/////////
    SFirstTriplet = new G4Tubs("SolidTQuad", InnerRadiusTriplet, ExternalRadiusQuad, (FirstQuadThickness/2.)+1*mm,
                               startAngleQuad, spanningAngleQuad);
    
    LFirstTriplet = new G4LogicalVolume(SFirstTriplet, internalChamberMaterial,"LogicTQuad", 0, 0, 0);
    
    PFirstTriplet = new G4PVPlacement(G4Transform3D(rm, G4ThreeVector(FirstQuadXPosition, QuadYPosition, QuadZPosition)),
                                      "PhysFirstTQuad", LFirstTriplet, PQuadChamber, false, 0);
    
    solidFirstQuad = new G4Tubs("SolidQuad", InnerRadiusQuad, ExternalRadiusQuad, FirstQuadThickness/2.,
                                startAngleQuad, spanningAngleQuad);
    
    logicFirstQuad = new G4LogicalVolume(solidFirstQuad, QuadMaterial, "LogicQuad", 0, 0, 0);
    
    physicFirstQuad = new G4PVPlacement(0, G4ThreeVector(FirstQXPosition, QYPosition, QZPosition),
                                        "PhysFirstQuad",logicFirstQuad, PFirstTriplet, false, 0);
    
    LFirstTriplet -> SetVisAttributes(yellow);
    logicFirstQuad -> SetVisAttributes(green);
}

/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::EnergySelectorChamber()
{
    // The whole energyselector is mounted inside a
    // a vacuum chamber  (called 'ExternalChamber')
    // inside which a vacuum box is inserted.
    
    solidExternalChamber = new G4Box("ExternalChamber",
                                     externalChamberXSize/2.0,
                                     externalChamberYSize/2.0,
                                     externalChamberZSize/2.0);
    
    logicExternalChamber = new G4LogicalVolume(solidExternalChamber,
                                               externalChamberMaterial,
                                               "ExternalChamber");
    
    physicExternalChamber = new G4PVPlacement(0,
                                              G4ThreeVector(externalChamberXPosition,
                                                            externalChamberYPosition,
                                                            externalChamberZPosition),
                                              "ExternalChamber",
                                              logicExternalChamber,
                                              physicTreatmentRoom,
                                              false,
                                              0);
    
    // Visualisation of the External part
    logicExternalChamber -> SetVisAttributes(red);
    
    // This is a vacuum box inside the steel box
    solidInternalChamber = new G4Box("SInternalChamber",
                                     internalChamberXSize/2.0,
                                     internalChamberYSize/2.0,
                                     internalChamberZSize/2.0);
    
    logicInternalChamber = new G4LogicalVolume(solidInternalChamber,
                                               internalChamberMaterial,
                                               "LInternalChamber");
    
    physicInternalChamber = new G4PVPlacement(0,
                                              G4ThreeVector(0,0,0),
                                              "InternalChamber",
                                              logicInternalChamber,
                                              physicExternalChamber,
                                              false,
                                              0);
    logicInternalChamber -> SetVisAttributes(white);
}

//////////////////////////////////////////////////// Entrance pipe ///////////////////////
void LaserDrivenBeamLine::EntrancePipe()
{
    // To rotate the EntrancePipe putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidEntrancePipe = new G4Tubs("EntrancePipe",
                                   InnerRadiusEntrancePipe,
                                   ExternalRadiusEntrancePipe,
                                   EntrancePipeheight/2.,
                                   startAngleEntrancePipe,
                                   spanningAngleEntrancePipe);
    
    logicEntrancePipe = new G4LogicalVolume(solidEntrancePipe,
                                            PipeMaterial,
                                            "EntrancePipe",
                                            0,
                                            0,
                                            0);
    
    physicEntrancePipe = new G4PVPlacement(G4Transform3D(rm,
                                                         G4ThreeVector(EntrancePipeXPosition,
                                                                       EntrancePipeYPosition,
                                                                       EntrancePipeZPosition)),
                                           "EntrancePipe",
                                           logicEntrancePipe,
                                           physicTreatmentRoom,
                                           false,
                                           0);
    
    logicEntrancePipe -> SetVisAttributes(red);
    
}

//////////////////////////////////////////////////// Entrance hole ///////////////////////
void LaserDrivenBeamLine::Entrancehole()
{
    // To rotate the ExitPipe putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidEntrancehole = new G4Tubs("Entrancehole",
                                   InnerRadiusEntrancehole,
                                   ExternalRadiusEntrancehole,
                                   EntranceholeThickness/2.,
                                   startAngleEntrancehole,
                                   spanningAngleEntrancehole);
    
    logicEntrancehole = new G4LogicalVolume(solidEntrancehole,
                                            internalChamberMaterial,
                                            "Entrancehole",
                                            0,
                                            0,
                                            0);
    //the hole in the energy selector chamber
    physicEntranceholeESSChamber = new G4PVPlacement(G4Transform3D(rm,
                                                                   G4ThreeVector(EntranceholeXPosition,
                                                                                 EntranceholeYPosition,
                                                                                 EntranceholeZPosition)),
                                                     "Entrancehole",
                                                     logicEntrancehole,
                                                     physicExternalChamber,
                                                     false,
                                                     0);
    //the hole in the quadrupoles chamber
    physicEntrancehole = new G4PVPlacement(G4Transform3D(rm,
                                                         G4ThreeVector(EntranceholeQuadXPosition,
                                                                       EntranceholeYPosition,
                                                                       EntranceholeZPosition)),
                                           "EntranceholeQuad",
                                           logicEntrancehole,
                                           PQuadChamberWall,
                                           false,
                                           0);
    
    logicEntrancehole -> SetVisAttributes(skyBlue);
    
    
}
/////////////////////////////////////////////////////////////////////////////
void LaserDrivenBeamLine::Collimator()
{
    // To rotate the collimator putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    //8x82x210 mm are the collimator default dimensions
    solidCollimator = new G4Box("collimator",
                                thicknessCollimator/2.0,
                                collimatorBoxYSize/2.0,
                                collimatorBoxZSize/2.0);
    
    logicCollimator = new G4LogicalVolume(solidCollimator,
                                          collimatorMaterial,
                                          "collimator");
    
    physicCollimator = new G4PVPlacement(0,
                                         G4ThreeVector(collimatorBox_XPosition,
                                                       collimatorBox_YPosition,
                                                       collimatorBox_ZPosition),
                                         "collimator",
                                         logicCollimator,
                                         physicInternalChamber,
                                         false,
                                         0);
    
    logicCollimator -> SetVisAttributes(darkOrange3);
    
    solidCollimatorHole = new G4Tubs("CollimatorHole",
                                     innerRadiusCollimator,
                                     outerRadiusCollimator,
                                     thicknessCollimator/2.,
                                     startAngleCollimator,
                                     spanningAngleCollimator);
    
    logicCollimatorHole = new G4LogicalVolume(solidCollimatorHole,
                                              collimatorHoleMaterial,
                                              "CollimatorHole",
                                              0,
                                              0,
                                              0);
    
    physicCollimatorHole = new G4PVPlacement(G4Transform3D(rm,
                                                           G4ThreeVector(collimatorXPosition,
                                                                         collimatorYPosition,
                                                                         collimatorZPosition)),
                                             "CollimatorHole",
                                             logicCollimatorHole,
                                             physicCollimator,
                                             false,
                                             0);
    
    logicCollimatorHole -> SetVisAttributes(skyBlue);
}

/////////////////////////////////////////////////////////////////////////////
// Magnet number 1
void LaserDrivenBeamLine::Magnet_1()
{  // The positions of the external and internal partes are given as respect the external chamber.
    solidExternalMagnet_1 = new G4Box("SolidExternalMagnet_1",
                                      externalMagnet_1XSize/2.0,
                                      externalMagnet_1YSize/2.0,
                                      externalMagnet_1ZSize/2.0);
    
    logicExternalMagnet_1 = new G4LogicalVolume(solidExternalMagnet_1,
                                                externalMagnet_1Material,
                                                "LogicExternalMagnet_1");
    
    physicExternalMagnet_1 = new G4PVPlacement(0,
                                               G4ThreeVector(externalMagnet_1XPosition,
                                                             externalMagnet_2YPosition,
                                                             externalMagnet_2ZPosition),
                                               "PhysicExternalMagnet_1",
                                               logicExternalMagnet_1,
                                               physicInternalChamber,
                                               false,
                                               0);
    physicExternalMagnet_1Down = new G4PVPlacement(0,
                                                   G4ThreeVector(externalMagnet_1XPosition,
                                                                 -externalMagnet_2YPosition,
                                                                 externalMagnet_2ZPosition),
                                                   "PhysicExternalMagnet_1Down",
                                                   logicExternalMagnet_1,
                                                   physicInternalChamber,
                                                   false,
                                                   0);
    
   	
    logicExternalMagnet_1 -> SetVisAttributes(gray);
    
    // The right and left part of the magnet
    solidMagnet_1 = new G4Box("SolidMagnet_1",
                              Magnet_1XSize/2.0,
                              Magnet_1YSize/2.0,
                              Magnet_1ZSize/2.0);
    
    logicMagnet_1 = new G4LogicalVolume(solidMagnet_1,
                                        externalMagnet_1Material,
                                        "LogicMagnet_1");
    
    physicMagnet_1Right = new G4PVPlacement(0,
                                            G4ThreeVector(Magnet_1XPosition,Magnet_1YPosition,
                                                          Magnet_1ZPosition),
                                            "PhysicMagnet_1Right",
                                            logicMagnet_1,
                                            physicInternalChamber,
                                            false,
                                            0);
    physicMagnet_1Left = new G4PVPlacement(0,
                                           G4ThreeVector(Magnet_1XPosition,Magnet_1YPosition,
                                                         -Magnet_1ZPosition),
                                           "PhysicMagnet_1Left",
                                           logicMagnet_1,
                                           physicInternalChamber,
                                           false,
                                           0);
    
    logicMagnet_1 -> SetVisAttributes(gray);
}

/////////////////////////////////////////////////////////////////////////////
// Magnet number 2
void LaserDrivenBeamLine::Magnet_2()
{ // The position of the external part are given as respect the external chamber.
    
    solidExternalMagnet_2 = new G4Box("SolidExternalMagnet_2",
                                      externalMagnet_2XSize/2.0,
                                      externalMagnet_2YSize/2.0,
                                      externalMagnet_2ZSize/2.0);
    
    logicExternalMagnet_2 = new G4LogicalVolume(solidExternalMagnet_2,
                                                externalMagnet_2Material,
                                                "LogicExternalMagnet_2");
    
    physicExternalMagnet_2 = new G4PVPlacement(0,
                                               G4ThreeVector(externalMagnet_2XPosition,
                                                             externalMagnet_2YPosition,
                                                             (externalMagnet_2ZPosition+32*mm)),
                                               "PhysicExternalMagnet_2",
                                               logicExternalMagnet_2,
                                               physicInternalChamber,
                                               false,
                                               0);
    
    physicExternalMagnet_2Down = new G4PVPlacement(0,
                                                   G4ThreeVector(externalMagnet_2XPosition,
                                                                 -externalMagnet_2YPosition,
                                                                 (externalMagnet_2ZPosition+32*mm)),
                                                   "PhysicExternalMagnet_2Down",
                                                   logicExternalMagnet_2,
                                                   physicInternalChamber,
                                                   false,
                                                   0);
    
   	
    logicExternalMagnet_2 -> SetVisAttributes(gray);
    
    // The right and left part of the magnet
    solidMagnet_2 = new G4Box("SolidMagnet_2",
                              Magnet_2XSize/2.0,
                              Magnet_2YSize/2.0,
                              Magnet_2ZSize/2.0);
    
    logicMagnet_2 = new G4LogicalVolume(solidMagnet_2,
                                        externalMagnet_2Material,
                                        "LogicMagnet_2");
    
    physicMagnet_2Right = new G4PVPlacement(0,
                                            G4ThreeVector(Magnet_2XPosition,Magnet_2YPosition,
                                                          (Magnet_2ZPosition)+32*mm),
                                            "PhysicMagnet_2Right",
                                            logicMagnet_2,
                                            physicInternalChamber,
                                            false,
                                            0);
    physicMagnet_2Left = new G4PVPlacement(0,
                                           G4ThreeVector(Magnet_2XPosition,Magnet_2YPosition,
                                                         (-(Magnet_2ZPosition)+32*mm)),
                                           "PhysicMagnet_2Left",
                                           logicMagnet_2,
                                           physicInternalChamber,
                                           false,
                                           0);
    logicMagnet_2 -> SetVisAttributes(gray);
}

/////////////////////////////////////////////////////////////////////////////
// Magnet number 3
void LaserDrivenBeamLine::Magnet_3()
{ // The position of the external part are given as respect the external chamber.
    
    solidExternalMagnet_3 = new G4Box("SolidExternalMagnet_3",
                                      externalMagnet_3XSize/2.0,
                                      externalMagnet_3YSize/2.0,
                                      externalMagnet_3ZSize/2.0);
    
    logicExternalMagnet_3 = new G4LogicalVolume(solidExternalMagnet_3,
                                                externalMagnet_3Material,
                                                "LogicExternalMagnet_3");
    
    physicExternalMagnet_3 = new G4PVPlacement(0,
                                               G4ThreeVector((externalMagnet_3XPosition),
                                                             externalMagnet_3YPosition,
                                                             (externalMagnet_3ZPosition+32*mm)),
                                               "PhysicExternalMagnet_3",
                                               logicExternalMagnet_3,
                                               physicInternalChamber,
                                               false,
                                               0);
    
    physicExternalMagnet_3Down = new G4PVPlacement(0,
                                                   G4ThreeVector((externalMagnet_3XPosition),
                                                                 -externalMagnet_3YPosition,
                                                                 (externalMagnet_3ZPosition+32*mm)),
                                                   "PhysicExternalMagnet_3Down",
                                                   logicExternalMagnet_3,
                                                   physicInternalChamber,
                                                   false,
                                                   0);
    
    logicExternalMagnet_3 -> SetVisAttributes(gray);
    
    // The right and left part of the magnet
    solidMagnet_3 = new G4Box("SolidMagnet_3",
                              Magnet_3XSize/2.0,
                              Magnet_3YSize/2.0,
                              Magnet_3ZSize/2.0);
    
    logicMagnet_3 = new G4LogicalVolume(solidMagnet_3,
                                        externalMagnet_3Material,
                                        "LogicMagnet_3");
    
    physicMagnet_3Right = new G4PVPlacement(0,
                                            G4ThreeVector(Magnet_3XPosition,Magnet_3YPosition,
                                                          (Magnet_3ZPosition+32*mm)),
                                            "PhysicMagnet_3Right",
                                            logicMagnet_3,
                                            physicInternalChamber,
                                            false,
                                            0);
    physicMagnet_3Left = new G4PVPlacement(0,
                                           G4ThreeVector(Magnet_3XPosition,Magnet_3YPosition,
                                                         (-(Magnet_3ZPosition)+32*mm)),
                                           "PhysicMagnet_3Left",
                                           logicMagnet_3,
                                           physicInternalChamber,
                                           false,
                                           0);
    logicMagnet_3 -> SetVisAttributes(gray);
    
}

/////////////////////////////////////////////////////////////////////////////
// Magnet number 4
void LaserDrivenBeamLine::Magnet_4()
{ // The position of the external part are given as respect the external chamber.
    
    solidExternalMagnet_4 = new G4Box("SolidExternalMagnet_4",
                                      externalMagnet_4XSize/2.0,
                                      externalMagnet_4YSize/2.0,
                                      externalMagnet_4ZSize/2.0);
    
    logicExternalMagnet_4 = new G4LogicalVolume(solidExternalMagnet_4,
                                                externalMagnet_4Material,
                                                "LogicExternalMagnet_4");
    
    physicExternalMagnet_4 = new G4PVPlacement(0,
                                               G4ThreeVector(externalMagnet_4XPosition,
                                                             externalMagnet_4YPosition,
                                                             externalMagnet_4ZPosition),
                                               "PhysicExternalMagnet_4",
                                               logicExternalMagnet_4,
                                               physicInternalChamber,
                                               false,
                                               0);
    
    physicExternalMagnet_4Down = new G4PVPlacement(0,
                                                   G4ThreeVector(externalMagnet_4XPosition,
                                                                 -externalMagnet_4YPosition,
                                                                 externalMagnet_4ZPosition),
                                                   "PhysicExternalMagnet_4Down",
                                                   logicExternalMagnet_4,
                                                   physicInternalChamber,
                                                   false,
                                                   0);
   	
    logicExternalMagnet_4 -> SetVisAttributes(gray);
    
    // The right and left part of the magnet
    solidMagnet_4 = new G4Box("SolidMagnet_4",
                              Magnet_4XSize/2.0,
                              Magnet_4YSize/2.0,
                              Magnet_4ZSize/2.0);
    
    logicMagnet_4 = new G4LogicalVolume(solidMagnet_4,
                                        externalMagnet_4Material,
                                        "LogicMagnet_4");
    
    physicMagnet_4Right = new G4PVPlacement(0,
                                            G4ThreeVector(Magnet_4XPosition,Magnet_4YPosition,
                                                          Magnet_4ZPosition),
                                            "PhysicMagnet_4Right",
                                            logicMagnet_4,
                                            physicInternalChamber,
                                            false,
                                            0);
    physicMagnet_4Left = new G4PVPlacement(0,
                                           G4ThreeVector(Magnet_4XPosition,Magnet_4YPosition,
                                                         -Magnet_4ZPosition),
                                           "PhysicMagnet_4Left",
                                           logicMagnet_4,
                                           physicInternalChamber,
                                           false,
                                           0);
    logicMagnet_4 -> SetVisAttributes(gray);
}

/////////////////////////////////////////////////////////////////////////////
// Slit
void LaserDrivenBeamLine::Slit()
{
    solidExternalSlit = new G4Box("ExternalSlit",
                                  externalSlitXSize/2.0,
                                  externalSlitYSize/2.0,
                                  externalSlitZSize/2.0);
    
    logicExternalSlit = new G4LogicalVolume(solidExternalSlit,
                                            externalSlitMaterial,
                                            "ExternalSlit");
    
    physicExternalSlit = new G4PVPlacement(0,
                                           G4ThreeVector(externalSlitXPosition,
                                                         externalSlitYPosition,
                                                         externalSlitZPosition),
                                           "ExternalSlit",
                                           logicExternalSlit,
                                           physicInternalChamber,
                                           false,
                                           0);
   	
    logicExternalSlit -> SetVisAttributes(green);
    // The hole
    solidInternalSlit = new G4Box("InternalSlit",
                                  internalSlitXSize/2.0,
                                  internalSlitYSize/2.0,
                                  internalSlitZSize/2.0);
    
    logicInternalSlit = new G4LogicalVolume(solidInternalSlit,
                                            internalSlitMaterial,
                                            "InternalSlit");
    
    physicInternalSlit = new G4PVPlacement(0,
                                           G4ThreeVector(internalSlitXPosition,
                                                         internalSlitYPosition,
                                                         internalSlitZPosition),
                                           "InternalSlit",
                                           logicInternalSlit,
                                           physicExternalSlit,
                                           false,
                                           0);
    
    logicInternalSlit -> SetVisAttributes(skyBlue);
    
}
////////////////////////////////////// Final collimator ////////////////////////////////////////////
void LaserDrivenBeamLine::FinalCollimator()
{
    // To rotate the collimator putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidFinalCollimator = new G4Box("collimatorFinal",
                                     collimatorFinalBoxXSize/2.0,
                                     collimatorFinalBoxYSize/2.0,
                                     collimatorFinalBoxZSize/2.0);
    
    logicFinalCollimator = new G4LogicalVolume(solidFinalCollimator,
                                               FinalcollimatorMaterial,
                                               "collimatorFinal");
    
    physicFinalCollimator = new G4PVPlacement(0,
                                              G4ThreeVector(collimatorFinalBox_XPosition,
                                                            collimatorFinalBox_YPosition,
                                                            collimatorFinalBox_ZPosition),
                                              "collimatorFinal",
                                              logicFinalCollimator,
                                              physicInternalChamber,
                                              false,
                                              0);
    logicFinalCollimator -> SetVisAttributes(darkOrange3);
    
    solidFinalCollimatorHole= new G4Tubs("FinalCollimatorHole",
                                         innerRadiusFinalCollimator,
                                         outerRadiusFinalCollimator,
                                         FinalCollimatorThickness/2.,
                                         startAngleFinalCollimator,
                                         spanningAngleFinalCollimator);
    
    logicFinalCollimatorHole = new G4LogicalVolume(solidFinalCollimatorHole,
                                                   FinalcollimatorHoleMaterial,
                                                   "FinalCollimatorHole",
                                                   0,
                                                   0,
                                                   0);
    
    physicFinalCollimatorHole = new G4PVPlacement(G4Transform3D(rm,
                                                                G4ThreeVector(FinalcollimatorXPosition,
                                                                              FinalcollimatorYPosition,
                                                                              FinalcollimatorZPosition)),
                                                  "FinalCollimatorHole",
                                                  logicFinalCollimatorHole,
                                                  physicFinalCollimator,
                                                  false,
                                                  0);
    logicFinalCollimatorHole -> SetVisAttributes(skyBlue);
}
//////////////////////////// Exit Window ////////////////////////////////////////////
void LaserDrivenBeamLine::ExitWindow()
{
    // To rotate the ExitWindow putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidExitWindow = new G4Tubs("ExitWindow",
                                 InnerRadiusExitWindow,
                                 ExternalRadiusExitWindow,
                                 ExitWindowThickness/2.,
                                 startAngleExitWindow,
                                 spanningAngleExitWindow);
    
    logicExitWindow = new G4LogicalVolume(solidExitWindow,
                                          WindowMaterial,
                                          "ExitWindow",
                                          0,
                                          0,
                                          0);
    
    physicExitWindow = new G4PVPlacement(G4Transform3D(rm,
                                                       G4ThreeVector(ExitWindowXPosition,
                                                                     ExitWindowYPosition,
                                                                     ExitWindowZPosition)),
                                         "ExitWindow",
                                         logicExitWindow,
                                         physicTreatmentRoom,
                                         false,
                                         0);
    
    logicExitWindow -> SetVisAttributes(skyBlue);
    
}

//////////////////////////////////////////////////// Exit pipe ///////////////////////
void LaserDrivenBeamLine::ExitPipe()
{
    // To rotate the ExitPipe putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidExitPipe = new G4Tubs("ExitPipe",
                               InnerRadiusExitPipe,
                               ExternalRadiusExitPipe,
                               ExitPipeheight/2.,
                               startAngleExitPipe,
                               spanningAngleExitPipe);
    
    logicExitPipe = new G4LogicalVolume(solidExitPipe,
                                        PipeMaterial,
                                        "ExitPipe",
                                        0,
                                        0,
                                        0);
    
    physicExitPipe = new G4PVPlacement(G4Transform3D(rm,
                                                     G4ThreeVector(ExitPipeXPosition,
                                                                   ExitPipeYPosition,
                                                                   ExitPipeZPosition)),
                                       "ExitPipe",
                                       logicExitPipe,
                                       physicTreatmentRoom,
                                       false,
                                       0);
    
    logicExitPipe -> SetVisAttributes(red);
    
}

///////////////////////////////////// Exit hole ///////////////////////
void LaserDrivenBeamLine::Exithole()
{
    // To rotate the ExitPipe putting its axis (along X direction) parallel to the beam axis
    G4double phi = 90. *deg;
    G4RotationMatrix rm;
    rm.rotateY(phi);
    
    solidExithole = new G4Tubs("Exithole",
                               InnerRadiusExithole,
                               ExternalRadiusExithole,
                               ExitholeThickness/2.,
                               startAngleExithole,
                               spanningAngleExithole);
    
    logicExithole = new G4LogicalVolume(solidExithole,
                                        internalChamberMaterial,
                                        "Exithole",
                                        0,
                                        0,
                                        0);
    
    physicExithole = new G4PVPlacement(G4Transform3D(rm,
                                                     G4ThreeVector(ExitholeXPosition,
                                                                   ExitholeYPosition,
                                                                   ExitholeZPosition)),
                                       "Exithole",
                                       logicExithole,
                                       physicExternalChamber,
                                       false,
                                       0);
    
    logicExithole -> SetVisAttributes(skyBlue);
    
}
/////////////////////////// MESSENGER ///////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
// Disable via external macro command the Energy Selector System
void LaserDrivenBeamLine::RemoveESS()
{
    if(physicMagnet_1Left) {delete physicMagnet_1Left; delete physicMagnet_1Right; delete logicMagnet_1; delete solidMagnet_1;}
    if(physicExternalMagnet_1Down){delete physicExternalMagnet_1Down; delete physicExternalMagnet_1; delete logicExternalMagnet_1; delete solidExternalMagnet_1;}
    if(physicMagnet_2Left){delete physicMagnet_2Left; delete physicMagnet_2Right; delete logicMagnet_2; delete solidMagnet_2;}
    if(physicExternalMagnet_2Down){ delete physicExternalMagnet_2Down; delete physicExternalMagnet_2; delete logicExternalMagnet_2; delete solidExternalMagnet_2; }
    if(physicMagnet_3Left){delete physicMagnet_3Left; delete physicMagnet_3Right; delete logicMagnet_3; delete solidMagnet_3; }
    if(physicExternalMagnet_3Down){delete physicExternalMagnet_3Down; delete physicExternalMagnet_3; delete logicExternalMagnet_3; delete solidExternalMagnet_3; }
    if(physicMagnet_4Left) {delete physicMagnet_4Left; delete physicMagnet_4Right; delete logicMagnet_4; delete solidMagnet_4; }
    if(physicExternalMagnet_4Down){delete physicExternalMagnet_4Down; delete physicExternalMagnet_4; delete logicExternalMagnet_4; delete solidExternalMagnet_4; }
    if(physicCollimatorHole){delete physicCollimatorHole; delete logicCollimatorHole; delete solidCollimatorHole; }
    if(physicCollimator) {delete physicCollimator; delete logicCollimator; delete solidCollimator; }
    if(physicFinalCollimatorHole) {delete physicFinalCollimatorHole; delete logicFinalCollimatorHole; delete solidFinalCollimatorHole; }
    if(physicFinalCollimator){delete physicFinalCollimator; delete logicFinalCollimator; delete solidFinalCollimator; }
    if(physicInternalSlit){ delete physicInternalSlit; delete logicInternalSlit; delete solidInternalSlit; }
    if(physicExternalSlit){delete physicExternalSlit; delete logicExternalSlit; delete solidExternalSlit; }
    if(physicExithole){delete physicExithole; delete logicExithole; delete solidExithole;}
    if(physicExitWindow){delete physicExitWindow; delete logicExitWindow; delete solidExitWindow;}
    if(physicExitPipe){delete physicExitPipe; delete logicExitPipe; delete solidExitPipe;}
    if(physicEntranceholeESSChamber){delete physicEntranceholeESSChamber;}
    if(physicInternalChamber){delete physicInternalChamber; delete logicInternalChamber; delete solidInternalChamber;}
    if(physicExternalChamber) {delete physicExternalChamber; delete logicExternalChamber; delete solidExternalChamber;}
    if(pFieldMgr) {delete pFieldMgr;}
    
    
    
    G4cout << "****************************************************" << G4endl;
    G4cout << "************ The ESS has been disabled *************"  << G4endl;
    G4cout << "****************************************************" << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
}
// Change via external macro command the diameter of the first collimator
void LaserDrivenBeamLine::SetFirstCollimatorRadius(G4double valueR)
{
    G4double radius = valueR;
    solidCollimatorHole -> SetOuterRadius(radius);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The first collimator aperture has been modified to "<< valueR/mm <<"mm in diameter" << G4endl;
}
/////////////////////////////////////////////////////////////////////////////
// Change via external macro command the thickness of the first collimator
void LaserDrivenBeamLine::SetFirstCollimatorThickness(G4double valueC)
{
    G4double thickness = valueC/2;
    solidCollimator -> SetXHalfLength(thickness);
    solidCollimatorHole -> SetZHalfLength(thickness);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The first collimator thickness has been modified to "<< valueC/mm <<" mm in thickness" << G4endl;
}

// Change via external macro command the Z position of the first collimator hole
void LaserDrivenBeamLine::SetFirstCollimatorPositionZ(G4double valueQ)
{
    physicCollimatorHole -> SetTranslation(G4ThreeVector(0., 0., valueQ));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The first collimator has been translated to "<< valueQ/mm <<"mm (along the z axis)" << G4endl;
}

// Change via external macro command the diameter of the second collimator
void LaserDrivenBeamLine::SetSecondCollimatorRadius(G4double value)
{
    G4double radius = value;
    solidFinalCollimatorHole -> SetOuterRadius(radius);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The second collimator aperture has been modified to "<< value/mm <<"mm in diameter" << G4endl;
}

/////////////////////////////////////////////////////////////////////////////
// Change via external macro command the thickness of the second collimator
void LaserDrivenBeamLine::SetSecondCollimatorThickness(G4double value)
{
    G4double thickness = value/2;
    solidFinalCollimator -> SetXHalfLength(thickness);
    solidFinalCollimatorHole -> SetZHalfLength(thickness);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The second collimator thickness has been modified to "<< value/mm <<" mm in thickness" << G4endl;
}

// Change via external macro command the Z position of the second collimator hole
void LaserDrivenBeamLine::SetSecondCollimatorPositionZ(G4double value)
{
    physicFinalCollimatorHole -> SetTranslation(G4ThreeVector(0., 0., value));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The second collimator has been translated to "<< value/mm <<"mm (along the z axis)" << G4endl;
}
// THE SLIT MESSENGERS
/////////////////////////////////////////////////////////////////////////////
// Change the thickness of the Slit
void LaserDrivenBeamLine::SetThicknessSlit(G4double value)
{
    if (value >(10.0*mm)) {
        G4cout <<"***************************************"<< G4endl;
        G4cout <<"******This is a warning messenger******"<< G4endl;
        G4cout <<"***************************************"<< G4endl;
        G4cout <<"The maximum value of the thickness of the slit is 10 mm, your value is >10 mm." << G4endl;
        G4cout <<"The default thickness value is used, it is: " << ((solidExternalSlit -> GetXHalfLength())*2.)/mm
        << G4endl;
        G4cout <<"***************************************"<< G4endl;
        
    }
    else  {
        G4double dimension = value/2;
        solidExternalSlit -> SetXHalfLength(dimension);
        solidInternalSlit -> SetXHalfLength(dimension);
        G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
        G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
        G4cout <<"The thickness of the slit is:" << ((solidExternalSlit -> GetXHalfLength())*2.)/mm
        << G4endl;
    }
}
/////////////////////////////////////////////////////////////////////////////
// Change the hole size (in Y direction) of the Slit
void LaserDrivenBeamLine::SetSlitHoleDimensionY(G4double value)
{
    G4double hole = value/2;
    solidInternalSlit -> SetYHalfLength(hole);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The hole of the Slit has been changed in the Y direction to "<< value/mm <<" mm" <<G4endl;
}

/////////////////////////////////////////////////////////////////////////////
// Change the hole size (in Z direction) of the Slit
void LaserDrivenBeamLine::SetSlitHoleDimensionZ(G4double value)
{
    G4double hole = value/2;
    solidInternalSlit -> SetZHalfLength(hole);
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The hole of the Slit has been changed in the Z direction to "<< value/mm <<" mm" <<G4endl;
}
/////////////////////////////////////////////////////////////////////////////
// Change the Z position of the hole of the Slit
void LaserDrivenBeamLine::SetSlitHolePositionZ(G4double value)
{
    physicInternalSlit -> SetTranslation(G4ThreeVector(0., 0., value));
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
    G4cout << "The hole of the slit has been translated to "<< value/mm <<" mm (along the Z axis)" <<G4endl;
}

// QUADRUPOLES

// Disable via external macro command all quadrupoles
void LaserDrivenBeamLine::RemoveQuads()
{
    if(physicFirstQuad)
    {delete solidFirstQuad; delete logicFirstQuad; delete physicFirstQuad;delete SFirstTriplet; delete LFirstTriplet; delete PFirstTriplet;}
    if(physicSecondQuad)
    {delete solidSecondQuad; delete logicSecondQuad; delete physicSecondQuad;delete SSecondTriplet; delete LSecondTriplet; 	delete PSecondTriplet;}
    if(physicThirdQuad)
    {delete solidThirdQuad; delete logicThirdQuad; delete physicThirdQuad;delete SThirdTriplet; delete LThirdTriplet; 	delete PThirdTriplet;}
    if(physicFourthQuad)
    {delete solidFourthQuad; delete logicFourthQuad; delete physicFourthQuad;delete SFourthTriplet; delete LFourthTriplet; 	delete PFourthTriplet;}
    if(pFieldMgrQuadFourth) {delete pFieldMgrQuadFourth;}
    if(pFieldMgrQuadThird) {delete pFieldMgrQuadThird;}
    if(pFieldMgrQuadSecond) {delete pFieldMgrQuadSecond;}
    if(pFieldMgrQuadFirst) {delete pFieldMgrQuadFirst;}
    
    
    G4cout << "******************************************************************" << G4endl;
    G4cout << "************ The Quadrupoles system has been disabled *************"  << G4endl;
    G4cout << "******************************************************************" << G4endl;
    G4RunManager::GetRunManager() -> GeometryHasBeenModified();
#ifdef G4VIS_USE
    G4UImanager::GetUIpointer() -> ApplyCommand("/vis/viewer/flush");
#endif
}

