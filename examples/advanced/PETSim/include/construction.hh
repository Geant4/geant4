#ifndef CONSTRUCTION_HH
#define CONSTRUCTION_HH

#include "G4SystemOfUnits.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4PVPlacement.hh"
#include "G4NistManager.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"

#include "detector.hh"

class MyDetectorConstruction : public G4VUserDetectorConstruction
{
public:
    MyDetectorConstruction();
    ~MyDetectorConstruction();

    virtual G4VPhysicalVolume *Construct();

private:
    G4Box *solidWorld;
    G4Box *solidScintillator;
    G4Box *solidDetector;
    G4Tubs *solidBody;

    G4LogicalVolume *logicWorld;
    G4LogicalVolume *logicScintillator;
    G4LogicalVolume *logicDetector;
    G4LogicalVolume *logicBody;

    G4VPhysicalVolume *physWorld;
    G4VPhysicalVolume *physScintillator;
    G4VPhysicalVolume *physDetector;
    G4VPhysicalVolume *physBody;
    
    G4OpticalSurface *mirrorSurface;
    
    G4double xWorld;
    G4double yWorld;
    G4double zWorld;

    G4double length;
    G4double gap;
    G4double radius;

    G4double xCryst;
    G4double yCryst;
    G4double zCryst;

    G4double xDet;
    G4double yDet;
    G4double zDet;

    G4Material *worldMat;
    G4Material *NaI;
    G4Material *softTissue;

    G4Element *Na;
    G4Element *I;
    
    void DefineMaterials();
    virtual void ConstructSDandField();
};

#endif