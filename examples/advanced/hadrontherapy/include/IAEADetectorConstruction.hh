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
// This is the *BASIC* version of Hadrontherapy, a Geant4-based application
// See more at: http://g4advancedexamples.lngs.infn.it/Examples/hadrontherapy
//
// Visit the Hadrontherapy web site (http://www.lns.infn.it/link/Hadrontherapy) to request 
// the *COMPLETE* version of this program, together with its documentation;
// Hadrontherapy (both basic and full version) are supported by the Italian INFN
// Institute in the framework of the MC-INFN Group
//

#ifndef IAEADetectorConstruction_H
#define IAEADetectorConstruction_H 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VisAttributes.hh"
class G4VPhysicalVolume;
class G4LogicalVolume;
class HadrontherapyDetectorROGeometry;
class PassiveProtonBeamLine;
class IAEADetectorMessenger;
class HadrontherapyModulator;
class HadrontherapyDetectorSD;

/**
 * Geometry for the IAEA benchmark
 *
 * This geometry includes two main components: water target and a
 * detector that counts the particles that go through the phantom.
 *
 */
class IAEADetectorConstruction : public G4VUserDetectorConstruction
{
public:

  IAEADetectorConstruction();

  ~IAEADetectorConstruction();

  G4VPhysicalVolume* Construct();  

private: 
  void ConstructPassiveProtonBeamLine();
  void ConstructDetector();
 
  void ConstructSensitiveDetector();

 
  //  G4VisAttributes* redWire;
  
public: 
  void setWaterThickness(G4double); //< sets thickness of water phantom, zero or negative value removes phantom and plexiedges from the simulation.
  G4double ComputeVoxelSize() {return detectorSizeX/numberOfVoxelsAlongX;};
  //<Returns the size of the voxel along the X axis
 
private:
  G4VisAttributes* skyBlue;
  G4VisAttributes* red;

  G4String emName;
  HadrontherapyDetectorSD* detectorSD; //<Pointer to sensitive detector

  HadrontherapyDetectorROGeometry* detectorROGeometry; //<Pointer to ROGeometry 

  PassiveProtonBeamLine* passiveProtonBeamLine; //<Pointer to the beam line 
  // geometry component

  HadrontherapyModulator* modulator; // Pointer to the modulator 
                                     // geometry component

  G4VPhysicalVolume* physicalTreatmentRoom;
  G4VPhysicalVolume* phantomPhysicalVolume;
  G4VPhysicalVolume* phantomEdge1PhysicalVolume;
  G4VPhysicalVolume* phantomEdge2PhysicalVolume;
  
  G4LogicalVolume* detectorLogicalVolume;
  G4LogicalVolume* beamWindowLogicalVolume; //<Logical volume for beam source window
  G4LogicalVolume* NewDetectorLogicalVolume; //<Logical volume for end-detector
  
  G4VPhysicalVolume* detectorPhysicalVolume;
  G4VPhysicalVolume* beamWindowPhysicalVolume; ///<Logical volume for end-detector
  G4VPhysicalVolume* NewDetectorPhysicalVolume; ///<Physical volume for end-detector

  G4double startDetectorThickness;
  G4double phantomCenter;
  G4double phantomDepth;
  G4double plexiThickness;
  G4double aluWindowThickness;
  G4double endDetectorThickness;
  G4double endDetectorPosition;
  G4double moveEndDetectorForward;

  IAEADetectorMessenger* detectorMessenger; 

  G4double detectorSizeX; 
  G4double detectorSizeY; 
  G4double detectorSizeZ;
   
  G4int numberOfVoxelsAlongX; 
  G4int numberOfVoxelsAlongY;
  G4int numberOfVoxelsAlongZ; 
 
  G4bool noPhantom; //<If true no water-phantom is constructed
};
#endif
