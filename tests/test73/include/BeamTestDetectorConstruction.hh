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
//
// $Id:$
// GEANT4 tag $Name:$
// 
#ifndef BEAMTESTDETECTORCONSTRUCTION_HH
#define BEAMTESTDETECTORCONSTRUCTION_HH

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"
//#include "BeamTestParameters.hh"
#include "G4ios.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VPVParameterisation;
class G4UserLimits;
class BeamTestDetectorMessenger;
class G4VSensitiveDetector;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class BeamTestDetectorConstruction : public G4VUserDetectorConstruction
{
	public:
		
		// Constructor
		explicit BeamTestDetectorConstruction(/*Parameters* parameter*/);
		// Destructor
		virtual ~BeamTestDetectorConstruction();

	public:

		G4VPhysicalVolume* Construct();

		/*const 
			G4VPhysicalVolume* GetTracker() {return physiTracker;};
		G4double GetTrackerFullLength() {return fTrackerLength;};
		G4double GetTargetFullLength()  {return fTargetLength;};
		G4double GetWorldFullLength()   {return fWorldLength;}; 
		*/

	private:
		// World logical and physical volumes


		G4Box*             solidWorld;    // pointer to the solid envelope 
		G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
		G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope

	
		G4Box*             solidTracker;  // pointer to the solid Tracker
		G4LogicalVolume*   logicTracker;  // pointer to the logical Tracker
		G4VPhysicalVolume* physiTracker;  // pointer to the physical Tracker

		G4Box*             solidChamber;  // pointer to the solid Chamber
		G4LogicalVolume*   logicChamber;  // pointer to the logical Chamber
		G4VPhysicalVolume* physChamber;  // pointer to the physical Chamber

		G4VPhysicalVolume* physiChamber;  // pointer to the physical Chamber

		G4Material*         ChamberMater; // pointer to the chamber material

		G4VPVParameterisation* chamberParam; // pointer to chamber parameterisation
		G4UserLimits* stepLimit;             // pointer to user step limits

		//ExN02MagneticField* fpMagField;   // pointer to the magnetic field 

		//ExN02DetectorMessenger* detectorMessenger;  // pointer to the Messenger

		G4double fWorldLength;            // Full length of the world volume
		G4double fTargetLength;           // Full length of Target
		G4double fTrackerLength;          // Full length of Tracker
		G4int    NbOfChambers;            // Nb of chambers in the tracker region
		G4double ChamberWidth;            // width of the chambers
		G4double ChamberSpacing;	       // distance between chambers
		G4double normalise;	       // distance between chambers
        BeamTestDetectorMessenger* messenger;   //Commands messenger
        G4VSensitiveDetector* monitor;        //The sensitive detector
    public:
        //Setters
        void SetNumberChambers(G4int num) { NbOfChambers = num; }
        void SetChamberSpacing(G4double num) { ChamberSpacing = num; }
        void SetChamberWidth(G4double num) { ChamberWidth = num; }
        void UpdateGeometry();  
        
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#endif

