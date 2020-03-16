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

// (adapted from B1DetectorConstruction)
// A.Knaian, N.MacFadden

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "fastAerosol.hh"

#include "G4UserLimits.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class DetectorConstructionMessenger;

/// Detector construction class to define materials and geometry.

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
	DetectorConstruction();
	virtual ~DetectorConstruction();

	virtual G4VPhysicalVolume* Construct();
	
	G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

	// Physics
	G4double stepLim = DBL_MAX;			// global step limit
	G4UserLimits* stepLimits;			// physics implementation of limit

	// Cloud build choice
	G4bool fastAerosolCloud = false;
	G4bool parameterisedCloud = false;
	G4bool smoothCloud = false;

	// Cloud droplet details
	G4double dropletR = 1;

	// fastAerosol cloud details
	fastAerosol* cloud = NULL;			// the cloud bulk and droplet positions
	G4bool prePopulate = false;			// whether to pre-load droplet positions (true) or not (false)
	G4double dropletNumDens = 0;		// numver density of droplets
	G4double minSpacing = 0.0;			// minimum spacing between droplets
	G4double gridPitch = 0.0;			// width of a voxel in the fastAerosol object
	int cloudSeed = 0;					// random seed dictating droplet distribution

	// parameterised cloud details
	G4double smartless = 2.0;			// control the 'smartless' property of parameterised solid. Roughly how many voxels the volume is split into for geometry optimization

	G4String cloudShapeStr = "box";		// cloud bulk shape
	G4String dropletShapeStr = "sphere";// droplet shape

	G4VSolid* cloudShape;				// actual cloud bulk shape
	G4VSolid* dropletShape;				// actual droplet solid

  protected:
	G4LogicalVolume*  fScoringVolume;

  private:
	DetectorConstructionMessenger* messenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

