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
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "globals.hh"

#include "FastAerosol.hh"

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
	G4double fStepLim = DBL_MAX;			// global step limit
	G4UserLimits* fStepLimits;			// physics implementation of limit

	// Cloud build choice
	G4bool fFastAerosolCloud = false;
	G4bool fParameterisedCloud = false;
	G4bool fSmoothCloud = false;

	// Cloud droplet details
	G4double fDropletR = 1;

	// FastAerosol cloud details
	FastAerosol* fCloud = NULL;				// the cloud bulk and droplet positions
	G4bool fPrePopulate = false;			// whether to pre-load droplet positions (true) or not (false)
	G4double fDropletNumDens = 0;			// number density of droplets
	G4double fMinSpacing = 0.0;				// minimum spacing between droplets
	int fCloudSeed = 0;						// random seed dictating droplet distribution

	// parameterised cloud details
	G4double fSmartless = 2.0;				// control the 'fSmartless' property of parameterised solid. Roughly how many voxels the volume is split into for geometry optimization

	G4String fCloudShapeStr = "box";		// cloud bulk shape
	G4String fDropletShapeStr = "sphere";	// droplet shape

	G4VSolid* fCloudShape;					// actual cloud bulk shape
	G4VSolid* fDropletShape;				// actual droplet solid

  protected:
	G4LogicalVolume*  fScoringVolume;

  private:
	DetectorConstructionMessenger* fMessenger;
};

#endif

