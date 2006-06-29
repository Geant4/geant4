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
// FredVoxelTest.hh
//
// Declaration of Fred's voxel solid tester
//

#ifndef FredVoxelTest_hh
#define FredVoxelTest_hh

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

class G4VSolid;
class G4VVisManager;

class FredVoxelTest {
	public:
	FredVoxelTest(  );
	~FredVoxelTest();
	
	void Test( const EAxis axis, const G4VSolid *solid );
	void Draw();
	
	void SetExtent( const EAxis axis, const G4double min, const G4double max );
	void SetOrigin( const G4ThreeVector origin );
	void Rotate( const EAxis axis, const G4double value );
	void ResetRotation();
	
	protected:
	
	void PlotMarker( G4ThreeVector point, 
			 G4VVisManager *visManager, G4bool test=false );
	void PlotLine( G4ThreeVector start, G4ThreeVector end, 
		       G4VVisManager *visManager,G4bool test=false );
	
	G4AffineTransform transform,
			  inverseTransform;
	G4VoxelLimits	  voxelLimits;
	G4RotationMatrix  rotation;
	
	typedef struct {
		const G4VSolid *solid;
		EAxis	 axis;
		G4bool	 result;
		G4double min, max;
	} FredVoxelTestResult;
	
	FredVoxelTestResult test;
};	

#endif
