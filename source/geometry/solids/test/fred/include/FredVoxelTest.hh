//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
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
