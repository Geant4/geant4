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

// --------------------------------------------------------------------
// GEANT 4 class header file
//
// 
// FastAerosolSolid
//
// Class description:
//
//   A FastAerosolSolid is a collection of fDroplet solids
//   with positions set randomly by FastAerosol in a
//   volume given by a bulk shape given by FastAerosol member
//   
//   The FastAerosol member fCloud handles the optimization
//   (finding nearest droplet, populating the cloud) needed for
//   efficient simulations.
//
//   This class is heavily based on system solids.
//
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)
// --------------------------------------------------------------------

#ifndef FastAerosolSolid_HH
#define FastAerosolSolid_HH

#include "FastAerosol.hh"

#include "G4Polyhedron.hh"

// rotations
#include <functional>
#include "G4RotationMatrix.hh"

class FastAerosolSolid : public G4VSolid
{
	public: 
		FastAerosolSolid(const G4String& pName,
							   FastAerosol* pCloud,
							   G4VSolid* pDroplet,
							   std::function<G4RotationMatrix (G4ThreeVector)> pRotation);

		FastAerosolSolid(const G4String& pName,
							   FastAerosol* pCloud,
							   G4VSolid* pDroplet);
		~FastAerosolSolid();

		// Access functions
		inline G4double GetCubicVolume();
		inline G4double GetSurfaceArea();

		// Solid standard methods
		G4bool CalculateExtent(const EAxis pAxis,
							   const G4VoxelLimits& pVoxelLimit,
							   const G4AffineTransform& pTransform,
							   G4double& pmin, G4double& pmax) const;
		EInside Inside(const G4ThreeVector& p) const;
		G4ThreeVector SurfaceNormal( const G4ThreeVector& p) const;
		G4double DistanceToIn(const G4ThreeVector& p,
													const G4ThreeVector& v) const;
		G4double DistanceToIn(const G4ThreeVector& p) const;
		G4double DistanceToOut(const G4ThreeVector& p,
													 const G4ThreeVector& v,
													 const G4bool calcNorm=G4bool(false),
																 G4bool *validNorm=0,
																 G4ThreeVector *n=0) const;
		G4double DistanceToOut(const G4ThreeVector& p) const;

		G4GeometryType GetEntityType() const;

		G4VSolid* Clone() const;

		std::ostream& StreamInfo(std::ostream& os) const;

		G4ThreeVector GetPointOnSurface() const;

		G4Polyhedron* GetPolyhedron () const;
		void DescribeYourselfTo(G4VGraphicsScene& scene) const;
		G4VisExtent GetExtent() const;
		G4Polyhedron* CreatePolyhedron() const;


	public:  // without description
	 
		FastAerosolSolid(__void__&);
			//
			// Fake default constructor for usage restricted to direct object
			// persistency for clients requiring preallocation of memory for
			// persistifiable objects.

		FastAerosolSolid(const FastAerosolSolid& rhs);
		FastAerosolSolid& operator=(const FastAerosolSolid& rhs); 
			// Copy constructor and assignment operator.

		inline void SetStepLim(G4double newLim);

	private:

		inline void Initialize();
			//
			// Reset relevant values to zero

		G4double fStepLim = DBL_MAX;		// Maximum step length. Allows speed up in droplet search

		FastAerosol* fCloud;				// FastAerosol which handles brunt of work
		G4VSolid* fDroplet;					// Droplet shape
		G4VSolid* fBulk;					// Aerosol bulk

		G4double fR = 0.0;					// Droplet bounding radius

		G4double fVisDx, fVisDy, fVisDz;	// Visual extent

		G4double fCubicVolume = 0.0;		// Cubic volume of all droplets
		G4double fSurfaceArea = 0.0;		// Surface area of all droplets

		G4double farFromCloudDist;

		std::function<G4RotationMatrix (G4ThreeVector)> fRotation;	// rotation function

	protected:  // without description

		mutable G4bool fRebuildPolyhedron;
		mutable G4Polyhedron* fpPolyhedron;
};

#include "FastAerosolSolid.icc"

#endif
