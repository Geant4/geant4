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
// FastAerosol
//
// Class description:
//
//   A FastAerosol is a collection of points in a voxelized 
//   arbitrarily-shaped volume with methods implementing population of
//   grids/voxels and for efficiently finding the nearest point.
//
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)
// --------------------------------------------------------------------

#ifndef FastAerosol_h
#define FastAerosol_h

#include "globals.hh"
#include "Randomize.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"

// rotations and number density distribution
#include <functional>
#include "G4RotationMatrix.hh"
#include <atomic>

//using namespace std;

class FastAerosol {
	public:
		// Constructor; creates a random cloud of droplets
		FastAerosol(const G4String& pName, G4VSolid* pCloud,
					   G4double pR, G4double pMinD, G4double pAvgNumDens, G4double pdR,
					   std::function<G4double (G4ThreeVector)> pNumDensDistribution);

		FastAerosol(const G4String& pName, G4VSolid* pCloud,
					G4double pR, G4double pMinD, G4double pNumDens, G4double pdR);

		FastAerosol(const G4String& pName, G4VSolid* pCloud,
					G4double pR, G4double pMinD, G4double pNumDens);
		
		// Destructor; frees memory
		~FastAerosol();

		// Populate all grids. Otherwise, they are populated on-the-fly
		void PopulateAllGrids();

		// Save locations of droplets to a file for visualization/analysis purposes
		void SaveToFile(const char *filename);

		// Get absolutely nearest droplet - must be public as FastAerosolSolid uses it
		bool GetNearestDroplet(const G4ThreeVector &p, G4ThreeVector &center, G4double &closestDistance, G4double stepLim, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation);
		
		// Get nearest droplet along a vector - must be public as FastAerosolSolid uses it
		bool GetNearestDroplet(const G4ThreeVector &p, const G4ThreeVector &v, G4ThreeVector &center, G4double &closestDistance, G4double stepLim, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation);

		// ======
		// Inline
		// ======
		// Input quantities
		inline G4String GetName() const;		// aerosol name
		inline G4VSolid* GetBulk() const;		// bulk shape
		inline G4double GetRadius() const;		// droplet radius
		inline G4double GetAvgNumDens() const;	// droplet number density
		//inline G4double GetPitch() const;		// grid pitch

		inline G4int GetNumDroplets() const;	// in case the absolute number is more relevant than density

		// Bulk quantities
		inline G4double GetXHalfLength() const;
		inline G4double GetYHalfLength() const;
		inline G4double GetZHalfLength() const;
		inline void GetBoundingLimits(G4ThreeVector &pMin, G4ThreeVector &pMax) const;
		inline G4double GetCubicVolume() const;

		inline G4double DistanceToCloud(const G4ThreeVector &p);
		inline G4double DistanceToCloud(const G4ThreeVector &p, const G4ThreeVector &v);

		// Misc getters and setters
		inline long GetSeed();
		inline void SetSeed(long seed);

		inline G4int GetNumPlacementTries();
		inline void SetNumPlacementTries(G4int numTries);

		inline G4int GetPreSphereR();
		inline void SetPreSphereR(G4int fPreSphereRIn);

		inline std::function<G4double (G4ThreeVector)> GetDistribution(); // the droplet number density distribution

		inline G4double GetDropletsPerVoxel();
		inline void SetDropletsPerVoxel(G4double newDropletsPerVoxel);

		// Printing diagnostic tool
		inline void PrintPopulationReport();
		
	private:
		G4double kCarTolerance;

		// Parameters, set in constructor
		G4String fName;

		G4VSolid* fCloud;					// Solid volume of the cloud
		G4double fDx, fDy, fDz;				// Half widths

		G4double fR;						// Bounding radius of each droplet
		G4double fR2;						// Bounding radius squared of each droplet

		G4double fdR;						// Uncertainty in DistanceToIn droplet when just using knowledge of droplet center

		G4double fMinD;						// Minimum distance allowed between faces of droplets when constructing random array of droplets

		std::function<G4double (G4ThreeVector)> fDistribution;
		G4double fAvgNumDens;				// Average droplet number density
		long int fNumDroplets = 0;			// Number of droplets that have been created

		G4double fGridPitch;				// Pitch of collision detection grid.  Must be greater than diameter of droplets for correctness of collision detection.

		// Ramdom engine
		CLHEP::HepJamesRandom fCloudEngine;
		long fSeed = 0;						// Global random seed

		G4double fDropletsPerVoxel = 4.0;	// Expected number of droplets per voxel

		// How far the voxel center must be inside the bulk order for there to be no risk of placing a droplet outside
		G4double fEdgeDistance;

		// Grid variables
		std::vector<std::vector<G4ThreeVector>> fGrid;// Grid of lists of inidices to grid points, used for fast collsion checking
		std::vector<G4double> fGridMean;			// Array listing mean count for each voxel
		std::atomic<bool> *fGridValid;		// Array listing validity of each grid. uses atomic variables

		G4int fNx, fNy, fNz;				// Number of x, y, and z elements in fGrid
		G4int fNxy;							// Cached fNx*fNy
		long int fNumGridCells;				// Cached fNx*fNy*fNz


		G4double fCollisionLimit2;			// Threshold distance squared when checking for collsion
		G4int fNumNewPointTries = 100;		// How many times we try to place droplets
		G4double fMaxDropPercent = 1.0;		// The maximal percentage of skipped droplets before crashing0
		G4int fMaxDropCount;				// The maximal number of skipped droplets before crashing
		G4int fNumDropped = 0;				// Number of skipped droplets due to collisions/out of bulk placement

		G4int fNumCollisions = 0;			// How many collisions occured when attempting to place

		// Droplet search variables
		G4int fVectorSearchRadius;			// maximum vector search radius

		// Droplet placement functions
		// ===========================
		void InitializeGrid();

		G4bool FindNewPoint(G4bool edgeVoxel, G4double dX, G4double dY, G4double dZ, G4double minX, G4double minY, G4double minZ, G4ThreeVector &foundVec);

		G4double VoxelOverlap(G4ThreeVector voxelCenter, G4int nStat, G4double epsilon);

		bool CheckCollision(G4double x, G4double y, G4double z);
		bool CheckCollisionInsideGrid(G4double x, G4double y, G4double z, unsigned int xi, unsigned int yi, unsigned int zi);
		bool CheckCollisionWithDroplet(G4double x, G4double y, G4double z, G4ThreeVector p);

		// Droplet distance functions
		// ==========================
		void SearchSphere(G4int searchRad, G4double &minDistance, std::vector<G4ThreeVector> &candidates, std::vector<G4double> &distances2, G4int xGrid, G4int yGrid, G4int zGrid, const G4ThreeVector &p);
		void GetNearestDropletInsideRegion(G4double &minDistance, G4ThreeVector &center, int xGrid, int yGrid, int zGrid, int xWidth, int yWidth, int zWidth, const G4ThreeVector &p, const G4ThreeVector &v, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation);
		void GetNearestDropletInsideGrid(G4double &minDistance, std::vector<G4ThreeVector> &candidates, std::vector<G4double> &distances2, unsigned int xGrid, unsigned int yGrid, unsigned int zGrid, const G4ThreeVector &p);
		void GetNearestDropletInsideGrid(G4double &minDistance, G4ThreeVector &center, unsigned int xGrid, unsigned int yGrid, unsigned int zGrid, const G4ThreeVector &p, const G4ThreeVector &v, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation);

		// Voxelized sphere methods
		// ========================
		// a collection of points as in {{x1,y1},{x2,y2},...}
		typedef std::vector<std::vector<int>> fCircleType;

		// a collection of points describing a spherical shell
		// with points (x,y,z)=(i-R,j-R,sphere[i][j][k])
		// that is, first index gives x-position, second index
		// gives y-position, and the value gives z-position
		//
		// this is done so that searching may be optimized:
		// if searching some x=i-R that is outside the
		// aerosol's bounding box, immediately increment i
		// (similar for y).
		typedef std::vector<std::vector<std::vector<int>>> fSphereType;

		G4int fMaxCircleR;
		G4int fMaxSphereR;
		G4int fPreSphereR = 20;
		std::vector<fCircleType> fCircleCollection;
		std::vector<fSphereType> fSphereCollection;
		
		fSphereType MakeSphere(G4int R);
		fCircleType MakeCircle(G4int R);
		fCircleType MakeHalfCircle(G4int R);

		void PopulateGrid(unsigned int xi, unsigned int yi, unsigned int zi, unsigned int& gi); 

		// ======
		// Inline
		// ======
		inline bool GetGrid(const G4ThreeVector &p, G4int &xGrid, G4int &yGrid, G4int &zGrid);
		inline bool AnyIndexOutOfBounds(G4int xGrid, G4int yGrid, G4int zGrid);

		inline unsigned int GetGridIndex(unsigned int xi, unsigned int yi, unsigned int zi);
		inline G4ThreeVector GetIndexCoord(G4int index);

		inline std::pair<G4int, G4int> GetMinMaxSide(G4int index, G4int numGrids);  
};

#include "FastAerosol.icc"

#endif
