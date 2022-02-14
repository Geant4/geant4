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
// Implementation for FastAerosol class
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)
// --------------------------------------------------------------------

#include "FastAerosol.hh"

#include "G4SystemOfUnits.hh"
#include "G4GeometryTolerance.hh"

#include <numeric>	// for summing vectors with accumulate

// multithreaded safety
//#include <atomic>
#include "G4AutoLock.hh"
namespace
{
	G4Mutex gridMutex = G4MUTEX_INITIALIZER;
	G4Mutex sphereMutex = G4MUTEX_INITIALIZER;
}

///////////////////////////////////////////////////////////////////////////////
//
// Constructor - check inputs
//             - initialize grid
//             - cache values
//
FastAerosol::FastAerosol(const G4String& pName,
							   G4VSolid* pCloud,
							   G4double pR,
							   G4double pMinD,
							   G4double pAvgNumDens,
							   G4double pdR,
							   std::function<G4double (G4ThreeVector)> pDistribution)
	: fName(pName), fCloud(pCloud), fMinD(pMinD), fDistribution(pDistribution)
{
	kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();


	// get and set bounding box dimensions
	G4ThreeVector minBounding, maxBounding;
	fCloud->BoundingLimits(minBounding, maxBounding);
	G4ThreeVector halfSizeVec = 0.5*(maxBounding - minBounding);

	G4double pX = halfSizeVec[0];
	G4double pY = halfSizeVec[1];
	G4double pZ = halfSizeVec[2];

	if (pX < 2*kCarTolerance ||
		pY < 2*kCarTolerance ||
		pZ < 2*kCarTolerance)  // limit to thickness of surfaces
	{
		std::ostringstream message;
		message << "Dimensions too small for cloud: " << GetName() << "!" << G4endl
				<< "     fDx, fDy, fDz = " << pX << ", " << pY << ", " << pZ;
		G4Exception("FastAerosol::FastAerosol()", "GeomSolids0002",
					FatalException, message);
	}
	fDx = pX;
	fDy = pY;
	fDz = pZ;


	// check and set droplet radius
	if (pR < 0.0)
	{
		std::ostringstream message;
		message << "Invalid radius for cloud: " << GetName() << "!" << G4endl
				<< "        Radius must be positive."
				<< "        Inputs: pR = " << pR;
		G4Exception("FastAerosol::FastAerosol()", "GeomSolids0002",
					FatalErrorInArgument, message);
	}
	fR = pR;
	fR2 = fR*fR;


	// check and set droplet radius safety
	if (pdR < 0.0 || pdR > fR)
	{
		std::ostringstream message;
		message << "Invalid sphericality measure for cloud: " << GetName() << "!" << G4endl
				<< "        Radius uncertainty must be between 0 and fR."
				<< "        Inputs: pdR = " << pdR
				<< "        Inputs: pR = " << pR;
		G4Exception("FastAerosol::FastAerosol()", "GeomSolids0002",
					FatalErrorInArgument, message);
	}
	fdR = pdR;


	// check and set number density
	if (pAvgNumDens <= 0)
	{
		std::ostringstream message;
		message << "Invalid average number density for cloud: " << GetName() << "!" << G4endl
				<< "     pAvgNumDens = " << pAvgNumDens;
		G4Exception("FastAerosol::FastAerosol()", "GeomSolids0002",
					FatalException, message);
	}
	fAvgNumDens = pAvgNumDens;


	// Set collision limit for collsion between equal sized balls with fMinD between them
	// no droplets will be placed closer than this distance forom each other
	fCollisionLimit2 = (2*fR + fMinD)*(2*fR + fMinD);

	// set maximum number of droplet that we are allowed to skip before halting with an error
	fMaxDropCount = (G4int)floor(fAvgNumDens*(fCloud->GetCubicVolume())*(0.01*fMaxDropPercent)); 

	// initialize grid variables
	InitializeGrid(); 

	
	// begin cache of voxelized circles and spheres
	// see header for more details on these data structures
	G4AutoLock lockSphere(&sphereMutex);	// lock for multithreaded safety. Likely not needed here, but doesn't hurt

	fCircleCollection.push_back({{0, 0}});	// the R=0 circle only has one point {{x1,y1}} = {{0,0}}
	fSphereCollection.push_back({{{0}}});	// the R=0 sphere only has one point {{{z1}}}  = {{{0}}}
	
	fMaxCircleR = 0;
	fMaxSphereR = 0;

	MakeSphere(fPreSphereR);

	lockSphere.unlock();					// unlock

	// vector search radius. In terms of voxel width, how far do you search for droplets in vector search
	// you need to search a larger area if fR is larger than one grid (currently disabled)
	fVectorSearchRadius = (G4int)ceill(fR/fGridPitch);
}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative Constructor
//
// Same as standard constructor except with a uniform droplet distribution
//
FastAerosol::FastAerosol(const G4String& pName,
							   G4VSolid* pCloud,
							   G4double pR,
							   G4double pMinD,
							   G4double pNumDens,
							   G4double pdR):
	FastAerosol(pName, pCloud, pR, pMinD, pNumDens, pdR,
				[](G4ThreeVector) {return 1.0;})
{}

///////////////////////////////////////////////////////////////////////////////
//
// Alternative Constructor
//
// Same as standard constructor except with a uniform droplet distribution and
// assuming no uncertainty in sphericality
//
FastAerosol::FastAerosol(const G4String& pName,
							   G4VSolid* pCloud,
							   G4double pR,
							   G4double pMinD,
							   G4double pNumDens):
	FastAerosol(pName, pCloud, pR, pMinD, pNumDens, 0.0,
				[](G4ThreeVector) {return 1.0;})
{}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
FastAerosol::~FastAerosol() {
}

///////////////////////////////////////////////////////////////////////////////
//
// Initialize grid
//
// Sets grids to initial values and calculates expected number of droplets
// for each voxel
//
void FastAerosol::InitializeGrid() {
	// set pitch so, on average, fDropletsPerVoxel droplets are in each voxel
	fGridPitch = std::pow(fDropletsPerVoxel/fAvgNumDens,1.0/3.0);

	// if a voxel has center farther than this distance from the bulk outside,
	// we know it is fully contained in the bulk
	fEdgeDistance = fGridPitch*std::sqrt(3.0)/2.0 + fR;

	// set number of grid cells
	fNx = (G4int)ceill(2*fDx / fGridPitch);
	fNy = (G4int)ceill(2*fDy / fGridPitch);
	fNz = (G4int)ceill(2*fDz / fGridPitch);
	fNumGridCells = (long) fNx*fNy*fNz;
	fNxy = fNx*fNy;

	// try... catch... because vectors can only be so long
	// thus you can't set grid too fine
	try
	{
		std::vector<G4ThreeVector> emptyVoxel{};
		fGrid.resize(fNumGridCells, emptyVoxel);
		fGridMean.resize(fNumGridCells, 0);
		fGridValid = new std::atomic<bool>[fNumGridCells];
	}
	catch ( const std::bad_alloc&)
	{
		std::ostringstream message;
		message << "Out of memory! Grid pitch too small for cloud: " << GetName() << "!" << G4endl
				<< "        Asked for fNumGridCells = " << fNumGridCells << G4endl
				<< "        each with element of size " << sizeof(std::vector<G4ThreeVector>) << G4endl
				<< "        each with element of size " << sizeof(bool) << G4endl
				<< "        but the max size is " << fGrid.max_size() << "!";
		G4Exception("FastAerosol::FastAerosol()", "GeomSolids0002",
					FatalErrorInArgument, message);
	}
	
	// initialize grid cache
	G4ThreeVector voxelCenter;
	G4bool valid;

	for (int i=0; i<fNumGridCells; i++)
	{
		voxelCenter = fGridPitch*(GetIndexCoord(i)+G4ThreeVector(0.5,0.5,0.5)); // center of grid at index i with assuming minimum voxel is at [0,fGridPitch]x[0,fGridPitch]x[0,fGridPitch]
		voxelCenter -= G4ThreeVector(fDx,fDy,fDz);								// shift all voxels so that aerosol center is properly at the origin

		// whether or not the grid is 'finished'
		// fGridValid[i] is true if we don't plan on populating more droplets in the voxel
		// false if otherwise
		valid = (fCloud->Inside(voxelCenter) != kInside && fCloud->DistanceToIn(voxelCenter) >= fEdgeDistance);

		if (valid)	// will not populate the voxel
		{
			// have to store 'valid' in this way so that the value is atomic and respects multithreadedness
			fGridValid[i].store(true, std::memory_order_release);
		}
		else		// might need to populate the voxel
		{
			fGridValid[i].store(false, std::memory_order_release);

			// find the overlap of the voxel with the bulk
			G4double volScaling=1.0;
			G4bool edgeVoxel =  ( kInside != fCloud->Inside(voxelCenter) || fCloud->DistanceToOut(voxelCenter) < fEdgeDistance );
			if (edgeVoxel)
			{
				volScaling = VoxelOverlap(voxelCenter, 100, 0.0);
			}

			// calculates number of droplets based off of voxel center - not ideal
			fGridMean[i] = std::max(0.0, volScaling*fDistribution(voxelCenter));
		}
	}

	// must scale fGridMean[i] so that the desired number densities are actualy achieved
	// this is because no restrictions are applied to fDistribution
	G4double tempScaling = (fCloud->GetCubicVolume())*fAvgNumDens/accumulate(fGridMean.begin(), fGridMean.end(), 0.0);
	for (int i=0; i<fNumGridCells; i++) {fGridMean[i] *= tempScaling;}
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate the overlap of the voxel with the aerosol bulk
//
// The method is largely copied from G4VSolid::EstimateCubicVolume
//
G4double FastAerosol::VoxelOverlap(G4ThreeVector voxelCenter, G4int nStat, G4double epsilon) {
	G4double cenX = voxelCenter.x();
	G4double cenY = voxelCenter.y();
	G4double cenZ = voxelCenter.z();
	
	G4int iInside=0;
	G4double px,py,pz;
	G4ThreeVector p;
	G4bool in;

	if(nStat < 100)    nStat   = 100;
	if(epsilon > 0.01) epsilon = 0.01;

	for(G4int i = 0; i < nStat; i++ )
	{
		px = cenX+(fGridPitch-epsilon)*(G4UniformRand()-0.5);
		py = cenY+(fGridPitch-epsilon)*(G4UniformRand()-0.5);
		pz = cenZ+(fGridPitch-epsilon)*(G4UniformRand()-0.5);
		p  = G4ThreeVector(px,py,pz);
		in = (fCloud->Inside(p) == kInside) && (fCloud->DistanceToOut(p) >= fR);
		if(in) iInside++;    
	}

	return (double)iInside/nStat;
}


///////////////////////////////////////////////////////////////////////////////
//
// Populate all voxels
//
// Allows for partially populated clouds.
//
void FastAerosol::PopulateAllGrids() {
	unsigned int gi;

	for (int xi=0; xi<fNx; xi++)
	{
		for (int yi=0; yi<fNy; yi++)
		{
			for (int zi=0; zi<fNz; zi++)
			{
				// PopulateGrid takes unsigned arguments
				// fNx, fNy, and fNz  are signed so that some other algebra works
				// thus we iterate with xi, yi, and zi signed and then cast them
				PopulateGrid((unsigned)xi,(unsigned)yi,(unsigned)zi,gi);
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Populate a single grid
//
// If grid is already populated, does nothing.
//
void FastAerosol::PopulateGrid(unsigned int xi, unsigned int yi, unsigned int zi, unsigned int& gi) {
	gi = GetGridIndex(xi,yi,zi);
	
	// Check if this grid needs update
	bool tmpValid = fGridValid[gi].load(std::memory_order_acquire);	// have to do this weirdly because we are using atomic variables so that multithreadedness is safe
	std::atomic_thread_fence(std::memory_order_acquire);
	
	if (!tmpValid) // if true then either outside the bulk or in a voxel that is already populated
	{
		G4AutoLock lockGrid(&gridMutex);

		// Check again now that we have the lock - in case grid became valid after first check and before lock
		tmpValid = fGridValid[gi].load(std::memory_order_acquire);

		if (!tmpValid)
		{
			// uniquely set the seed to randomly place droplets
			// changing global seed gaurantees a totally new batch of seeds
			fCloudEngine.setSeed((long)gi + (long)fNumGridCells*(long)fSeed);

			// find if the voxel is near the bulk edge. In that case, we need to check whether a placed droplet is inside the bulk
			// if not an edge voxel, we can place droplets by only considering interference between with other droplets
			G4ThreeVector voxelCenter = fGridPitch*G4ThreeVector(xi+0.5,yi+0.5,zi+0.5)-G4ThreeVector(fDx,fDy,fDz);
			G4bool edgeVoxel = ( kInside != fCloud->Inside(voxelCenter) || fCloud->DistanceToOut(voxelCenter) < fEdgeDistance );
			
			// number of droplets to place
			unsigned int numDropletsToPlace = CLHEP::RandPoisson::shoot(&fCloudEngine, fGridMean[gi]);

			// actually add the points
			G4ThreeVector point;

			while (numDropletsToPlace > 0)
			{
				// find a new point not overlapping with the others
				if (FindNewPoint(edgeVoxel, fGridPitch, fGridPitch, fGridPitch, (G4double)xi*fGridPitch-fDx, (G4double)yi*fGridPitch-fDy, (G4double)zi*fGridPitch-fDz, point) )
				{
					fGrid[gi].push_back(point);
					fNumDroplets++;
				}

				
				numDropletsToPlace--;
			}

			// Memory fence to ensure sequential consistency,
			// because fGridValid is read without a lock
			// Voxel data update must be complete before updating fGridValid
			std::atomic_thread_fence(std::memory_order_release);

			fGridValid[gi].store(true, std::memory_order_release);	// set the grid as populated
		}

		lockGrid.unlock();
	}
}


///////////////////////////////////////////////////////////////////////////////
//
// Find a new droplet position in a box of full (not half) widths dX, dY, and dZ
// with minimum corner at (minX, minY, minZ).
//
// ----------------------------------------------------------------------------
//
// Ex: FindNewPoint(fGridPitch, fGridPitch, fGridPitch,
//                   (G4double)xi*fGridPitch-fDx, (G4double)yi*fGridPitch-fDy,
//                   (G4double)zi*fGridPitch-fDz);
//
// 1) CLHEP::RandFlat::shoot(&fCloudEngine)*fGridPitch shoots a point in range
//    [0, fGridPitch]
//
// 1.5) Note that the minimum x value of the (xi, yi, zi) grid is
//      xi*fGridPitch-fDx
//
// 2) xi*fGridPitch-fDx + (1) adds the minimum x value of the (xi, yi, zi)
//    voxel
//
// ----------------------------------------------------------------------------
//
// Ex: FindNewPoint(2.0*fDx, 2.0*fDy, 2.0*fDz, -fDx, -fDy, -fDz);
//
// 1) CLHEP::RandFlat::shoot(&fCloudEngine)*2.0*fDx shoots a point in
//     range [0, 2.0*fDx]
//
// 2) add -fDx to change range into [-fDx, fDx]
//
G4bool FastAerosol::FindNewPoint(G4bool edgeVoxel, G4double dX, G4double dY, G4double dZ, G4double minX, G4double minY, G4double minZ, G4ThreeVector &foundPt) {
	G4int tries = 0;	// counter of tries. Give up after fNumNewPointTries (you likely asked for too dense in this case)

	G4double x, y, z;

	G4ThreeVector point;
	G4bool placedOutside = false;	// variable whether or not the point was placed outside. Used for edgeVoxel checking

	// Generate a droplet and check if it overlaps with existing droplets
	do {
		tries++;

		if (tries > fNumNewPointTries)		// skip if we tried more than fNumNewPointTries
		{
			fNumDropped++;
			if (fNumDropped < fMaxDropCount)	// return error if we skipped more than fMaxDropCount droplets
			{
				return false;
			}

			std::ostringstream message;
			message << "Threw out too many droplets for cloud: " << GetName()  << G4endl
					<< "        Tried to place individual droplest " << fNumNewPointTries << " times." << G4endl
					<< "        This failed for " << fMaxDropCount << " droplets.";
			G4Exception("FastAerosol::FindNewPoint()", "GeomSolids0002",
						FatalErrorInArgument, message);
		}

		x = minX + CLHEP::RandFlat::shoot(&fCloudEngine)*dX;
		y = minY + CLHEP::RandFlat::shoot(&fCloudEngine)*dY;
		z = minZ + CLHEP::RandFlat::shoot(&fCloudEngine)*dZ;

		if (edgeVoxel)
		{
			point = G4ThreeVector(x,y,z);
			placedOutside = (fCloud->Inside(point) != kInside) || (fCloud->DistanceToOut(point) < fR);
		}
	}
	while (CheckCollision(x,y,z) || placedOutside);

	foundPt = G4ThreeVector(x,y,z);
	return true;  
}

///////////////////////////////////////////////////////////////////////////////
//
// Get absolutely nearest droplet
//
// Maximum search radius is maxSearch
//
// Start in grid of point p, check if any speres are in there
// and then go out one more (spherically-shaped) level to the surrounding grids
// and so on, until the whole volume has been checked or a droplet has been found
// in general as you go out, you always need to check one more level out to make sure that
// the one you have (at the closer level) is the actually the nearest one
//
bool FastAerosol::GetNearestDroplet(const G4ThreeVector &p, G4ThreeVector &center, G4double &minRealDistance, G4double maxSearch, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation)
{
	G4double cloudDistance = fCloud->DistanceToIn(p);

	// starting radius/diam of voxel layer
	G4int searchRad = (G4int)floor(0.5+cloudDistance/fGridPitch);	// no reason to search shells totally outside cloud

	// find the voxel containing p
	int xGrid, yGrid, zGrid;
	GetGrid(p, xGrid, yGrid, zGrid);    // it is OK if the starting point is outside the volume

	// initialize the containers for all candidate closest droplets and their respective distances
	std::vector<G4ThreeVector> candidates;
	std::vector<G4double> distances;

	// initialize distances to indicate that no droplet has yet been found
	G4double minDistance = kInfinity;
	G4double maxCheckDistance = kInfinity;

	bool unsafe = true;


	while (unsafe)
	{
		// limit maximum search radius
		if (searchRad*fGridPitch > maxSearch+fGridPitch) { return false; }

		SearchSphere(searchRad, minDistance, candidates, distances, xGrid, yGrid, zGrid, p);
		maxCheckDistance = minDistance+fdR;

		// theory says that, to safely have searched for centers up to some maxCheckDistance, we must search voxelized spheres at least up to ceil(0.25+maxCheckDistance/fGridPitch)
		// *** unsure if fully accurate. Calculations for 2D, not 3D... ***
		unsafe = searchRad < std::ceil(0.25+maxCheckDistance/fGridPitch);

		searchRad++;
	}


	// delete any collected droplets that are too far out. Check all candidates to see what is closest
	unsigned int index = 0;

	G4ThreeVector tempCenter;
	G4double tempDistance;

	minRealDistance = kInfinity;


	while (index < distances.size())
	{

		if (distances[index]>maxCheckDistance)
		{
			candidates.erase(candidates.begin()+index);
			distances.erase(distances.begin()+index);
		}
		else
		{
			tempCenter = candidates[index];
			tempDistance = droplet->DistanceToIn( rotation(tempCenter).inverse()*(p - tempCenter) );

			if (tempDistance < minRealDistance)
			{
				minRealDistance = tempDistance;
				center = tempCenter;
			}

			index++; // move onto next droplet
		}
	}

	return true;
}

///////////////////////////////////////////////////////////////////////////////
//
// Search sphere of radius R for droplets
//
// Saves any droplets to "center" if they are closer than minDistance
// (updating this minimum found distance at each interval).
//
// Returns if minDistance<infinity (i.e., if we have yet found a droplet)
//
void FastAerosol::SearchSphere(G4int searchRad, G4double &minDistance, std::vector<G4ThreeVector> &candidates, std::vector<G4double> &distances, G4int xGrid, G4int yGrid, G4int zGrid, const G4ThreeVector &p)
{
	// search variables
	G4int xSearch, ySearch, zSearch;		// x, y, and z coordinates of the currently searching voxel in the shell

	// voxel layer variables
	fSphereType shell;						// the shell that we are searching for droplets

	// we pre-calculate spheres up to radius fPreSphereR to speed up calculations
	// any sphere smaller than that does not need to use locks
	if (searchRad > fPreSphereR)
	{
		// need to consider making a sphere
		G4AutoLock lockSphere(&sphereMutex);

		if (searchRad > fMaxSphereR) {
			// actually need to make a sphere
			shell = MakeSphere(searchRad);
		}
		else
		{
			// we have previously made a sphere large enough
			shell = fSphereCollection[searchRad];
		}
	
		lockSphere.unlock();
	}
	else
	{
		shell = fSphereCollection[searchRad];
	}
		
	// search spherical voxel layer
	for (int i=0; i<(2*searchRad+1); i++) {
		xSearch = xGrid+i-searchRad;

		if (0<=xSearch && xSearch<fNx)
		{
			for (int j=0; j<(2*searchRad+1); j++) {
				ySearch = yGrid+j-searchRad;

				if (0<=ySearch && ySearch<fNy)
				{
					for (auto it = shell[i][j].begin(); it != shell[i][j].end(); ++it) {
						zSearch = *it+zGrid;

						if (0<=zSearch && zSearch<fNz)
						{
							GetNearestDropletInsideGrid(minDistance, candidates, distances, (unsigned)xSearch, (unsigned)ySearch, (unsigned)zSearch, p);
						}
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Make voxelized spheres up to radius R for fSphereCollection
//
// These spheres are used for searching for droplets
//
// To create them, we first create a helper half-circle "zr" for which each
// voxel/point in zr represents a layer of our sphere. The 2nd coordinate of
// zr represents the radius of the layer and the 1st coordinate denotes the
// layer's displacement from the origin along the z-axis (layers are perp to
// z in our convention).
//
std::vector<std::vector<std::vector<int>>> FastAerosol::MakeSphere(G4int R) {
	// inductively make smaller spheres since fSphereCollection organizes by index
	if (fMaxSphereR<(R-1)) {
		MakeSphere(R-1);
	}

	std::vector<int> z;								// z-coordinates of voxels in (x,y) of sphere
	std::vector<std::vector<int>> y(2*R+1, z);			// plane YZ of sphere
	fSphereType x(2*R+1, y);					// entire sphere

	x[R][R].push_back(R);						// add (0,0,R)
	x[R][R].push_back(-R);						// add (0,0,-R)

	fCircleType zr = MakeHalfCircle(R);			// break sphere into a collection of circles centered at (0,0,z) for different -R<=z<=R

	for (auto it1 = zr.begin(); it1 != zr.end(); ++it1) {
		std::vector<int> pt1 = *it1;
		fCircleType circ = MakeCircle(pt1[1]);	// make the circle

		for (auto it2 = circ.begin(); it2 != circ.end(); ++it2) {
			std::vector<int> pt2 = *it2;
			x[pt2[0]+R][pt2[1]+R].push_back(pt1[0]);
		}
	}

	fSphereCollection.push_back(x);
	fMaxSphereR = R;
	return x;
}

///////////////////////////////////////////////////////////////////////////////
//
// Make voxelized circles up to radius R for fCircleCollection
//
// These circles are used to create voxelized spheres (these spheres are then
// used for searching for droplets). This just ports the calculation to making
// a half-circle, adding the points (R,0) and (0,R), and reflecting the half-
// circle along the x-axis.
//
std::vector<std::vector<int>> FastAerosol::MakeCircle(G4int R) {
	// inductively make smaller circles since fCircleCollection organizes by index
	if (fMaxCircleR<R) {
		if (fMaxCircleR<(R-1)) {
			MakeCircle(R-1);
		}

		fCircleType voxels = {{R, 0}, {-R, 0}};		// add (R,0) and (-R,0) since half circle excludes them
		fCircleType halfVoxels = MakeHalfCircle(R);

		// join the voxels, halfVoxels, and -halfVoxels
		for (auto it = halfVoxels.begin(); it != halfVoxels.end(); ++it) {
			std::vector<int> pt = *it;
			voxels.push_back(pt);
			voxels.push_back({-pt[0], -pt[1]});
		}

		fCircleCollection.push_back(voxels);
		fMaxCircleR++;
	}
	
	return fCircleCollection[R];
}

///////////////////////////////////////////////////////////////////////////////
//
// Make half circle by a modified midpoint circle algorithm
//
// Modifies the algorithm so it doesn't have 'holes' when making many circles
// 
// See B. Roget, J. Sitaraman, "Wall distance search algorithim using 
// voxelized marching spheres," Journal of Computational Physics, Vol. 241,
// May 2013, pp. 76-94, https://doi.org/10.1016/j.jcp.2013.01.035
//
std::vector<std::vector<int>> FastAerosol::MakeHalfCircle(G4int R) {
	// makes an octant of a voxelized circle in the 1st quadrant starting at y-axis
	fCircleType voxels = {{0, R}};	// hard code inclusion of (0,R)
	int x = 1;
	int y = R;

	int dx = 4-R;		// measure of whether (x+1,y) will be outside the sphere
	int dxup = 1;		// measure of whether (x+1,y+1) will be outside the sphere
	while (x<y) {
		// mirror the points to be added to change octant->upper half plane
		voxels.push_back({ x, y});
		voxels.push_back({ y, x});
		voxels.push_back({-x, y});
		voxels.push_back({-y, x});

		// if next point outside the sphere, de-increment y
		if (dx>0) {
			// if dxup<=0, the layer above just moves to the right
			// this would create a hole at (x+1,y) unless we add it
			dxup = dx + 2*y -2*R-1;
			if ((dxup<=0) && (x+1<y)) {
				voxels.push_back({ 1+x,   y});
				voxels.push_back({   y, 1+x});
				voxels.push_back({-1-x,   y});
				voxels.push_back({  -y, 1+x});
			}

			dx += 4-2*y;
			y--;
		}

		dx += 1+2*x;
		x++;
	}

	// add points on y=+-x
	if (x==y || (x==y+1 && dxup<=0)) {
		voxels.push_back({x,x});
		voxels.push_back({-x,x});
	}

	return voxels;
}

///////////////////////////////////////////////////////////////////////////////
//
// Get nearest droplet inside a voxel at (xGrid, yGrid, zGrid)
//
void FastAerosol::GetNearestDropletInsideGrid(G4double &minDistance, std::vector<G4ThreeVector> &candidates, std::vector<G4double> &distances, unsigned int xGrid, unsigned int yGrid, unsigned int zGrid, const G4ThreeVector &p) {
	unsigned int gi;
	PopulateGrid(xGrid, yGrid, zGrid, gi);

	// initialize values
	G4double foundDistance;
//	std::vector<G4ThreeVector>::iterator bestPt;
	
	// find closest droplet
	for (auto it = fGrid[gi].begin(); it != fGrid[gi].end(); ++it) {
		foundDistance = std::sqrt((p-*it).mag2());

		if (foundDistance < minDistance+fdR) {
			minDistance = std::min(minDistance, foundDistance);
			
			candidates.push_back(*it);
			distances.push_back(foundDistance);
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Get nearest droplet along vector
//
// Maximum search radius is fVectorSearchRadius
// Maximum distance travelled along vector is maxSearch
//
// 1) Find the closest voxel along v that is inside the bulk
// 2) Search for droplets in a box width (2*fVectorSearchRadius+1) around this voxel
// that our ray would  intersect
// 3) If so, save the droplet that is closest along the ray to p and return true
// 4) If not, step 0.99*fGridPitch along v until in a new voxel
// 5) If outside bulk, jump until inside bulk and then start again from 2
// 6) If inside bulk, search the new voxels in a box centered of width (2*fVectorSearchRadius+1)
// centered at our new point
// 7) If we find a point, do as in (3). If not, do as in (4)
//
// This is searching box-shaped regions of voxels, an approach that is only valid
// for droplets which fit inside a voxel, which is one reason why the code checks
// to make sure that fR < fGridPitch at initialization.
//
bool FastAerosol::GetNearestDroplet(const G4ThreeVector &p, const G4ThreeVector &normalizedV, G4ThreeVector &center, G4double &minDistance, G4double maxSearch, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation)
{
	// get the grid box this initial position is inside
	int xGrid, yGrid, zGrid;
	GetGrid(p, xGrid, yGrid, zGrid);

	// initialize some values
	G4ThreeVector currentP = p;
	G4double distanceAlongVector = 0.0;
	G4double deltaDistanceAlongVector;

	int newXGrid = 0; int newYGrid = 0; int newZGrid = 0;
	int deltaX = 0; int deltaY = 0; int deltaZ = 0;
	int edgeX = 0 ; int edgeY = 0; int edgeZ = 0;
	G4bool boxSearch = true;
	
	G4double distanceToCloud;
	
	minDistance = kInfinity;

	// actual search
	while(true) {
		deltaDistanceAlongVector = 0.0;

		// Jump gaps
		distanceToCloud = DistanceToCloud(currentP,normalizedV);
		if (distanceToCloud == kInfinity)
		{
			return(minDistance != kInfinity);
		}
		else if (distanceToCloud > fGridPitch)
		{
			deltaDistanceAlongVector = distanceToCloud-fGridPitch;
			distanceAlongVector += deltaDistanceAlongVector;
			boxSearch = true;									// if jumping gap, use box search
			currentP += deltaDistanceAlongVector*normalizedV;	// skip gaps
		}

		// Search for droplets
		if (distanceAlongVector > maxSearch)					// quit if will jump too far
		{
			return(minDistance != kInfinity);
		}
		else
		{
			if (boxSearch)	// we jumped
			{
				GetGrid(currentP, xGrid, yGrid, zGrid);
				GetNearestDropletInsideRegion(minDistance, center,
											  xGrid, yGrid, zGrid,
											  fVectorSearchRadius, fVectorSearchRadius, fVectorSearchRadius,
											  p, normalizedV,
											  droplet, rotation);
			}
			else			// didn't jump
			{
				// Searching endcaps
				// =================
				if (deltaX != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  edgeX, yGrid, zGrid,			
												  0, fVectorSearchRadius, fVectorSearchRadius,
												  p, normalizedV,
												  droplet, rotation);
				}

				if (deltaY != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  xGrid, edgeY, zGrid,
												  fVectorSearchRadius, 0, fVectorSearchRadius,
												  p, normalizedV,
												  droplet, rotation);
				}

				if (deltaZ != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  xGrid, yGrid, edgeZ,
												  fVectorSearchRadius, fVectorSearchRadius, 0,
												  p, normalizedV,
												  droplet, rotation);
				}

				// Search bars
				// ===========
				// (a bar joins two endcaps)
				if (deltaX != 0 && deltaY != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  edgeX, edgeY, zGrid,
												  0, 0, fVectorSearchRadius,
												  p, normalizedV,
												  droplet, rotation);
				}

				if (deltaX != 0 && deltaZ != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  edgeX, yGrid, edgeZ,
												  0, fVectorSearchRadius, 0,
												  p, normalizedV,
												  droplet, rotation);
				}

				if (deltaY != 0 && deltaZ != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  xGrid, edgeY, edgeZ,
												  fVectorSearchRadius, 0, 0,
												  p, normalizedV,
												  droplet, rotation);
				}

				// Search corner
				// =============
				if (deltaX != 0 && deltaY != 0 && deltaZ != 0)
				{
					GetNearestDropletInsideRegion(minDistance, center,
												  edgeX, edgeY, edgeZ,
												  0, 0, 0,
												  p, normalizedV,
												  droplet, rotation);
				}

				// Update position
				// ================================
				xGrid = newXGrid; yGrid = newYGrid; zGrid = newZGrid;
			}
			
			// Check if we are done
			// ====================
			if (distanceAlongVector>minDistance+fdR) {
				return(true);
			}
			

			// walk along the grid
			// advance by 0.99 fGridPitch  (0.99 rather than 1.0 to avoid issues caused by 
			// rounding errors for lines parallel to the grid and based on a grid boundary)
			while (true) {
				deltaDistanceAlongVector = fGridPitch*0.99;
				distanceAlongVector += deltaDistanceAlongVector;

				// exit returning false if we have traveled more than MaxVectorFollowingDistance
				if (distanceAlongVector > maxSearch) { return(false); }

				currentP += deltaDistanceAlongVector*normalizedV;
				GetGrid(currentP, newXGrid, newYGrid, newZGrid);

				if ((newXGrid != xGrid) || (newYGrid != yGrid) || (newZGrid != zGrid)) {
					deltaX = newXGrid - xGrid; edgeX = xGrid+deltaX*(1+fVectorSearchRadius);
					deltaY = newYGrid - yGrid; edgeY = yGrid+deltaY*(1+fVectorSearchRadius);
					deltaZ = newZGrid - zGrid; edgeZ = zGrid+deltaZ*(1+fVectorSearchRadius);

					break;
				}
			}
			boxSearch = false;
		}
	}

	return false;
}

///////////////////////////////////////////////////////////////////////////////
//
// Get nearest droplet (along a vector normalizedV) inside a voxel centered at
// (xGridCenter, yGridCenter, zGridCenter) of width (xWidth, yWidth, zWidth)
//
// This searches box-shaped regions.
//
void FastAerosol::GetNearestDropletInsideRegion(G4double &minDistance, G4ThreeVector &center, int xGridCenter, int yGridCenter, int zGridCenter, int xWidth, int yWidth, int zWidth, const G4ThreeVector &p, const G4ThreeVector &normalizedV, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation) {
	for (int xGrid = (xGridCenter-xWidth); xGrid<=(xGridCenter+xWidth); xGrid++)
	{
		if (0<=xGrid && xGrid<fNx)
		{
			for (int yGrid = (yGridCenter-yWidth); yGrid<=(yGridCenter+yWidth); yGrid++)
			{
				if (0<=yGrid && yGrid<fNy)
				{
					for (int zGrid = (zGridCenter-zWidth); zGrid<=(zGridCenter+zWidth); zGrid++)
					{
						if (0<=zGrid && zGrid<fNz)
						{
							GetNearestDropletInsideGrid(minDistance, center, (unsigned)xGrid, (unsigned)yGrid, (unsigned)zGrid, p, normalizedV, droplet, rotation);
						}
					}
				}
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Get nearest droplet inside a voxel at (xGrid, yGrid, zGrid) along a vector
//
// Return the closest one, as measured along the line
//
void FastAerosol::GetNearestDropletInsideGrid(G4double &minDistance, G4ThreeVector &center, unsigned int xGrid, unsigned int yGrid, unsigned int zGrid, const G4ThreeVector &p, const G4ThreeVector &normalizedV, G4VSolid* droplet, std::function<G4RotationMatrix (G4ThreeVector)> rotation) {
	unsigned int gi;
	PopulateGrid(xGrid, yGrid, zGrid, gi);

	G4double foundDistance;
	
	// find closest droplet
	for (auto it = fGrid[gi].begin(); it != fGrid[gi].end(); ++it) {
		// could have the following check to see if the ray pierces the bounding sphere. Currently seems like unnecessary addition
		/*
		deltaP = *it-p;
		foundDistance = normalizedV.dot(deltaP);

		if ((0<=foundDistance) && ((deltaP - normalizedV*foundDistance).mag2() < fR2)) {
			
		}
		*/

		G4RotationMatrix irotm = rotation(*it).inverse();

		G4ThreeVector relPos = irotm*(p - *it);

		if (droplet->Inside(relPos) == kInside)
		{
			minDistance = 0;
			center = *it;
		}
		else
		{
			foundDistance = droplet->DistanceToIn( relPos , irotm*normalizedV);

			if (foundDistance<minDistance)
			{
				minDistance = foundDistance;
				center = *it;
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////
//
// Check if a droplet at the indicated point would collide with an existing droplet
//
// Search neighboring voxels, too. Returns true if there is a collision
//
bool FastAerosol::CheckCollision(G4double x, G4double y, G4double z) {
	G4ThreeVector p(x,y,z);
	
	std::pair<int, int> minMaxXGrid, minMaxYGrid, minMaxZGrid;

	int xGrid, yGrid, zGrid;
	GetGrid(p, xGrid, yGrid, zGrid);

	minMaxXGrid = GetMinMaxSide(xGrid, fNx);
	minMaxYGrid = GetMinMaxSide(yGrid, fNy);
	minMaxZGrid = GetMinMaxSide(zGrid, fNz);

	for (int xi=minMaxXGrid.first; xi <= minMaxXGrid.second; xi++) {
		for (int yi = minMaxYGrid.first; yi <= minMaxYGrid.second; yi++) {
			for (int zi = minMaxZGrid.first; zi <= minMaxZGrid.second; zi++) {
				if (CheckCollisionInsideGrid(x, y, z, (unsigned)xi, (unsigned)yi, (unsigned)zi)) {
					fNumCollisions++;  // log number of collisions for statistics print
					return(true);
				}
			}
		}
	}
	return(false);
}

///////////////////////////////////////////////////////////////////////////////
//
// Check if a droplet at the indicated point would collide with an existing
// droplet inside a specific grid
//
// Note that you don't need to lock the mutex since this is only called by code
// that already has the mutex (always called by PopulateGrid).
//
bool FastAerosol::CheckCollisionInsideGrid(G4double x, G4double y, G4double z, unsigned int xi, unsigned int yi, unsigned int zi) {
	std::vector<G4ThreeVector> *thisGrid = &(fGrid[GetGridIndex(xi, yi, zi)]);
	unsigned int numel = (unsigned int)thisGrid->size();

	for (unsigned int i=0; i < numel; i++) {
		if (CheckCollisionWithDroplet(x, y, z, (*thisGrid)[i]))
			return(true);
	}

	return(false);
}

///////////////////////////////////////////////////////////////////////////////
//
// Check for collsion with a specific droplet
//
bool FastAerosol::CheckCollisionWithDroplet(G4double x, G4double y, G4double z, G4ThreeVector p ) {
	return( std::pow(x-p.x(), 2.0) + std::pow(y-p.y(), 2.0) + std::pow(z-p.z(), 2.0) < fCollisionLimit2 );
}

///////////////////////////////////////////////////////////////////////////////
//
// Save droplet positions to a file for visualization and analysis
//
void FastAerosol::SaveToFile(const char* filename) {
	G4cout << "Saving droplet positions to " << filename << "..." << G4endl;
	std::ofstream file;
	file.open(filename);

	std::vector<G4ThreeVector> voxel;
	G4ThreeVector pt;

	for (auto it1 = fGrid.begin(); it1 != fGrid.end(); ++it1) {
		voxel = *it1;

		for (auto it2 = voxel.begin(); it2 != voxel.end(); ++it2) {
			pt = *it2;
			
			double x = pt.x();
			double y = pt.y();
			double z = pt.z();

			file << std::setprecision(15) << x/mm << "," << y/mm << "," << z/mm << "\n";
		}
	}

	file.close();
}
