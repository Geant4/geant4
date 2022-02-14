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
// Implementation for FastAerosolSolid class
// Author: A.Knaian (ara@nklabs.com), N.MacFadden (natemacfadden@gmail.com)
// --------------------------------------------------------------------

#include "FastAerosolSolid.hh"

#include "G4SystemOfUnits.hh"

// calculate extent
#include "G4BoundingEnvelope.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

// visualization
#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

// polyhedron
#include "G4AutoLock.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"

namespace
{
	G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}


///////////////////////////////////////////////////////////////////////////////
//
// Constructor
//
FastAerosolSolid::FastAerosolSolid(const G4String& pName,
										 FastAerosol* pCloud,
										 G4VSolid* pDroplet,
										 std::function<G4RotationMatrix (G4ThreeVector)> pRotation)
	: G4VSolid(pName), fCloud(pCloud), fDroplet(pDroplet), fRotation(pRotation), fRebuildPolyhedron(false), fpPolyhedron(0)
{
	// Get cloud size from fCloud
	G4ThreeVector cloudPMin, cloudPMax;
	fCloud->GetBoundingLimits(cloudPMin, cloudPMax);

	fVisDx = cloudPMax.x();
	fVisDy = cloudPMax.y();
	fVisDz = cloudPMax.z();
	
	// Check and set droplet radius
	G4double pR = fCloud->GetRadius();
	// would be nice to add a check to make sure pDroplet fits in sphere of radius pR
	fR = pR;

	fBulk = fCloud->GetBulk(); 

	farFromCloudDist = fCloud->GetPreSphereR()*fR;
}


///////////////////////////////////////////////////////////////////////////////
//
// Alternative constructor (constant rotation function)
//
FastAerosolSolid::FastAerosolSolid(const G4String& pName,
										 FastAerosol* pCloud,
										 G4VSolid* pDroplet):
	FastAerosolSolid(pName, pCloud, pDroplet,
					 [](G4ThreeVector) {return G4RotationMatrix();})
{}

///////////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
FastAerosolSolid::FastAerosolSolid( __void__& a )
	: G4VSolid(a), fCloud(nullptr), fDroplet(nullptr),
	fBulk(nullptr), fR(0.),
	fVisDx(0.), fVisDy(0.), fVisDz(0.),
	fCubicVolume(0.), fSurfaceArea(0.),
	farFromCloudDist(0.),
	fRotation([](G4ThreeVector) {return G4RotationMatrix();}),
	fRebuildPolyhedron(false), fpPolyhedron(0)
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Destructor
//
FastAerosolSolid::~FastAerosolSolid() {
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
FastAerosolSolid::FastAerosolSolid(const FastAerosolSolid &rhs)
	: G4VSolid(rhs), fCloud(rhs.fCloud), fDroplet(rhs.fDroplet),
	fBulk(rhs.fBulk), fR(rhs.fR),
	fVisDx(rhs.fVisDx), fVisDy(rhs.fVisDy), fVisDz(rhs.fVisDz),
	fCubicVolume(rhs.fCubicVolume), fSurfaceArea(rhs.fSurfaceArea),
	farFromCloudDist(rhs.farFromCloudDist),
	fRotation(rhs.fRotation),
	fRebuildPolyhedron(rhs.fRebuildPolyhedron), fpPolyhedron(rhs.fpPolyhedron)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
FastAerosolSolid &FastAerosolSolid::operator = (const FastAerosolSolid &rhs)
{
	// Check assignment to self
	//
	if (this == &rhs)
	{
		return *this;
	}

	// Copy base class data
	//
	G4VSolid::operator=(rhs);

	// Copy data
	//
	fCloud = rhs.fCloud;
	fDroplet = rhs.fDroplet;
	fBulk = rhs.fBulk;
	fR = rhs.fR;
	fVisDx = rhs.fVisDx;
	fVisDy = rhs.fVisDy;
	fVisDz = rhs.fVisDz;
	fCubicVolume = rhs.fCubicVolume;
	fSurfaceArea = rhs.fSurfaceArea;
	farFromCloudDist = rhs.farFromCloudDist;
	fRotation = rhs.fRotation;
	fRebuildPolyhedron = rhs.fRebuildPolyhedron;
	fpPolyhedron = rhs.fpPolyhedron;


	return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
//
G4bool FastAerosolSolid::CalculateExtent(const EAxis pAxis,
											const G4VoxelLimits &pVoxelLimit,
											const G4AffineTransform &pTransform,
											G4double &pMin, G4double &pMax) const
{
	// Get smallest box to fully contain the cloud of objects, not just the centers
	//
	G4ThreeVector bmin, bmax;
	fCloud->GetBoundingLimits(bmin, bmax);

	// Find extent
	//
	G4BoundingEnvelope bbox(bmin, bmax);
	return bbox.CalculateExtent(pAxis, pVoxelLimit, pTransform, pMin, pMax);
}

///////////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
//
// This function assumes the cloud has at least 1 droplet
//
EInside FastAerosolSolid::Inside(const G4ThreeVector &p) const
{
	G4ThreeVector center;
	G4double closestDistance;

	fCloud->GetNearestDroplet(p, center, closestDistance, fR, fDroplet, fRotation);

	if (closestDistance==0.0)
	{
		G4RotationMatrix irotm = fRotation(center).inverse();

		return fDroplet->Inside( irotm*(p - center) );
	}
	else
	{
		return kOutside;
	}
}

/////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
//
// This function assumes the cloud has at least 1 droplet
//
G4ThreeVector FastAerosolSolid::SurfaceNormal(const G4ThreeVector &p) const
{
	G4ThreeVector center;
	G4double closestDistance;

	fCloud->GetNearestDroplet(p, center, closestDistance, DBL_MAX, fDroplet, fRotation);

	G4RotationMatrix rotm = fRotation(center);

	return rotm*( fDroplet->SurfaceNormal( rotm.inverse()*(p - center) ) );
}

///////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
//
// This CANNOT be an underestimate
//
G4double FastAerosolSolid::DistanceToIn(const G4ThreeVector &p, const G4ThreeVector &v) const
{
	G4ThreeVector center;
	G4double closestDistance;

	if (fCloud->GetNearestDroplet(p, v, center, closestDistance, fStepLim, fDroplet, fRotation))	// if we found a droplet within fStepLim of query
	{
		return closestDistance;
	}
	else if (fCloud->DistanceToCloud(p,v)<DBL_MAX)			// if there is cloud in front of us
	{
		return 1.1*fStepLim;
	}
	else													// flying away from cloud
	{
		return kInfinity;
	}
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from outside
//
// This function assumes the cloud has at least 1 droplet
//
// This can be an underestimate
//
G4double FastAerosolSolid::DistanceToIn(const G4ThreeVector &p) const
{
	G4ThreeVector center;
	G4double closestDistance;

	G4double distanceToCloud = fBulk->DistanceToIn(p);

	if (fBulk->Inside(p)==kOutside && distanceToCloud>=farFromCloudDist)
	{
		return distanceToCloud;
	}
	else if (fCloud->GetNearestDroplet(p, center, closestDistance, fStepLim, fDroplet, fRotation))		// if we found a droplet within fStepLim of query
	{
		return closestDistance;
	}
	else
	{
		return 1.1*fStepLim;
	}
}

//////////////////////////////////////////////////////////////////////
//
// Calculate distance (<= actual) to closest surface of shape from inside
//
// Despite being a vector distance, we find the absolutely closest
// droplet to our point since we assume that p is in a droplet and p
// could be past the center
//
// This CANNOT be an underestimate
//
G4double FastAerosolSolid::DistanceToOut(const G4ThreeVector &p,
										 const G4ThreeVector &v,
										 const G4bool calcNorm,
										 G4bool *validNorm,
										 G4ThreeVector *n) const
{
	G4ThreeVector center;
	G4double distanceToIn; // should be 0

	fCloud->GetNearestDroplet(p, center, distanceToIn, fR, fDroplet, fRotation);	// if we call this function, must be inside and thus must have a droplet within fR

	G4RotationMatrix rotm = fRotation(center);
	G4RotationMatrix irotm = rotm.inverse();

	G4ThreeVector relPos = irotm*(p-center);

	if (fDroplet->Inside(relPos) == kOutside) // something went wrong... we should be inside
	{
		std::ostringstream message;
		message << std::setprecision(15) << "The particle at point p = " << p/mm << "mm"
				<< std::setprecision(15) << " called DistanceToOut(p,v) and found the closest droplet to be at center = " << center/mm << "mm"
				<< " but p is outside the droplet!";
		G4Exception("FastAerosolSolid::DistanceToOut()", "GeomSolids0002",
					FatalErrorInArgument, message);
	}

	G4double dist = fDroplet->DistanceToOut(relPos, irotm*v, calcNorm, validNorm, n);
	*n = rotm*(*n);
	*validNorm = false; // even if droplet is convex, the aerosol isn't

	return dist;
}

/////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside
//
// This can be an underestimate
//
G4double FastAerosolSolid::DistanceToOut(const G4ThreeVector &p) const
{
	G4ThreeVector center;
	G4double distanceToIn; // should be 0

	fCloud->GetNearestDroplet(p, center, distanceToIn, fR, fDroplet, fRotation);	// if we call this function, must be inside and thus must have a droplet within fR

	G4RotationMatrix irotm = fRotation(center).inverse();
	G4ThreeVector relPos = irotm*(p-center);

	if (fDroplet->Inside(relPos) == kOutside) // something went wrong... we should be inside
	{
		std::ostringstream message;
		message << "The particle at point p = " << p/mm << "mm"
				<< " called DistanceToOut(p) and found the closest droplet to be at center = " << center/mm << "mm"
				<< " but p is outside the droplet!";
		G4Exception("FastAerosolSolid::DistanceToOut()", "GeomSolids0002", 
					FatalErrorInArgument, message);
	}

	return fDroplet->DistanceToOut(relPos);
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType
//
G4GeometryType FastAerosolSolid::GetEntityType() const
{
	return G4String("FastAerosolSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// G4EntityType
//
G4VSolid* FastAerosolSolid::Clone() const
{
	return new FastAerosolSolid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream
//
std::ostream &FastAerosolSolid::StreamInfo(std::ostream &os) const
{
	os << "-----------------------------------------------------------\n"
	   << "    *** Dump for solid - " << GetName() << " ***\n"
	   << "    ===================================================\n"
	   << " Solid type: FastAerosolSolid\n"
	   << " Parameters: \n"
	   << "    numDroplets: " << fCloud->GetNumDroplets() << "\n"
	   << "    fDroplet type: " << fDroplet->GetName() << "\n"
	   << "    fDroplet parameters: \n";
	   fDroplet->StreamInfo(os);
	  os  << "-----------------------------------------------------------\n";
	return os;
}

////////////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface
//
// Currently hardcoded to look at all droplets, not just the populated ones
//
G4ThreeVector FastAerosolSolid::GetPointOnSurface() const
{
	G4ThreeVector center;
	G4double closestDistance;

	G4double fDx = fCloud->GetXHalfLength();
	G4double fDy = fCloud->GetYHalfLength();
	G4double fDz = fCloud->GetZHalfLength();

	G4ThreeVector p(2.0*fDx*G4UniformRand(),2.0*fDy*G4UniformRand(),2.0*fDz*G4UniformRand());
	p -= G4ThreeVector(fDx, fDy, fDz);

	fCloud->GetNearestDroplet(p, center, closestDistance, DBL_MAX, fDroplet, fRotation);

	return(center + fRotation(center)*fDroplet->GetPointOnSurface());
}


/////////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation
//
void FastAerosolSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
	scene.AddSolid(*this);
}

G4VisExtent FastAerosolSolid::GetExtent() const
{
	return G4VisExtent (-fVisDx, fVisDx, -fVisDy, fVisDy, -fVisDz, fVisDz);
}

G4Polyhedron* FastAerosolSolid::CreatePolyhedron () const
{
	return fBulk->CreatePolyhedron();
}


// copied from G4Ellipsoid
G4Polyhedron* FastAerosolSolid::GetPolyhedron () const
{
	if (!fpPolyhedron ||
		fRebuildPolyhedron ||
		fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
		fpPolyhedron->GetNumberOfRotationSteps())
	{
		G4AutoLock l(&polyhedronMutex);
		delete fpPolyhedron;
		fpPolyhedron = CreatePolyhedron();
		fRebuildPolyhedron = false;
		l.unlock();
	}
	return fpPolyhedron;
}