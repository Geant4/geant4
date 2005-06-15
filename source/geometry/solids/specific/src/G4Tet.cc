//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of     *
// * the Vanderbilt University Free Electron Laser Center             *
// * Vanderbilt University, Nashville, TN, USA                        *
// * development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P                        *
// * and was written by Marcus H. Mendenhall and Robert A. Weller     *
// *                                                                  *
// * contributed to the Geant4 Core, January, 2005                    *
// *                                                                  *
// ********************************************************************
////
// Implementation for G4Tet class
//  20040903 - Marcus Mendenhall, created G4Tet
//  20041101 - Marcus Mendenhall, optimized constant dot products with fCdotNijk values
//  20041101 - MHM removed tracking error by clipping DistanceToOut to 0 for surface cases
//  20041101 - MHM many speed optimizations in if statements
//  20041101 - MHM changed vdotn comparisons to 1e-12 instead of 0.0 to avoid nearly-parallel problems
//  20041102 - MHM Added extra distance into solid to DistanceToIn(p,v) hit testing
//  20041102 - MHM added ability to check for degeneracy without throwing G4Exception
//  20041103 - MHM removed many unused variables from class
// --------------------------------------------------------------------
#include "G4Tet.hh"

const char G4Tet::CVSVers[]="$Id: G4Tet.cc,v 1.1 2005-06-15 14:45:56 japost Exp $";

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

// #include "G4Trd.hh"
// #include "G4Trap.hh"

#include "G4ThreeVector.hh"

#include <math.h>

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4Tets are created
// a Tet has all of its geometrical infomation precomputed
G4Tet::G4Tet(const G4String& pName,
			 G4ThreeVector anchor,
			 G4ThreeVector p2,
			 G4ThreeVector p3,
			 G4ThreeVector p4, G4bool *degeneracyFlag)
: G4VSolid(pName), warningFlag(0)
{
	
	// fV<x><y> is vector from vertex <y> to vertex <x>
	G4ThreeVector fV21=p2-anchor;
	G4ThreeVector fV31=p3-anchor;
	G4ThreeVector fV41=p4-anchor;
	
	// make sure this is a correctly oriented set of points for the tetrahedron
	G4double signed_vol=fV21.cross(fV31).dot(fV41);
	if(signed_vol<0.0) {
		G4ThreeVector temp(p4);
		p4=p3;
		p3=temp;
		temp=fV41;
		fV41=fV31;
		fV31=temp; 
	}
	
	G4ThreeVector fV24=p2-p4;
	G4ThreeVector fV43=p4-p3;
	G4ThreeVector fV32=p3-p2;
	
	fXMin=fmin(fmin(fmin(anchor.x(), p2.x()),p3.x()),p4.x());
	fXMax=fmax(fmax(fmax(anchor.x(), p2.x()),p3.x()),p4.x());
	fYMin=fmin(fmin(fmin(anchor.y(), p2.y()),p3.y()),p4.y());
	fYMax=fmax(fmax(fmax(anchor.y(), p2.y()),p3.y()),p4.y());
	fZMin=fmin(fmin(fmin(anchor.z(), p2.z()),p3.z()),p4.z());
	fZMax=fmax(fmax(fmax(anchor.z(), p2.z()),p3.z()),p4.z());
	
	fDx=(fXMax-fXMin)*0.5; fDy=(fYMax-fYMin)*0.5; fDz=(fZMax-fZMin)*0.5;

	fMiddle=G4ThreeVector(fXMax+fXMin, fYMax+fYMin, fZMax+fZMin)*0.5;
	fMaxSize=fmax(fmax(fmax((anchor-fMiddle).mag(), (p2-fMiddle).mag()),
				(p3-fMiddle).mag()),(p4-fMiddle).mag());
	
	G4bool degenerate=fabs(signed_vol) < 1e-9*fMaxSize*fMaxSize*fMaxSize;
	
	if(degeneracyFlag) *degeneracyFlag=degenerate;
	else if (degenerate) {
		G4Exception("G4Tet Constructor:  constructing degenerate tetrahedron");
	}
	
	fTol=1e-9*(fabs(fXMin)+fabs(fXMax)+fabs(fYMin)+fabs(fYMax)+fabs(fZMin)+fabs(fZMax));
	//fTol=kCarTolerance;
		
	fAnchor=anchor;
	fP2=p2;
	fP3=p3;
	fP4=p4;

	G4ThreeVector fCenter123=(anchor+p2+p3)*(1.0/3.0); // face center
	G4ThreeVector fCenter134=(anchor+p4+p3)*(1.0/3.0);
	G4ThreeVector fCenter142=(anchor+p4+p2)*(1.0/3.0);
	G4ThreeVector fCenter234=(p2+p3+p4)*(1.0/3.0);
	
	fNormal123=fV31.cross(fV21).unit();
	fNormal134=fV41.cross(fV31).unit();
	fNormal142=fV21.cross(fV41).unit();
	fNormal234=fV32.cross(fV43).unit();
	
	fCdotN123=fCenter123.dot(fNormal123);
	fCdotN134=fCenter134.dot(fNormal134);
	fCdotN142=fCenter142.dot(fNormal142);
	fCdotN234=fCenter234.dot(fNormal234);
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4Tet::~G4Tet()
{
}

//////////////////////////////////////////////////////////////////////////////

G4bool G4Tet::CheckDegeneracy(
			 G4ThreeVector anchor,
			 G4ThreeVector p2,
			 G4ThreeVector p3,
			 G4ThreeVector p4)
{
	G4bool result;
	G4Tet *object=new G4Tet("temp",anchor,p2,p3,p4,&result);
	delete object;
	return result;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.

void G4Tet::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Tet::CalculateExtent(const EAxis pAxis,
                              const G4VoxelLimits& pVoxelLimit,
                              const G4AffineTransform& pTransform,
							  G4double& pMin, G4double& pMax) const
{
		
    G4double xMin,xMax;
    G4double yMin,yMax;
    G4double zMin,zMax;
	
	if (pTransform.IsRotated())
	{
		G4ThreeVector pp0=pTransform.TransformPoint(fAnchor);
		G4ThreeVector pp1=pTransform.TransformPoint(fP2);
		G4ThreeVector pp2=pTransform.TransformPoint(fP3);
		G4ThreeVector pp3=pTransform.TransformPoint(fP4);
		
		xMin    = fmin(fmin(fmin(pp0.x(), pp1.x()),pp2.x()),pp3.x());
		xMax    = fmax(fmax(fmax(pp0.x(), pp1.x()),pp2.x()),pp3.x());
		yMin    = fmin(fmin(fmin(pp0.y(), pp1.y()),pp2.y()),pp3.y());
		yMax    = fmax(fmax(fmax(pp0.y(), pp1.y()),pp2.y()),pp3.y());		
		zMin    = fmin(fmin(fmin(pp0.z(), pp1.z()),pp2.z()),pp3.z());
		zMax    = fmax(fmax(fmax(pp0.z(), pp1.z()),pp2.z()),pp3.z());
		
	} else {
		G4double xoffset = pTransform.NetTranslation().x() ;
		xMin    = xoffset + fXMin;
		xMax    = xoffset + fXMax;
		G4double yoffset = pTransform.NetTranslation().y() ;
		yMin    = yoffset + fYMin;
		yMax    = yoffset + fYMax;
		G4double zoffset = pTransform.NetTranslation().z() ;
		zMin    = zoffset + fZMin;
		zMax    = zoffset + fZMax;		
	}
		
    if (pVoxelLimit.IsXLimited())
    {
		if ( xMin > pVoxelLimit.GetMaxXExtent()+fTol || 
			 xMax < pVoxelLimit.GetMinXExtent()-fTol    ) return false ;
		else
		{
			xMin = fmax(xMin, pVoxelLimit.GetMinXExtent());
			xMax = fmin(xMax, pVoxelLimit.GetMaxXExtent());
		}
    }

    if (pVoxelLimit.IsYLimited())
    {
		if ( yMin > pVoxelLimit.GetMaxYExtent()+fTol ||
			 yMax < pVoxelLimit.GetMinYExtent()-fTol   ) return false ;
		else
		{
			yMin = fmax(yMin, pVoxelLimit.GetMinYExtent()) ;
			yMax = fmin(yMax, pVoxelLimit.GetMaxYExtent()) ;
		}
    }

    if (pVoxelLimit.IsZLimited())
    {
		if ( zMin > pVoxelLimit.GetMaxZExtent()+fTol ||
			 zMax < pVoxelLimit.GetMinZExtent()-fTol   ) return false ;
		else
		{
			zMin = fmax(zMin, pVoxelLimit.GetMinZExtent()) ;
			zMax = fmin(zMax, pVoxelLimit.GetMaxZExtent()) ;
		}
    }
	
    switch (pAxis)
    {
		case kXAxis:
			pMin = xMin ;
			pMax = xMax ;
			break ;
		case kYAxis:
			pMin=yMin;
			pMax=yMax;
			break;
		case kZAxis:
			pMin=zMin;
			pMax=zMax;
			break;
		default:
			break;
    }
	
    return true;
} 

/////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

EInside G4Tet::Inside(const G4ThreeVector& p) const
{
	
	G4double r123, r134, r142, r234;
	// this is written to allow if-statement truncation so the outside test
	// (where most of the world is) can fail very quickly and efficiently
	if ( 
		(r123=p.dot(fNormal123)-fCdotN123) > fTol ||
		(r134=p.dot(fNormal134)-fCdotN134) > fTol ||
		(r142=p.dot(fNormal142)-fCdotN142) > fTol ||
		(r234=p.dot(fNormal234)-fCdotN234) > fTol ) return kOutside; // at least one is out!
	else if(r123 < -fTol && r134 < -fTol && r142 < -fTol && r234 < -fTol) 
		return kInside; // all are definitively inside
	else return kSurface; // too close to tell
}

///////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If two sides are equidistant, normal of first side (x/y/z) 
// encountered returned
// this assumes that we are looking from the inside!

G4ThreeVector G4Tet::SurfaceNormal( const G4ThreeVector& p) const
{
	
	G4double r123=fabs(p.dot(fNormal123)-fCdotN123);
	G4double r134=fabs(p.dot(fNormal134)-fCdotN134);
	G4double r142=fabs(p.dot(fNormal142)-fCdotN142);
	G4double r234=fabs(p.dot(fNormal234)-fCdotN234);
		
	if(r123 <= r134 && r123 <= r142 && r123 <= r234) return fNormal123;
	else if (r134 <=r142 && r134 <= r234) return fNormal134;
	else if (r142 <= r234) return fNormal142;
	else return fNormal234;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate distance to box from an outside point
// - return kInfinity if no intersection.
// all this is very unrolled, for speed
G4double G4Tet::DistanceToIn(const G4ThreeVector& p,const G4ThreeVector& v) const
{
		
		G4ThreeVector vu(v.unit()), hp;
		G4double vdotn, t, tmin=kInfinity;
		
		G4double extraDistance=10.0*fTol; // a little ways into the solid
		
		vdotn=-vu.dot(fNormal123);
		if(vdotn > 1e-12) { // this is a candidate face, since it is pointing at us
			t=(p.dot(fNormal123)-fCdotN123)/vdotn; // #  distance to intersection
			if(t>=-fTol && t < tmin) {  // if not true, we're going away from this face or it's not close
				hp=p+vu*(t+extraDistance); // a little beyond point of intersection
				if (
					hp.dot(fNormal134)-fCdotN134 < 0.0 &&
					hp.dot(fNormal142)-fCdotN142 < 0.0 &&
					hp.dot(fNormal234)-fCdotN234 < 0.0 ) tmin=t;
			}
		}
		
		vdotn=-vu.dot(fNormal134);
		if(vdotn > 1e-12) { // # this is a candidate face, since it is pointing at us
			t=(p.dot(fNormal134)-fCdotN134)/vdotn; // #  distance to intersection
			if(t>=-fTol && t < tmin) { // if not true, we're going away from this face
				hp=p+vu*(t+extraDistance); // a little beyond point of intersection
				if (
					hp.dot(fNormal123)-fCdotN123 < 0.0 && 
					hp.dot(fNormal142)-fCdotN142 < 0.0 &&
					hp.dot(fNormal234)-fCdotN234 < 0.0 ) tmin=t;
			}
		}
				
		vdotn=-vu.dot(fNormal142);
		if(vdotn > 1e-12) { // # this is a candidate face, since it is pointing at us
			t=(p.dot(fNormal142)-fCdotN142)/vdotn; // #  distance to intersection
			if(t>=-fTol && t<tmin) { // if not true, we're going away from this face
				hp=p+vu*(t+extraDistance); // a little beyond point of intersection
				if (
					hp.dot(fNormal123)-fCdotN123 < 0.0 &&
					hp.dot(fNormal134)-fCdotN134 < 0.0 &&
					hp.dot(fNormal234)-fCdotN234 < 0.0 ) tmin=t;
			}
		}

		vdotn=-vu.dot(fNormal234);
		if(vdotn > 1e-12) { // # this is a candidate face, since it is pointing at us
			t=(p.dot(fNormal234)-fCdotN234)/vdotn; // #  distance to intersection
			if(t>=-fTol && t<tmin) { // if not true, we're going away from this face
				hp=p+vu*(t+extraDistance); // a little beyond point of intersection
				if (
					hp.dot(fNormal123)-fCdotN123 < 0.0 &&
					hp.dot(fNormal134)-fCdotN134 < 0.0 &&
					hp.dot(fNormal142)-fCdotN142 < 0.0 ) tmin=t;
			}
		}

	return fmax(0.0,tmin);
}

//////////////////////////////////////////////////////////////////////////
// 
// Approximate distance to tet.
// returns distance to sphere centered on bounding box
// - If inside return 0

G4double G4Tet::DistanceToIn(const G4ThreeVector& p) const
{
	G4double dd=(p-fMiddle).mag() - fMaxSize - fTol;
	return fmax(0.0, dd);
}

/////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface of box from inside
// by calculating distances to box's x/y/z planes.
// Smallest distance is exact distance to exiting.

G4double G4Tet::DistanceToOut( const G4ThreeVector& p,const G4ThreeVector& v,
                               const G4bool calcNorm,
							   G4bool *validNorm,G4ThreeVector *n) const
{
		G4ThreeVector vu(v.unit());
		G4double t1=kInfinity,t2=kInfinity,t3=kInfinity,t4=kInfinity, vdotn, tt;
				
		vdotn=vu.dot(fNormal123);
		if(vdotn > 1e-12)  // #we're heading towards this face, so it is a candidate
			t1=(fCdotN123-p.dot(fNormal123))/vdotn; // #  distance to intersection
		
		vdotn=vu.dot(fNormal134);
		if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
			t2=(fCdotN134-p.dot(fNormal134))/vdotn; // #  distance to intersection

		vdotn=vu.dot(fNormal142);
		if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
			t3=(fCdotN142-p.dot(fNormal142))/vdotn; // #  distance to intersection

		vdotn=vu.dot(fNormal234);
		if(vdotn > 1e-12) // #we're heading towards this face, so it is a candidate
			t4=(fCdotN234-p.dot(fNormal234))/vdotn; // #  distance to intersection
		
		tt=fmin(fmin(fmin(t1,t2),t3),t4);
		
		if (warningFlag && (tt == kInfinity || tt < -fTol)) {
			DumpInfo();
			G4cout << "p = " << p / mm << "mm" << G4endl;
			G4cout << "v = " << v  << G4endl;
			G4cout << "t1, t2, t3, t4 (mm) " << t1/mm << ", " << t2/mm << ", " << t3/mm << ", " << t4/mm << G4endl;
			G4cout << G4endl;
			G4Exception("G4Tet::DistanceToOut(p, v, ...)", "Notification", JustWarning, 
						 "No good intersection found or we are already outside!?" );
			if(validNorm) *validNorm=false; // flag normal as meaningless		
		} else if(calcNorm && n) {
			static G4ThreeVector normal;
			if(tt==t1) normal=fNormal123;
			else if (tt==t2) normal=fNormal134;
			else if (tt==t3) normal=fNormal142;
			else if (tt==t4) normal=fNormal234;
			n=&normal;
			if(validNorm) *validNorm=true; 		
		} 	

		return fmax(tt,0.0); // avoid tt<0.0 by a tiny bit if we are right on a face
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - If outside return 0

G4double G4Tet::DistanceToOut(const G4ThreeVector& p) const
{
	G4double t1,t2,t3,t4;
	t1=fCdotN123-p.dot(fNormal123); //  distance to plane, must be positive if we are inside
	t2=fCdotN134-p.dot(fNormal134); //  distance to plane
	t3=fCdotN142-p.dot(fNormal142); //  distance to plane
	t4=fCdotN234-p.dot(fNormal234); //  distance to plane
	
	// if any one of these is negative, we are outside, so return zero in that case
	G4double tmin=fmin(fmin(fmin(t1,t2),t3),t4);
	return (tmin < fTol)? 0:tmin;
}

////////////////////////////////////////////////////////////////////////
//
// Create a List containing the transformed vertices
// Note:
//  Caller has deletion resposibility

G4ThreeVectorList*
G4Tet::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
	G4ThreeVectorList* vertices = new G4ThreeVectorList();
	vertices->reserve(4);
	
	if (vertices)
	{
		G4ThreeVector vertex0(fAnchor) ;
		G4ThreeVector vertex1(fP2) ;
		G4ThreeVector vertex2(fP3) ;
		G4ThreeVector vertex3(fP4) ;
		
		vertices->push_back(pTransform.TransformPoint(vertex0));
		vertices->push_back(pTransform.TransformPoint(vertex1));
		vertices->push_back(pTransform.TransformPoint(vertex2));
		vertices->push_back(pTransform.TransformPoint(vertex3));
	}
	else
	{
		DumpInfo();
		G4Exception("G4Tet::CreateRotatedVertices()",
					"FatalError", FatalException,
					"Error in allocation of vertices. Out of memory !");
	}
	return vertices;
}

//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4Tet::GetEntityType() const
{
	return G4String("G4Tet");
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4Tet::StreamInfo(std::ostream& os) const
{
	os << "-----------------------------------------------------------\n"
	<< "    *** Dump for solid - " << GetName() << " ***\n"
	<< "    ===================================================\n"
	<< " Solid type: G4Tet\n"
	<< " Parameters: \n"
	<< "    anchor: " << fAnchor/mm << " mm \n"
	<< "    p2: " << fP2/mm << " mm \n"
	<< "    p3: " << fP3/mm << " mm \n"
	<< "    p4: " << fP4/mm << " mm \n"
	<< "    normal123: " << fNormal123 << " \n"
	<< "    normal134: " << fNormal134 << " \n"
	<< "    normal142: " << fNormal142 << " \n"
	<< "    normal234: " << fNormal234 << " \n"
	<< "-----------------------------------------------------------\n";
	
	return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4Tet::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
	scene.AddSolid (*this);
}

G4VisExtent G4Tet::GetExtent() const 
{
	return G4VisExtent (fXMin, fXMax, fYMin, fYMax, fZMin, fZMax);
}

G4Polyhedron* G4Tet::CreatePolyhedron () const 
{
	G4Polyhedron *ph=new G4Polyhedron;
	G4double xyz[4][3];
	static G4int faces[4][4]={{1,2,3,0},{1,3,4,0},{1,4,2,0},{2,3,4,0}};
	xyz[0][0]=fAnchor.x(); xyz[0][1]=fAnchor.y(); xyz[0][2]=fAnchor.z();
	xyz[1][0]=fP2.x(); xyz[1][1]=fP2.y(); xyz[1][2]=fP2.z();
	xyz[2][0]=fP3.x(); xyz[2][1]=fP3.y(); xyz[2][2]=fP3.z();
	xyz[3][0]=fP4.x(); xyz[3][1]=fP4.y(); xyz[3][2]=fP4.z();
	
	ph->createPolyhedron(4,4,xyz,faces);
	
	return ph;
}

G4NURBS* G4Tet::CreateNURBS () const 
{
	return new G4NURBSbox (fDx, fDy, fDz);
}
