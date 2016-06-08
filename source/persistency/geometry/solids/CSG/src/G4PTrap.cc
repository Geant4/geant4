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
// $Id: G4PTrap.cc,v 1.3 2001/07/11 10:02:24 gunter Exp $
// GEANT4 tag $Name: geant4-04-01 $
//
// class G4PTrap
//
// Implementation for G4PTrap class
//
// History:
// 19.06.98 A.Kimura Converted G4Trap.cc

#include <math.h>

#include "G4VSolid.hh"
#include "G4PTrap.hh"
#include "G4Trap.hh"

const G4double kCoplanar_Tolerance=1E-4;

// Destructor

G4PTrap::~G4PTrap()
{;}

// make a transient object
G4VSolid* G4PTrap::MakeTransientObject() const {
    G4VSolid* transientObject = new G4Trap(GetName(),
					 fDx1, fDx2,
					 fDy1, fDy2,
					 fDz);
    return transientObject;
}

// Constructor for G4Trd

G4PTrap::G4PTrap(const G4Trap* theTrap)
 : G4PCSGSolid(theTrap->GetName())
{

    G4double pDx1 = theTrap->GetXHalfLength1();
    G4double pDx2 = theTrap->GetXHalfLength2();
    G4double pDy1 = theTrap->GetYHalfLength1();
    G4double pDy2 = theTrap->GetYHalfLength2();
    G4double pDz  = theTrap->GetZHalfLength();

    G4bool good;

    if (pDz>0 && pDy1>0 && pDx1>0 && pDx2>0 && pDy2>0 )
	{
	   fDz = pDz;
	   fTthetaCphi = 0 ;
	   fTthetaSphi = 0 ;
	    
	   fDy1 = pDy1 ;
	   fDx1 = pDx1 ;
	   fDx2 = pDx1 ;
	   fTalpha1 = 0 ;
	   
	   fDy2 = pDy2 ;
	   fDx3 = pDx2 ;
	   fDx4 = pDx2 ;
	   fTalpha2 = 0 ;

	   G4ThreeVector pt[8] ;
	   
	   pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
           pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
	   pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
	   pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
	   pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
	   pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
// Bottom side with normal approx. -Y
	    good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
	    if (!good)
		{
		    G4Exception("G4PTrap::G4PTrap - face at ~-Y not planar");

		}
// Top side with normal approx. +Y
	    good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
	    if (!good)
		{
		    G4Exception("G4PTrap::G4PTrap - face at ~+Y not planar");
		}
// Front side with normal approx. -X
	    good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
	    if (!good)
		{
		    G4Exception("G4PTrap::G4PTrap - face at ~-X not planar");
		}
// Back side iwth normal approx. +X
	    good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
	    if (!good)
		{
		    G4Exception("G4PTrap::G4PTrap - face at ~+X not planar");
		}
	}
    else
	{
	    G4Exception("Error in G4PTrap::G4PTrap - Invalid Length G4PTrapmeters");
	}
}

// Set all parameters, as for constructor - check and set half-widths
// as well as angles: final check of coplanarity

void G4PTrap::SetAllParameters ( G4double pDz,
				G4double pTheta,
				G4double pPhi,
				G4double pDy1,
				G4double pDx1,
				G4double pDx2,
				G4double pAlp1,
				G4double pDy2,
				G4double pDx3,
				G4double pDx4,
				G4double pAlp2)
{
    if (pDz>0 && pDy1>0 && pDx1>0 && pDx2>0 && pDy2>0 && pDx3>0 && pDx4>0)
	{
	   fDz=pDz;
	   fTthetaCphi=tan(pTheta)*cos(pPhi);
	   fTthetaSphi=tan(pTheta)*sin(pPhi);
	    
	   fDy1=pDy1;
	   fDx1=pDx1;
	   fDx2=pDx2;
	   fTalpha1=tan(pAlp1);
	   
	   fDy2=pDy2;
	   fDx3=pDx3;
	   fDx4=pDx4;
	   fTalpha2=tan(pAlp2);

	   MakePlanes();
	}
    else
	{
	    G4Exception("Error in G4PTrap::SetAll - Invalid Length G4PTrapmeters");
	}
}

G4bool G4PTrap::MakePlanes()
{
  G4bool good = true;

  G4ThreeVector pt[8] ;
	   
  pt[0]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=G4ThreeVector(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,-fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=G4ThreeVector(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,-fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=G4ThreeVector(+fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,+fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=G4ThreeVector(+fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,+fDz*fTthetaSphi+fDy2,+fDz);
  // Bottom side with normal approx. -Y
  good=MakePlane(pt[0],pt[4],pt[5],pt[1],fPlanes[0]);
  if (!good)
    {
      G4Exception("G4PTrap::G4PTrap - face at ~-Y not planar");
    }
  // Top side with normal approx. +Y
  good=MakePlane(pt[2],pt[3],pt[7],pt[6],fPlanes[1]);
  if (!good)
    {
      G4Exception("G4PTrap::G4PTrap - face at ~+Y not planar");
    }
  // Front side with normal approx. -X
  good=MakePlane(pt[0],pt[2],pt[6],pt[4],fPlanes[2]);
  if (!good)
    {
      G4Exception("G4PTrap::G4PTrap - face at ~-X not planar");
    }
  // Back side iwth normal approx. +X
  good=MakePlane(pt[1],pt[5],pt[7],pt[3],fPlanes[3]);
  if (!good)
    {
      G4Exception("G4PTrap::G4PTrap - face at ~+X not planar");
    }
  return good;
}

// --------------------------------------------------------------------------

// Calculate the coef's of the plane p1->p2->p3->p4->p1
// where the ThreeVectors 1-4 are in anti-clockwise order when viewed from
// infront of the plane.
//
// Return true if the ThreeVectors are coplanar + set coef;s
//        false if ThreeVectors are not coplanar

G4bool G4PTrap::MakePlane( const G4ThreeVector& p1,
			  const G4ThreeVector& p2,
		          const G4ThreeVector& p3,
			  const G4ThreeVector& p4,
			  PTrapSidePlane& plane )
{
    G4double a,b,c,s;
    G4ThreeVector v12,v13,v14,Vcross;

    G4bool good;

    v12 = p2-p1;
    v13 = p3-p1;
    v14 = p4-p1;
    Vcross=v12.cross(v13);

    if (fabs(Vcross.dot(v14)/(Vcross.mag()*v14.mag()))>kCoplanar_Tolerance)
	{
	    good=false;
	}
    else
	{
    
// a,b,c correspond to the x/y/z components of the normal vector to the plane
	   
           // a=(p2.y()-p1.y())*(p1.z()+p2.z())+(p3.y()-p2.y())*(p2.z()+p3.z());
	   // a+=(p4.y()-p3.y())*(p3.z()+p4.z())+(p1.y()-p4.y())*(p4.z()+p1.z()); // may be delete ?
    
	   // b=(p2.z()-p1.z())*(p1.x()+p2.x())+(p3.z()-p2.z())*(p2.x()+p3.x());
	   //  b+=(p4.z()-p3.z())*(p3.x()+p4.x())+(p1.z()-p4.z())*(p4.x()+p1.x()); // ?
	    
	   // c=(p2.x()-p1.x())*(p1.y()+p2.y())+(p3.x()-p2.x())*(p2.y()+p3.y());
	   // c+=(p4.x()-p3.x())*(p3.y()+p4.y())+(p1.x()-p4.x())*(p4.y()+p1.y());  // ?
// Let create diagonals 4-2 and 3-1 than (4-2)x(3-1) provides vector perpendicular to the
// plane directed to outside !!! and a,b,c, = f(1,2,3,4)	   
	   a = +(p4.y() - p2.y())*(p3.z() - p1.z()) - (p3.y() - p1.y())*(p4.z() - p2.z()) ;
	   b = -(p4.x() - p2.x())*(p3.z() - p1.z()) + (p3.x() - p1.x())*(p4.z() - p2.z()) ; 
	   c = +(p4.x() - p2.x())*(p3.y() - p1.y()) - (p3.x() - p1.x())*(p4.y() - p2.y()) ;
	   s=sqrt(a*a+b*b+c*c);   // so now vector plane.(a,b,c) is unit 
	   plane.a=a/s;
	   plane.b=b/s;
	   plane.c=c/s;
                          // Calculate D: p1 in in plane so D=-n.p1.Vect()
	    
	    plane.d=-(plane.a*p1.x()+plane.b*p1.y()+plane.c*p1.z());
	    good=true;
	}
    return good;
}


// ********************************  End of G4PTrap.cc   ********************************





