// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Cons2.cc,v 1.1 1999-01-08 16:31:51 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Simple test of G4Cons
// Basic checks on each function + awkward cases for tracking / geom algorithms

#include <assert.h>
#include <math.h>
#include "G4ios.hh"

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Cons.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#define	DELTA 0.0001

// Returns false if actual is within wanted+/- DELTA
//         true if error
G4bool OutRange(G4double actual,G4double wanted)
{
    G4bool rng=false;
    if (actual<wanted-DELTA||actual>wanted+DELTA) rng=true;
    return rng;
}
G4bool OutRange(G4ThreeVector actual,G4ThreeVector wanted)
{
    G4bool rng=false;
    if (OutRange(actual.x(),wanted.x())
	||OutRange(actual.y(),wanted.y())
	||OutRange(actual.z(),wanted.z())  ) rng=true;
    return rng;
}

int main(void)
{
	double Dist,low,high;

	G4ThreeVector   pzero(0,0,0);
	
	G4ThreeVector   pplx(120,0,0),pply(0,120,0),pplz(0,0,120);
	
	G4ThreeVector   pmix(-120,0,0),pmiy(0,-120,0),pmiz(0,0,-120);
	
	G4ThreeVector   ponmiz(0,75,-50),ponplz(0,75,50);
	
	G4ThreeVector	ponr1(sqrt(50*50/2.0),sqrt(50*50/2.0),0),
	                ponr2(sqrt(100*100/2.0),sqrt(100*100/2.0),0),
	        	ponphi1(60*cos(M_PI/6),-60*sin(M_PI/6),0),
		        ponphi2(60*cos(M_PI/6),60*sin(M_PI/6),0),
	                ponr2b(150,0,0);
	
	G4ThreeVector pnearplz(45,45,45),pnearmiz(45,45,-45);
	G4ThreeVector pydx(60,150,0),pbigx(500,0,0);

	G4ThreeVector proot1(0,125,-1000),proot2(0,75,-1000);
	
	G4ThreeVector pparr1(0,25,-150);                // Test case parallel to both rs of c8
	G4ThreeVector pparr2(0,75,-50),pparr3(0,125,50);
	G4ThreeVector vparr(0,1./sqrt(5.),2./sqrt(5.)); 

	G4ThreeVector vnphi1(-sin(M_PI/6),-cos(M_PI/6),0),
	              vnphi2(-sin(M_PI/6),cos(M_PI/6),0);

	G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1),
	              vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1),
		      vxy(1./sqrt(2.),1./sqrt(2.),0),
	              vxmy(1./sqrt(2.),-1./sqrt(2.),0),
	              vmxmy(-1./sqrt(2.),-1./sqrt(2.),0),
	              vmxy(-1./sqrt(2.),1./sqrt(2.),0),
	              vxmz(1./sqrt(2.),0,-1./sqrt(2.));
	
	G4RotationMatrix r90X,r90Y,r90Z,r180X,r45X,r30Y;
	
	G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*M_PI),
	       c2("Hollow Full Cone",50,100,50,200,50,-1,2*M_PI),
	       c3("Hollow Cut Tube",50,100,50,100,50,-M_PI/6,M_PI/3),
	       c4("Hollow Cut Cone",50,100,50,200,50,-M_PI/6,M_PI/3),
	       c5("Hollow Cut Cone",25,50,75,150,50,0,3*M_PI/2),
 	       c6("Solid Full Cone",0,150,0,150,50,0,2*M_PI),
	       c7("Thin Tube",95,100,95,100,50,0,2*M_PI),
	       c8a("Solid Full Cone2",0,100,0,150,50,0,2*M_PI),
	       c8b("Hollow Full Cone2",50,100,100,150,50,0,2*M_PI),
	       c8c("Hollow Full Cone2inv",100,150,50,100,50,0,2*M_PI);
	
	G4ThreeVector norm,*pNorm;
	G4bool *pgoodNorm,goodNorm,calcNorm=true;
	
	pNorm=&norm;
	pgoodNorm=&goodNorm;

	r90X.rotateX(M_PI/2);
	r90Y.rotateY(M_PI/2);
	r90Z.rotateZ(M_PI/2);
	r45X.rotateX(M_PI/4);
	r30Y.rotateY(M_PI/6);
	//	G4cout << "G4Cons:"<< c4.GetName()
	//   << " ID=" << c4.GetIdentifier() << "\n";

	G4cout << "Testing G4Cons::Inside...\n";
	if (c1.Inside(pzero)!=kOutside)
		G4cout << "Error A" << endl;
	if (c6.Inside(pzero)!=kInside)
		G4cout << "Error A2" << endl;
	if (c1.Inside(pplx)!=kOutside)
	    G4cout << "Error B1" << endl;
	if (c2.Inside(pplx)!=kInside)
	    G4cout << "Error B2" << endl;
	if (c3.Inside(pplx)!=kOutside)
	    G4cout << "Error B3" << endl;
	if (c4.Inside(pplx)!=kInside)
	    G4cout << "Error B4" << endl;
	if (c1.Inside(ponmiz)!=kSurface)
	    G4cout << "Error C" << endl;
	if (c1.Inside(ponplz)!=kSurface)
	    G4cout << "Error D" << endl;
	if (c1.Inside(ponr1)!=kSurface)
	    G4cout << "Error E" << endl;
	if (c1.Inside(ponr2)!=kSurface)
	    G4cout << "Error F" << endl;
	if (c3.Inside(ponphi1)!=kSurface)
	    G4cout << "Error G" << endl;
	if (c3.Inside(ponphi2)!=kSurface)
	    G4cout << "Error H" << endl;

	if (c5.Inside(G4ThreeVector(70,1,0))!=kInside)
	    G4cout << "Error I" << endl;
	if (c5.Inside(G4ThreeVector(50,-50,0))!=kOutside)
	    G4cout << "Error I2" << endl;
	if (c5.Inside(G4ThreeVector(70,0,0))!=kSurface)
	    G4cout << "Error I3" << endl;
// on tolerant r, inside z, within phi
	if (c5.Inside(G4ThreeVector(100,0,0))!=kSurface)
	    G4cout << "Error I4" << endl;
	if (c3.Inside(G4ThreeVector(100,0,0))!=kSurface)
	    G4cout << "Error I5" << endl;
// on tolerant r, tolerant z, within phi
	if (c5.Inside(G4ThreeVector(100,0,50))!=kSurface)
	    G4cout << "Error I4" << endl;
	if (c3.Inside(G4ThreeVector(100,0,50))!=kSurface)
	    G4cout << "Error I5" << endl;
	

	G4cout << "Testing G4Cons::SurfaceNormal...\n";
	norm=c1.SurfaceNormal(ponplz);
	if (OutRange(norm,G4ThreeVector(0,0,1)))
	    G4cout << "Error A " << norm << endl;
	norm=c1.SurfaceNormal(ponmiz);
	if (OutRange(norm,G4ThreeVector(0,0,-1)))
	    G4cout << "Error B " << norm << endl;
	norm=c1.SurfaceNormal(ponr1);
	if (OutRange(norm,G4ThreeVector(-1.0/sqrt(2.0),-1.0/sqrt(2.0),0)))
	    G4cout << "Error C " << norm << endl;
	norm=c1.SurfaceNormal(ponr2);
	if (OutRange(norm,G4ThreeVector(1.0/sqrt(2.0),1.0/sqrt(2.0),0)))
	    G4cout << "Error D " << norm << endl;
	norm=c3.SurfaceNormal(ponphi1);
	if (OutRange(norm,vnphi1))
	    G4cout << "Error E " << norm << endl;
	norm=c3.SurfaceNormal(ponphi2);
	if (OutRange(norm,vnphi2))
	    G4cout << "Error F " << norm << endl;
	norm=c4.SurfaceNormal(ponr2b);
	if (OutRange(norm,vxmz))
	    G4cout << "Error G " << norm << endl;

	norm=c5.SurfaceNormal(G4ThreeVector(51,0,-50));
	if (OutRange(norm,vmz))
	    G4cout << "Errot H " << norm << endl;

	G4cout << "Testing G4Cons::DistanceToOut...\n";
	Dist=c4.DistanceToOut(ponphi1);
	if (OutRange(Dist,0))
		G4cout << "Error A " << Dist << endl;

	Dist=c1.DistanceToOut(ponphi1);
	if (OutRange(Dist,10))
		G4cout << "Error B " << Dist << endl;

	Dist=c1.DistanceToOut(pnearplz);
	if (OutRange(Dist,5))
		G4cout << "Error C " << Dist << endl;
	Dist=c1.DistanceToOut(pnearmiz);
	if (OutRange(Dist,5))
		G4cout << "Error D " << Dist << endl;

	Dist=c1.DistanceToOut(ponr1);
	if (OutRange(Dist,0))
	    G4cout << "Error E " << Dist << endl;
	Dist=c1.DistanceToOut(ponr2);
	if (OutRange(Dist,0))
	    G4cout << "Error F " << Dist << endl;

	Dist=c6.DistanceToOut(pzero);
	if (OutRange(Dist,50))
	    G4cout << "Error G " << Dist << endl;

	Dist=c5.DistanceToOut(G4ThreeVector(0,-70,0));
	if (OutRange(Dist,0))
	    G4cout << "Error H " << Dist << endl;
	
       	G4cout << "Testing G4Cons::DistanceToOut...\n";
	Dist=c4.DistanceToOut(pplx,vx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,30)||OutRange(*pNorm,vxmz)||!*pgoodNorm)
	    G4cout << "Error Rmax1 " << Dist << endl;

	Dist=c2.DistanceToOut(pplx,vx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,30)||OutRange(*pNorm,vxmz)||!*pgoodNorm)
	    G4cout << "Error Rmax2 " << Dist << endl;

	Dist=c4.DistanceToOut(pplx,vmx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,70)||*pgoodNorm)
	    G4cout << "Error Rmin1 " << Dist << endl;


	Dist=c2.DistanceToOut(pplx,vmx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,70)||*pgoodNorm)
	    G4cout << "Error Rmin2 " << Dist << endl;

	Dist=c3.DistanceToOut(ponphi1,vmy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    OutRange(*pNorm,vnphi1)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 1" << Dist << endl;
	Dist=c3.DistanceToOut(ponphi1,vy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,2*60*sin(M_PI/6))||
	    OutRange(*pNorm,vnphi2)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 2" << Dist << endl;

	Dist=c3.DistanceToOut(ponphi2,vy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    OutRange(*pNorm,vnphi2)||
	    !*pgoodNorm)
	    G4cout << "Error PhiE 1" << Dist << endl;
	Dist=c3.DistanceToOut(ponphi2,vmy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,2*60*sin(M_PI/6))||
	    OutRange(*pNorm,vnphi1)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 2" << Dist << endl;


	Dist=c6.DistanceToOut(ponplz,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100)||
	    OutRange(*pNorm,vmz)||
	    !*pgoodNorm)
	    G4cout << "Error Top Z 1" << Dist << endl;
	Dist=c6.DistanceToOut(ponplz,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    OutRange(*pNorm,vz)||
	    !*pgoodNorm)
	    G4cout << "Error Top Z 2" << Dist << endl;

	Dist=c6.DistanceToOut(ponmiz,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100)||
	    OutRange(*pNorm,vz)||
	    !*pgoodNorm)
	    G4cout << "Error Lower Z 1" << Dist << endl;
	Dist=c6.DistanceToOut(ponmiz,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    OutRange(*pNorm,vmz)||
	    !*pgoodNorm)
	    G4cout << "Error Lower Z 2" << Dist << endl;

// Test case for rmax root bug
	Dist=c7.DistanceToOut(ponr2,vmx,calcNorm,pgoodNorm,pNorm);
	if (OutRange(Dist,100/sqrt(2.)-sqrt(95*95-100*100/2.))||*pgoodNorm)
	    G4cout << "Error rmax root bug" << Dist << endl;

// Parallel radii test cases
	Dist=c8a.DistanceToOut(pparr2,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100.*sqrt(5.)/2.)||
                     !*pgoodNorm||
                     OutRange(*pNorm,vz))
	    G4cout << "Error solid parr2a " <<Dist << endl;
	Dist=c8a.DistanceToOut(pparr2,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr2b " <<Dist << endl;

	Dist=c8a.DistanceToOut(pparr2,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr2c " <<Dist << endl;
	Dist=c8a.DistanceToOut(pparr2,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr2d " <<Dist << endl;

	Dist=c8a.DistanceToOut(pparr3,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr3a " <<Dist << endl;
	Dist=c8a.DistanceToOut(pparr3,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100*sqrt(5.)/2.)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr3b " <<Dist << endl;
	Dist=c8a.DistanceToOut(pparr3,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr3c " <<Dist << endl;

	Dist=c8a.DistanceToOut(pparr3,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,50)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,2./sqrt(5.0),-1./sqrt(5.0))))
	    G4cout << "Error solid parr3d " <<Dist << endl;


	Dist=c8b.DistanceToOut(pparr2,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100*sqrt(5.)/2.)||
                     !*pgoodNorm||
                     OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr2a " <<Dist << endl;
	Dist=c8b.DistanceToOut(pparr2,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr2b " <<Dist << endl;

	Dist=c8b.DistanceToOut(pparr2,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,50)||*pgoodNorm)
	    G4cout << "Error hollow parr2c " <<Dist << endl;
	Dist=c8b.DistanceToOut(pparr2,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr2d " <<Dist << endl;


	Dist=c8b.DistanceToOut(pparr3,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr3a " <<Dist << endl;
	Dist=c8b.DistanceToOut(pparr3,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,100.*sqrt(5.)/2.)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr3b " <<Dist << endl;
	Dist=c8b.DistanceToOut(pparr3,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr3c " <<Dist << endl;

	Dist=c8b.DistanceToOut(pparr3,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(Dist,50)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,2./sqrt(5.),-1.0/sqrt(5.))))
	    G4cout << "Error hollow parr3d " <<Dist << endl;

	G4cout << "Testing G4Cons::DistanceToIn...\n";
	Dist=c1.DistanceToIn(pzero);
	if (OutRange(Dist,50))
	  G4cout << "Error A " << Dist << endl;

	Dist=c1.DistanceToIn(pplx);
	if (OutRange(Dist,20))
	  G4cout << "Error B " << Dist << endl;

	Dist=c1.DistanceToIn(pply);
	if (OutRange(Dist,20))
	  G4cout << "Error C " << Dist << endl;

	Dist=c4.DistanceToIn(pply);
	if (OutRange(Dist,120*sin(M_PI/3)))
	  G4cout << "Error D " << Dist << endl;

	Dist=c4.DistanceToIn(pmiy);
	if (OutRange(Dist,120*sin(M_PI/3)))
	  G4cout << "Error D " << Dist << endl;

	Dist=c1.DistanceToIn(pplz);
	if (OutRange(Dist,70))
	    G4cout << "Error E " << Dist << endl;
// Check with both rmins=0
	Dist=c5.DistanceToIn(pplx);
	if (OutRange(Dist,20./sqrt(2.)))
	  G4cout << "Error F " << Dist << endl;

	G4cout << "Testing G4Cons::DistanceToIn...\n";

	Dist=c1.DistanceToIn(pplz,vmz);
	if (OutRange(Dist,kInfinity))
	  G4cout << "Error A " << Dist << endl;
	Dist=c1.DistanceToIn(pplx,vmx);
	if (OutRange(Dist,20))
	  G4cout << "Error B " << Dist << endl;
	Dist=c1.DistanceToIn(pplz,vx);
	if (OutRange(Dist,kInfinity))
	  G4cout << "Error C " << Dist << endl;
	Dist=c4.DistanceToIn(pply,vmy);
	if (OutRange(Dist,kInfinity))
	  G4cout << "Error D " << Dist << endl;

	Dist=c1.DistanceToIn(pydx,vmy);
	if (OutRange(Dist,70))
	  G4cout << "Error E " << Dist << endl;
	Dist=c3.DistanceToIn(pydx,vmy);
	if (OutRange(Dist,150-60*tan(M_PI/6)))
	  G4cout << "Error F " << Dist << endl;

	Dist=c1.DistanceToIn(pplx,vmx);
	if (OutRange(Dist,20))
	  G4cout << "Error G " << Dist << endl;
	Dist=c1.DistanceToIn(pplx,vx);
	if (OutRange(Dist,kInfinity))
	  G4cout << "Error G2 " << Dist << endl;

	Dist=c4.DistanceToIn(pbigx,vmx);
	if (OutRange(Dist,350))
	    G4cout << "Error G3 " << Dist << endl;

	Dist=c4.DistanceToIn(pzero,vx);
	if (OutRange(Dist,50))
	  G4cout << "Error H " << Dist << endl;

	Dist=c1.DistanceToIn(ponr2,vx);
	if (OutRange(Dist,kInfinity))
	    G4cout << "Error I" << Dist << endl;
	Dist=c1.DistanceToIn(ponr2,vmx);
	if (OutRange(Dist,0))
	    G4cout << "Error I2" << Dist << endl;
	
	Dist=c1.DistanceToIn(ponr1,vx);
	if (OutRange(Dist,0))
	    G4cout << "Error J" << Dist << endl;
	Dist=c1.DistanceToIn(ponr1,vmx);
	if (OutRange(Dist,2.0*sqrt(50*50/2.)))
	    G4cout << "Error J2" << Dist << endl;

	Dist=c1.DistanceToIn(ponr2,vmxmy);
	if (OutRange(Dist,0))
	    G4cout << "Error K" << Dist << endl;

// Parallel test case -> parallel to both radii
	Dist=c8b.DistanceToIn(pparr1,vparr);
	if (OutRange(Dist,100*sqrt(5.)/2.))
	    G4cout << "Error parr1 " << Dist << endl;
	Dist=c8b.DistanceToIn(pparr2,-vparr);
	if (OutRange(Dist,kInfinity))
	    G4cout << "Error parr2 " << Dist << endl;
	Dist=c8b.DistanceToIn(pparr3,vparr);
	if (OutRange(Dist,kInfinity))
	    G4cout << "Error parr3a " << Dist << endl;
	Dist=c8b.DistanceToIn(pparr3,-vparr);
	if (OutRange(Dist,0))
	    G4cout << "Error parr3b " << Dist << endl;

// Check we don't Hit `shadow cone' at `-ve radius' on rmax or rmin
	Dist=c8a.DistanceToIn(proot1,vz);
	if (OutRange(Dist,1000))
	    G4cout << "Error shadow rmax root problem " << Dist << endl;

	Dist=c8c.DistanceToIn(proot2,vz);
	if (OutRange(Dist,1000))
	    G4cout << "Error shadow rmin root problem " << Dist << endl;

	G4cout << "NOT Checking G4Cons::ScopeCar...\n";
	G4cout << "NOT Checking G4Cons::ScopePhi...\n";
	G4cout << "NOT Checking G4Cons::ScopeRad...\n";

	return 0;
}











