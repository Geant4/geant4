// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Cons2.cc,v 1.5 2000-08-16 08:02:36 grichine Exp $
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
    G4bool rng = false ;
    if (actual < wanted-DELTA || actual > wanted + DELTA ) rng = true ;
    return rng ;
}
G4bool OutRange(G4ThreeVector actual,G4ThreeVector wanted)
{
    G4bool rng = false ;
    if (OutRange(actual.x(),wanted.x())
	||OutRange(actual.y(),wanted.y())
	||OutRange(actual.z(),wanted.z())  ) rng = true ;
    return rng ;
}

int main(void)
{
	double dist,low,high;

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
	
	G4ThreeVector pparr1(0,25,-150);   // Test case parallel to both rs of c8
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
                vx2mz(1.0/sqrt(5.0),0,-2.0/sqrt(5.0)),
	        vxmz(1./sqrt(2.),0,-1./sqrt(2.));
	
	G4RotationMatrix r90X,r90Y,r90Z,r180X,r45X,r30Y;
	
  G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*M_PI),
	 c2("Hollow Full Cone",50,100,50,200,50,-1,2*M_PI),
	 c3("Hollow Cut Tube",50,100,50,100,50,-M_PI/6,M_PI/3),
	 c4("Hollow Cut Cone",50,100,50,200,50,-M_PI/6,M_PI/3),
	 c5("Hollow Cut Cone",25,50,75,150,50,0,3*M_PI/2),
 	 c6("Solid Full Tube",0,150,0,150,50,0,2*M_PI),
	 c7("Thin Tube",95,100,95,100,50,0,2*M_PI),
	 c8a("Solid Full Cone2",0,100,0,150,50,0,2*M_PI),
	 c8b("Hollow Full Cone2",50,100,100,150,50,0,2*M_PI),
	 c8c("Hollow Full Cone2inv",100,150,50,100,50,0,2*M_PI),
	 c9("Excotic Cone",50,60,
	    0,           // 1.0e-7,   500*kRadTolerance,
                           10,50,0,2*M_PI), 
	 cms("cms cone",0.0,70.0,0.0,157.8,2949.0,0.0,6.2831853071796);

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

        G4cout.precision(16) ;
	G4cout << "Testing G4Cons::Inside...\n";
	if (c1.Inside(pzero)!=kOutside)
		G4cout << "Error A" << G4endl;
	if (c6.Inside(pzero)!=kInside)
		G4cout << "Error A2" << G4endl;
	if (c1.Inside(pplx)!=kOutside)
	    G4cout << "Error B1" << G4endl;
	if (c2.Inside(pplx)!=kInside)
	    G4cout << "Error B2" << G4endl;
	if (c3.Inside(pplx)!=kOutside)
	    G4cout << "Error B3" << G4endl;
	if (c4.Inside(pplx)!=kInside)
	    G4cout << "Error B4" << G4endl;
	if (c1.Inside(ponmiz)!=kSurface)
	    G4cout << "Error C" << G4endl;
	if (c1.Inside(ponplz)!=kSurface)
	    G4cout << "Error D" << G4endl;
	if (c1.Inside(ponr1)!=kSurface)
	    G4cout << "Error E" << G4endl;
	if (c1.Inside(ponr2)!=kSurface)
	    G4cout << "Error F" << G4endl;
	if (c3.Inside(ponphi1)!=kSurface)
	    G4cout << "Error G" << G4endl;
	if (c3.Inside(ponphi2)!=kSurface)
	    G4cout << "Error H" << G4endl;

	if (c5.Inside(G4ThreeVector(70,1,0))!=kInside)
	    G4cout << "Error I" << G4endl;
	if (c5.Inside(G4ThreeVector(50,-50,0))!=kOutside)
	    G4cout << "Error I2" << G4endl;
	if (c5.Inside(G4ThreeVector(70,0,0))!=kSurface)
	    G4cout << "Error I3" << G4endl;
// on tolerant r, inside z, within phi
	if (c5.Inside(G4ThreeVector(100,0,0))!=kSurface)
	    G4cout << "Error I4" << G4endl;
	if (c3.Inside(G4ThreeVector(100,0,0))!=kSurface)
	    G4cout << "Error I5" << G4endl;
// on tolerant r, tolerant z, within phi
	if (c5.Inside(G4ThreeVector(100,0,50))!=kSurface)
	    G4cout << "Error I4" << G4endl;
	if (c3.Inside(G4ThreeVector(100,0,50))!=kSurface)
	    G4cout << "Error I5" << G4endl;
	

	G4cout << "Testing G4Cons::SurfaceNormal...\n";
	norm=c1.SurfaceNormal(ponplz);
	if (OutRange(norm,G4ThreeVector(0,0,1)))
	    G4cout << "Error A " << norm << G4endl;
	norm=c1.SurfaceNormal(ponmiz);
	if (OutRange(norm,G4ThreeVector(0,0,-1)))
	    G4cout << "Error B " << norm << G4endl;
	norm=c1.SurfaceNormal(ponr1);
	if (OutRange(norm,G4ThreeVector(-1.0/sqrt(2.0),-1.0/sqrt(2.0),0)))
	    G4cout << "Error C " << norm << G4endl;
	norm=c1.SurfaceNormal(ponr2);
	if (OutRange(norm,G4ThreeVector(1.0/sqrt(2.0),1.0/sqrt(2.0),0)))
	    G4cout << "Error D " << norm << G4endl;
	norm=c3.SurfaceNormal(ponphi1);
	if (OutRange(norm,vnphi1))
	    G4cout << "Error E " << norm << G4endl;
	norm=c3.SurfaceNormal(ponphi2);
	if (OutRange(norm,vnphi2))
	    G4cout << "Error F " << norm << G4endl;
	norm=c4.SurfaceNormal(ponr2b);
	if (OutRange(norm,vxmz))
	    G4cout << "Error G " << norm << G4endl;

	norm=c5.SurfaceNormal(G4ThreeVector(51,0,-50));
	if (OutRange(norm,vmz))
	    G4cout << "Errot H " << norm << G4endl;

	G4cout << "Testing G4Cons::DistanceToOut...\n";
	dist=c4.DistanceToOut(ponphi1);
	if (OutRange(dist,0))
		G4cout << "Error A " << dist << G4endl;

	dist=c1.DistanceToOut(ponphi1);
	if (OutRange(dist,10))
		G4cout << "Error B " << dist << G4endl;

	dist=c1.DistanceToOut(pnearplz);
	if (OutRange(dist,5))
		G4cout << "Error C " << dist << G4endl;
	dist=c1.DistanceToOut(pnearmiz);
	if (OutRange(dist,5))
		G4cout << "Error D " << dist << G4endl;

	dist=c1.DistanceToOut(ponr1);
	if (OutRange(dist,0))
	    G4cout << "Error E " << dist << G4endl;
	dist=c1.DistanceToOut(ponr2);
	if (OutRange(dist,0))
	    G4cout << "Error F " << dist << G4endl;

	dist=c6.DistanceToOut(pzero);
	if (OutRange(dist,50))
	    G4cout << "Error G " << dist << G4endl;

	dist=c5.DistanceToOut(G4ThreeVector(0,-70,0));
	if (OutRange(dist,0))
	    G4cout << "Error H " << dist << G4endl;
	
       	G4cout << "Testing G4Cons::DistanceToOut...\n";
	dist=c4.DistanceToOut(pplx,vx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,30)||OutRange(*pNorm,vxmz)||!*pgoodNorm)
	    G4cout << "Error Rmax1 " << dist << G4endl;

	dist=c2.DistanceToOut(pplx,vx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,30)||OutRange(*pNorm,vxmz)||!*pgoodNorm)
	    G4cout << "Error Rmax2 " << dist << G4endl;

	dist=c4.DistanceToOut(pplx,vmx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,70)||*pgoodNorm)
	    G4cout << "Error Rmin1 " << dist << G4endl;


	dist=c2.DistanceToOut(pplx,vmx,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,70)||*pgoodNorm)
	    G4cout << "Error Rmin2 " << dist << G4endl;

	dist=c3.DistanceToOut(ponphi1,vmy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    OutRange(*pNorm,vnphi1)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 1" << dist << G4endl;
	dist=c3.DistanceToOut(ponphi1,vy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,2*60*sin(M_PI/6))||
	    OutRange(*pNorm,vnphi2)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 2" << dist << G4endl;

	dist=c3.DistanceToOut(ponphi2,vy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    OutRange(*pNorm,vnphi2)||
	    !*pgoodNorm)
	    G4cout << "Error PhiE 1" << dist << G4endl;
	dist=c3.DistanceToOut(ponphi2,vmy,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,2*60*sin(M_PI/6))||
	    OutRange(*pNorm,vnphi1)||
	    !*pgoodNorm)
	    G4cout << "Error PhiS 2" << dist << G4endl;


	dist=c6.DistanceToOut(ponplz,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100)||
	    OutRange(*pNorm,vmz)||
	    !*pgoodNorm)
	    G4cout << "Error Top Z 1" << dist << G4endl;
	dist=c6.DistanceToOut(ponplz,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    OutRange(*pNorm,vz)||
	    !*pgoodNorm)
	    G4cout << "Error Top Z 2" << dist << G4endl;

	dist=c6.DistanceToOut(ponmiz,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100)||
	    OutRange(*pNorm,vz)||
	    !*pgoodNorm)
	    G4cout << "Error Lower Z 1" << dist << G4endl;
	dist=c6.DistanceToOut(ponmiz,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    OutRange(*pNorm,vmz)||
	    !*pgoodNorm)
	    G4cout << "Error Lower Z 2" << dist << G4endl;

// Test case for rmax root bug
	dist=c7.DistanceToOut(ponr2,vmx,calcNorm,pgoodNorm,pNorm);
	if (OutRange(dist,100/sqrt(2.)-sqrt(95*95-100*100/2.))||*pgoodNorm)
	    G4cout << "Error rmax root bug" << dist << G4endl;

// Parallel radii test cases
	dist=c8a.DistanceToOut(pparr2,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100.*sqrt(5.)/2.)||
                     !*pgoodNorm||
                     OutRange(*pNorm,vz))
	    G4cout << "Error solid parr2a " <<dist << G4endl;
	dist=c8a.DistanceToOut(pparr2,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr2b " <<dist << G4endl;

	dist=c8a.DistanceToOut(pparr2,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr2c " <<dist << G4endl;
	dist=c8a.DistanceToOut(pparr2,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr2d " <<dist << G4endl;

	dist=c8a.DistanceToOut(pparr3,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr3a " <<dist << G4endl;
	dist=c8a.DistanceToOut(pparr3,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100*sqrt(5.)/2.)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error solid parr3b " <<dist << G4endl;
	dist=c8a.DistanceToOut(pparr3,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error solid parr3c " <<dist << G4endl;

	dist=c8a.DistanceToOut(pparr3,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,50)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,2./sqrt(5.0),-1./sqrt(5.0))))
	    G4cout << "Error solid parr3d " <<dist << G4endl;


	dist=c8b.DistanceToOut(pparr2,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100*sqrt(5.)/2.)||
                     !*pgoodNorm||
                     OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr2a " <<dist << G4endl;
	dist=c8b.DistanceToOut(pparr2,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr2b " <<dist << G4endl;

	dist=c8b.DistanceToOut(pparr2,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,50)||*pgoodNorm)
	    G4cout << "Error hollow parr2c " <<dist << G4endl;
	dist=c8b.DistanceToOut(pparr2,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr2d " <<dist << G4endl;


	dist=c8b.DistanceToOut(pparr3,vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr3a " <<dist << G4endl;
	dist=c8b.DistanceToOut(pparr3,-vparr,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,100.*sqrt(5.)/2.)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vmz))
	    G4cout << "Error hollow parr3b " <<dist << G4endl;
	dist=c8b.DistanceToOut(pparr3,vz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,0)||
	    !*pgoodNorm||
	    OutRange(*pNorm,vz))
	    G4cout << "Error hollow parr3c " <<dist << G4endl;

	dist=c8b.DistanceToOut(pparr3,vmz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,50)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,2./sqrt(5.),-1.0/sqrt(5.))))
	    G4cout << "Error hollow parr3d " <<dist << G4endl;

	dist=c9.DistanceToOut(G4ThreeVector(1e3*kRadTolerance,0,50),
                              vx2mz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,111.8033988)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,0,-1.0)))
G4cout<<"Error:c9.Out((1e3*kRadTolerance,0,50),vx2mz,...) = " <<dist << G4endl;

	dist=c9.DistanceToOut(G4ThreeVector(5,0,50),
                              vx2mz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,111.8033988)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,0,-1.0)))
	    G4cout << "Error:c9.Out((5,0,50),vx2mz,...) = " <<dist << G4endl;

	dist=c9.DistanceToOut(G4ThreeVector(10,0,50),
                              vx2mz,calcNorm,pgoodNorm,pNorm);
	*pNorm=pNorm->unit();
	if (OutRange(dist,111.8033988)||
	    !*pgoodNorm||
	    OutRange(*pNorm,G4ThreeVector(0,0,-1.0)))
	    G4cout << "Error:c9.Out((10,0,50),vx2mz,...) = " <<dist << G4endl;

	dist=cms.DistanceToOut(
        G4ThreeVector(0.28628920024909,-0.43438111004815,-2949.0),
        G4ThreeVector(6.0886686196674e-05,-9.2382200635766e-05,0.99999999387917),
        calcNorm,pgoodNorm,pNorm);
	if (OutRange(dist,5898.0))
	G4cout << "Error:cms.DistToOut() =  " <<dist << G4endl;

	dist=cms.DistanceToOut(
        G4ThreeVector(0.28628920024909,-0.43438111004815,
                     -2949.0 + kCarTolerance*0.25),
        G4ThreeVector(6.0886686196674e-05,-9.2382200635766e-05,0.99999999387917),
        calcNorm,pgoodNorm,pNorm);
	if (OutRange(dist,5898.0))
	G4cout << "Error:cms.DistToOut(+) =  " <<dist << G4endl;

	dist=cms.DistanceToOut(
        G4ThreeVector(0.28628920024909,-0.43438111004815,
                      -2949.0 - kCarTolerance*0.25),
        G4ThreeVector(6.0886686196674e-05,-9.2382200635766e-05,0.99999999387917),
        calcNorm,pgoodNorm,pNorm);
	if (OutRange(dist,5898.0))
	G4cout << "Error:cms.DistToOut(-) =  " <<dist << G4endl;

	/////////////////////////////////////////////
	//

	G4cout << "Testing G4Cons::DistanceToIn(p) ...\n";


	dist=c1.DistanceToIn(pzero);
	if (OutRange(dist,50))
	  G4cout << "Error A " << dist << G4endl;

	dist=c1.DistanceToIn(pplx);
	if (OutRange(dist,20))
	  G4cout << "Error B " << dist << G4endl;

	dist=c1.DistanceToIn(pply);
	if (OutRange(dist,20))
	  G4cout << "Error C " << dist << G4endl;

	dist=c4.DistanceToIn(pply);
	if (OutRange(dist,120*sin(M_PI/3)))
	  G4cout << "Error D " << dist << G4endl;

	dist=c4.DistanceToIn(pmiy);
	if (OutRange(dist,120*sin(M_PI/3)))
	  G4cout << "Error D " << dist << G4endl;

	dist=c1.DistanceToIn(pplz);
	if (OutRange(dist,70))
	    G4cout << "Error E " << dist << G4endl;
// Check with both rmins=0
	dist=c5.DistanceToIn(pplx);
	if (OutRange(dist,20./sqrt(2.)))
	  G4cout << "Error F " << dist << G4endl;

	/////////////////////////////////////////////////////
	//

	G4cout << "Testing G4Cons::DistanceToIn(p,v,...) ...\n";

	dist=c1.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	  G4cout << "Error A " << dist << G4endl;

	dist=c2.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c2.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c3.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c3.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c4.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c4.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c5.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c5.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c6.DistanceToIn(pplz,vmz);
	if (OutRange(dist,70.0))
	G4cout << "Error:c6.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c7.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c7.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c8a.DistanceToIn(pplz,vmz);
	if (OutRange(dist,70.0))
	G4cout << "Error:c8a.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c8b.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c8b.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c8c.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c8c.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c9.DistanceToIn(pplz,vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c9.DistanceToIn(pplz,vmz) = " << dist << G4endl;

	dist=c9.DistanceToIn(G4ThreeVector(0,0,50),vmz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c9.DistanceToIn((0,0,50),vmz) = " << dist << G4endl;

	///////////////

	dist=c1.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error A " << dist << G4endl;

	dist=c2.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c2.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c3.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c3.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c4.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c4.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c5.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c5.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c6.DistanceToIn(pmiz,vz);
	if (OutRange(dist,70.0))
	G4cout << "Error:c6.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c7.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c7.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c8a.DistanceToIn(pmiz,vz);
	if (OutRange(dist,70.0))
	G4cout << "Error:c8a.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c8b.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c8b.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c8c.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c8c.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	dist=c9.DistanceToIn(pmiz,vz);
	if (OutRange(dist,kInfinity))
	G4cout << "Error:c9.DistanceToIn(pmiz,vz) = " << dist << G4endl;

	//////////////

	dist=c1.DistanceToIn(pplx,vmx);
	if (OutRange(dist,20))
	  G4cout << "Error B " << dist << G4endl;
	dist=c1.DistanceToIn(pplz,vx);
	if (OutRange(dist,kInfinity))
	  G4cout << "Error C " << dist << G4endl;
	dist=c4.DistanceToIn(pply,vmy);
	if (OutRange(dist,kInfinity))
	  G4cout << "Error D " << dist << G4endl;

	dist=c1.DistanceToIn(pydx,vmy);
	if (OutRange(dist,70))
	  G4cout << "Error E " << dist << G4endl;
	dist=c3.DistanceToIn(pydx,vmy);
	if (OutRange(dist,150-60*tan(M_PI/6)))
	  G4cout << "Error F " << dist << G4endl;

	dist=c1.DistanceToIn(pplx,vmx);
	if (OutRange(dist,20))
	  G4cout << "Error G " << dist << G4endl;
	dist=c1.DistanceToIn(pplx,vx);
	if (OutRange(dist,kInfinity))
	  G4cout << "Error G2 " << dist << G4endl;

	dist=c4.DistanceToIn(pbigx,vmx);
	if (OutRange(dist,350))
	    G4cout << "Error G3 " << dist << G4endl;

	dist=c4.DistanceToIn(pzero,vx);
	if (OutRange(dist,50))
	  G4cout << "Error H " << dist << G4endl;

	dist=c1.DistanceToIn(ponr2,vx);
	if (OutRange(dist,kInfinity))
	    G4cout << "Error I" << dist << G4endl;
	dist=c1.DistanceToIn(ponr2,vmx);
	if (OutRange(dist,0))
	    G4cout << "Error I2" << dist << G4endl;
	
	dist=c1.DistanceToIn(ponr1,vx);
	if (OutRange(dist,0))
	    G4cout << "Error J" << dist << G4endl;
	dist=c1.DistanceToIn(ponr1,vmx);
	if (OutRange(dist,2.0*sqrt(50*50/2.)))
	    G4cout << "Error J2" << dist << G4endl;

	dist=c1.DistanceToIn(ponr2,vmxmy);
	if (OutRange(dist,0))
	    G4cout << "Error K" << dist << G4endl;

// Parallel test case -> parallel to both radii
	dist=c8b.DistanceToIn(pparr1,vparr);
	if (OutRange(dist,100*sqrt(5.)/2.))
	    G4cout << "Error parr1 " << dist << G4endl;
	dist=c8b.DistanceToIn(pparr2,-vparr);
	if (OutRange(dist,kInfinity))
	    G4cout << "Error parr2 " << dist << G4endl;
	dist=c8b.DistanceToIn(pparr3,vparr);
	if (OutRange(dist,kInfinity))
	    G4cout << "Error parr3a " << dist << G4endl;
	dist=c8b.DistanceToIn(pparr3,-vparr);
	if (OutRange(dist,0))
	    G4cout << "Error parr3b " << dist << G4endl;

// Check we don't Hit `shadow cone' at `-ve radius' on rmax or rmin
	dist=c8a.DistanceToIn(proot1,vz);
	if (OutRange(dist,1000))
	    G4cout << "Error shadow rmax root problem " << dist << G4endl;

	dist=c8c.DistanceToIn(proot2,vz);
	if (OutRange(dist,1000))
	    G4cout << "Error shadow rmin root problem " << dist << G4endl;

	G4cout << "NOT Checking G4Cons::ScopeCar...\n";
	G4cout << "NOT Checking G4Cons::ScopePhi...\n";
	G4cout << "NOT Checking G4Cons::ScopeRad...\n";

	return 0;
}











