// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Para2.cc,v 1.1 1999-01-08 16:31:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Test of G4Para
// Includes all/most of the tests Done for a box

#include <assert.h>
#include <math.h>
#include "G4ios.hh"

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Para.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"



//#include "G4ios.hh"
//#include "globals.hh"
//#include "G4Para.hh"

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
    G4double Dist,low,high;
    G4ThreeVector pzero(0,0,0),px(10,0,0),py(0,10,0),pz(0,0,10);
    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);
    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40);
    G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
    G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40);
    G4RotationMatrix runit;
    G4RotationMatrix r90X,r90Y,r90Z,r45X,r30Y;
    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1,1,0);
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    r90X.rotateX(M_PI_2);
    r90Y.rotateY(M_PI_2);
    r90Z.rotateZ(M_PI_2);
    r45X.rotateX(M_PI_4);
    r30Y.rotateY(M_PI/6);

    vxy=vxy.unit();

    G4Para p1("Box",20,30,40,0,0,0),
	p2("2",50,50,50,M_PI/6,0,0),
	p3("3",50,50,50,0,M_PI/6,0),
	p4("4",50,50,50,0,0,M_PI/6),
	p5("5",50,50,50,0,M_PI/6,M_PI/6),	
	p6("6",50,50,50,M_PI/6,M_PI/6,M_PI/6);	

    G4cout << "Name:"<< p1.GetName()
	 << " ID=" <<endl;
    
    G4cout << "Checking G4Para::Inside...\n";
    if (p1.Inside(pzero)!=kInside)
	G4cout << "Error A" << endl;
    if (p1.Inside(pbigz)!=kOutside)
	G4cout << "Error B" << endl;
    if (p1.Inside(ponxside)!=kSurface)
	G4cout << "Error C" << endl;
    if (p1.Inside(ponyside)!=kSurface)
	G4cout << "Error D" << endl;
    if (p1.Inside(ponzside)!=kSurface)
	G4cout << "Error E" << endl;



    G4cout << "Checking G4Para::SurfaceNormal...\n";
    norm=p1.SurfaceNormal(ponxside);
    if (OutRange(norm,G4ThreeVector(1,0,0)))
	G4cout << "Error A " << norm << endl;
    norm=p1.SurfaceNormal(ponmxside);
    if (OutRange(norm,G4ThreeVector(-1,0,0)))
	G4cout << "Error B " << norm << endl;
    norm=p1.SurfaceNormal(ponyside);
    if (OutRange(norm,G4ThreeVector(0,1,0)))
	G4cout << "Error C " << norm << endl;
    norm=p1.SurfaceNormal(ponmyside);
    if (OutRange(norm,G4ThreeVector(0,-1,0)))
	G4cout << "Error D " << norm << endl;
    norm=p1.SurfaceNormal(ponzside);
    if (OutRange(norm,G4ThreeVector(0,0,1)))
	G4cout << "Error E " << norm << endl;
    norm=p1.SurfaceNormal(ponmzside);
    if (OutRange(norm,G4ThreeVector(0,0,-1)))
	G4cout << "Error F " << norm << endl;
    norm=p1.SurfaceNormal(ponzsidey);
    if (OutRange(norm,G4ThreeVector(0,0,1)))
	G4cout << "Error G " << norm << endl;
    norm=p1.SurfaceNormal(ponmzsidey);
    if (OutRange(norm,G4ThreeVector(0,0,-1)))
	G4cout << "Error H " << norm << endl;


    G4cout << "Checking G4Para::DistanceToOut(P)...\n";
    Dist=p1.DistanceToOut(pzero);
    if (OutRange(Dist,20))
	G4cout << "Error A1 " << Dist << endl;
    Dist=p2.DistanceToOut(pzero);
    if (OutRange(Dist,50*cos(M_PI/6)))
	G4cout << "Error A2 " << Dist << endl;
     Dist=p3.DistanceToOut(pzero);
    if (OutRange(Dist,50*cos(M_PI/6)))
	G4cout << "Error A3 " << Dist << endl;
    Dist=p5.DistanceToOut(pzero);
    if (OutRange(Dist,2*50/sqrt(5.)))
	G4cout << "Error A4 " << Dist << endl;

     Dist=p1.DistanceToOut(px);
    if (OutRange(Dist,10))
	G4cout << "Error B " << Dist << endl;
    Dist=p1.DistanceToOut(py);
    if (OutRange(Dist,20))
	G4cout << "Error C " << Dist << endl;
    Dist=p1.DistanceToOut(pz);
    if (OutRange(Dist,20))
	G4cout << "Error D " << Dist << endl;




    G4cout << "Checking G4Para::DistanceToOut(P,V)...\n";

    Dist=p1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,20)||OutRange(*pNorm,vx)||!*pgoodNorm)
	G4cout << "Error A " << Dist << endl;
    Dist=p1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,20)||OutRange(norm,vmx)||!*pgoodNorm)
 	G4cout << "Error B " << Dist << endl;
    Dist=p1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,30)||OutRange(norm,vy)||!*pgoodNorm)
 	G4cout << "Error C " << Dist << endl;
    Dist=p1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,30)||OutRange(norm,vmy)||!*pgoodNorm)
 	G4cout << "Error D " << Dist << endl;
     Dist=p1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,40)||OutRange(norm,vz)||!*pgoodNorm)
 	G4cout << "Error E " << Dist << endl;
     Dist=p1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,40)||OutRange(norm,vmz)||!*pgoodNorm)
 	G4cout << "Error F " << Dist << endl;
    Dist=p1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,sqrt(800.))||!*pgoodNorm)
 	G4cout << "Error F " << Dist << endl;

    Dist=p1.DistanceToOut(ponxside,vx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(*pNorm,vx)||!*pgoodNorm)
	G4cout << "Error A2 " << Dist << endl;
    Dist=p1.DistanceToOut(ponmxside,vmx,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmx)||!*pgoodNorm)
 	G4cout << "Error B2 " << Dist << endl;
    Dist=p1.DistanceToOut(ponyside,vy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vy)||!*pgoodNorm)
 	G4cout << "Error C2 " << Dist << endl;
    Dist=p1.DistanceToOut(ponmyside,vmy,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmy)||!*pgoodNorm)
 	G4cout << "Error D2 " << Dist << endl;
     Dist=p1.DistanceToOut(ponzside,vz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vz)||!*pgoodNorm)
 	G4cout << "Error E2 " << Dist << endl;
    Dist=p1.DistanceToOut(ponmzside,vmz,calcNorm,pgoodNorm,pNorm);
    if (OutRange(Dist,0)||OutRange(norm,vmz)||!*pgoodNorm)
 	G4cout << "Error F2 " << Dist << endl;
  
    G4cout << "Checking G4Para::DistanceToIn(P)...\n";
    Dist=p1.DistanceToIn(pbigx);
    if (OutRange(Dist,80))
 	G4cout << "Error A " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmx);
    if (OutRange(Dist,80))
 	G4cout << "Error B " << Dist << endl;
    Dist=p1.DistanceToIn(pbigy);
    if (OutRange(Dist,70))
 	G4cout << "Error C " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmy);
    if (OutRange(Dist,70))
 	G4cout << "Error D " << Dist << endl;
    Dist=p1.DistanceToIn(pbigz);
    if (OutRange(Dist,60))
 	G4cout << "Error E " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmz);
    if (OutRange(Dist,60))
 	G4cout << "Error F " << Dist << endl;

    Dist=p3.DistanceToIn(pbigx);
    if (OutRange(Dist,50*cos(M_PI/6)))
	G4cout << "Error G1 " << Dist <<endl;
    Dist=p3.DistanceToIn(pbigy);
    if (OutRange(Dist,50))
	G4cout << "Error G2 " << Dist <<endl;

    G4cout << "Checking G4Para::DistanceToIn(P,V)...\n";

    Dist=p1.DistanceToIn(pbigx,vmx);
    if (OutRange(Dist,80))
 	G4cout << "Error A " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmx,vx);
    if (OutRange(Dist,80))
 	G4cout << "Error B " << Dist << endl;
    Dist=p1.DistanceToIn(pbigy,vmy);
    if (OutRange(Dist,70))
 	G4cout << "Error C " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmy,vy);
    if (OutRange(Dist,70))
 	G4cout << "Error D " << Dist << endl;
    Dist=p1.DistanceToIn(pbigz,vmz);
    if (OutRange(Dist,60))
 	G4cout << "Error E " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmz,vz);
    if (OutRange(Dist,60))
 	G4cout << "Error F " << Dist << endl;
    Dist=p1.DistanceToIn(pbigx,vxy);
    if (OutRange(Dist,kInfinity))
 	G4cout << "Error G " << Dist << endl;
    Dist=p1.DistanceToIn(pbigmx,vxy);
    if (OutRange(Dist,kInfinity))
 	G4cout << "Error H " << Dist << endl;

    return 0;    
}




