// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4Box.cc,v 1.1 1999-01-08 16:31:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// testG4Box
//             Ensure asserts are compiled in

#include <assert.h>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

G4bool testG4Box()
{
    G4ThreeVector pzero(0,0,0);
    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40);
    G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
    G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40);

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxy(-1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/sqrt(2.0),-1/sqrt(2.0),0);
    G4ThreeVector vxmy(1/sqrt(2.0),-1/sqrt(2.0),0);
    G4ThreeVector vxmz(1/sqrt(2.0),0,-1/sqrt(2.0));

    G4double Dist;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4Box b1("Test Box #1",20,30,40);
    G4Box b2("Test Box #2",10,10,10);

    G4ThreeVector pJohnXZ(9,0,12);
    Dist = b2.DistanceToIn(pJohnXZ,vxmz) ;
    G4cout<<"b2.DistanceToIn(pJohnXZ,vxmz) = "<<Dist<<endl ;

    G4ThreeVector pJohnXY(12,9,0);
    Dist = b2.DistanceToIn(pJohnXY,vmxy) ;
    G4cout<<"b2.DistanceToIn(pJohnXY,vmxy) = "<<Dist<<endl ;

    Dist = b2.DistanceToIn(pJohnXY,vmx) ;
    G4cout<<"b2.DistanceToIn(pJohnXY,vmx) = "<<Dist<<endl ;

    G4ThreeVector pMyXY(32,-11,0);
    Dist = b2.DistanceToIn(pMyXY,vmxy) ;
    G4cout<<"b2.DistanceToIn(pMyXY,vmxy) = "<<Dist<<endl ;



// Check name
    assert(b1.GetName()=="Test Box #1");

// Check Inside
    assert(b1.Inside(pzero)==kInside);
    assert(b1.Inside(pbigz)==kOutside);
    assert(b1.Inside(ponxside)==kSurface);
    assert(b1.Inside(ponyside)==kSurface);
    assert(b1.Inside(ponzside)==kSurface);

// Check Surface Normal
    G4ThreeVector normal;

    normal=b1.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));
    normal=b1.SurfaceNormal(ponmxside);
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));
    normal=b1.SurfaceNormal(ponyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));
    normal=b1.SurfaceNormal(ponmyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));
    normal=b1.SurfaceNormal(ponzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=b1.SurfaceNormal(ponmzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));
    normal=b1.SurfaceNormal(ponzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));
    normal=b1.SurfaceNormal(ponmzsidey);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));


// DistanceToOut(P)
    Dist=b1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,20));
    Dist=b1.DistanceToOut(vx);
    assert(ApproxEqual(Dist,19));
    Dist=b1.DistanceToOut(vy);
    assert(ApproxEqual(Dist,20));
    Dist=b1.DistanceToOut(vz);
    assert(ApproxEqual(Dist,20));

// DistanceToOut(P,V)
    Dist=b1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,20)&&ApproxEqual(norm,vmx)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,vy)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,30)&&ApproxEqual(norm,vmy)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(norm,vmz)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,sqrt(800.))&&*pgoodNorm);

    Dist=b1.DistanceToOut(ponxside,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,40)&&ApproxEqual(*pNorm,vmx)&&*pgoodNorm);
    Dist=b1.DistanceToOut(pbigx,vy,calcNorm,pgoodNorm,pNorm);
G4cout<<"b1.DistanceToOut(ponxside,vy) = "<<Dist<<endl;
//   assert(ApproxEqual(Dist,0)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponmxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmx)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponyside,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vy)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponmyside,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmy)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponzside,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vz)&&*pgoodNorm);
    Dist=b1.DistanceToOut(ponmzside,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,0)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

//DistanceToIn(P)
    Dist=b1.DistanceToIn(pbigx);
    assert(ApproxEqual(Dist,80));
    Dist=b1.DistanceToIn(pbigmx);
    assert(ApproxEqual(Dist,80));
    Dist=b1.DistanceToIn(pbigy);
    assert(ApproxEqual(Dist,70));
    Dist=b1.DistanceToIn(pbigmy);
    assert(ApproxEqual(Dist,70));
    Dist=b1.DistanceToIn(pbigz);
    assert(ApproxEqual(Dist,60));
    Dist=b1.DistanceToIn(pbigmz);
    assert(ApproxEqual(Dist,60));

// DistanceToIn(P,V)
    Dist=b1.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(Dist,80));
    Dist=b1.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(Dist,80));
    Dist=b1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,70));
    Dist=b1.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(Dist,70));
    Dist=b1.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(Dist,60));
    Dist=b1.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,60));
    Dist=b1.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));
    Dist=b1.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(Dist,kInfinity));

// CalculateExtent
    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    assert(b1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));
    assert(b1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(b1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);
    assert(b1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-120)&&ApproxEqual(max,-80));
    assert(b1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-140)&&ApproxEqual(max,-80));
    assert(b1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-160)&&ApproxEqual(max,-80));

    G4RotationMatrix r90Z;
    r90Z.rotateZ(M_PI/2);
    G4AffineTransform tRotZ(r90Z,pzero);
    assert(b1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(b1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));
    assert(b1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

// Check that clipped away
    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-50);
    assert(!b1.CalculateExtent(kXAxis,xClip,origin,min,max));

// Assert clipped to volume
    G4VoxelLimits allClip;
    allClip.AddLimit(kXAxis,-5,+5);
    allClip.AddLimit(kYAxis,-5,+5);
    allClip.AddLimit(kZAxis,-5,+5);
    G4RotationMatrix genRot;
    genRot.rotateX(M_PI/6);
    genRot.rotateY(M_PI/6);
    genRot.rotateZ(M_PI/6);
    G4AffineTransform tGen(genRot,vx);
    assert(b1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(b1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(b1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    G4VoxelLimits buggyClip2;
    buggyClip2.AddLimit(kXAxis,5,15);
    assert(b1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(b1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(b1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    buggyClip2.AddLimit(kYAxis,5,15);
    assert(b1.CalculateExtent(kXAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(b1.CalculateExtent(kYAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,5)&&ApproxEqual(max,15));
    assert(b1.CalculateExtent(kZAxis,buggyClip2,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4VoxelLimits buggyClip1;
    buggyClip1.AddLimit(kXAxis,-5,+5);
    assert(b1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(b1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));
    assert(b1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));
    buggyClip1.AddLimit(kYAxis,-5,+5);
    assert(b1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(b1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));
    assert(b1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    /* **********************************************************
    */ /////////////////////////////////////////////////////

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Box());
    return 0;
}

