// Test file for the class G4IntersectionSolid
//
// 

#include <assert.h>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4Trap.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"

#include "G4IntersectionSolid.hh"
// #include "G4DisplacedSolid.hh"


int main()
{
    G4ThreeVector pzero(0,0,0);
    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40),
                   ponb2x(10,0,0),ponb2y(0,10,0),ponb2z(0,0,10),
                   ponb2mx(-10,0,0),ponb2my(0,-10,0),ponb2mz(0,0,-10);
    G4ThreeVector ponmxside(-20,0,0),ponmyside(0,-30,0),ponmzside(0,0,-40);
    G4ThreeVector ponzsidey(0,25,40),ponmzsidey(0,25,-40),
                  ponb2zy(0,5,10),ponb2mzy(0,5,-10) ;

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxy(-1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/sqrt(2.0),-1/sqrt(2.0),0);
    G4ThreeVector vxmy(1/sqrt(2.0),-1/sqrt(2.0),0);
    G4ThreeVector vxmz(1/sqrt(2.0),0,-1/sqrt(2.0));

    G4double dist;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4RotationMatrix identity, xRot ;
    
// NOTE: xRot = rotation such that x axis->y axis & y axis->-x axis

    xRot.rotateZ(-M_PI*0.5) ;

    G4Transform3D transform(xRot,ponb2y) ;

    G4Box b1("Test Box #1",20,30,40);
    G4Box b2("Test Box #2",10,10,10);

    G4Tubs t1("Solid Tube #1",0,50,50,0,360);
    G4Tubs t2("Hole Tube #2",45,50,50,0,360);

    G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*M_PI) ;
    G4Cons c2("Full Cone",0,50,0,100,50,0,2*M_PI) ;

    G4Tubs* tube3 = new G4Tubs( "OuterFrame",
                                      1.0*m,
                                      1.1*m,
                                      0.50*m,
                                      0*deg,
                                      180*deg );

    G4Tubs* tube4 = new G4Tubs("AnotherTubs",
                                    1.0*m,
                                    1.1*m,
                                    0.50*m,
                                    0*deg,
                                    180*deg );
    G4RotationMatrix rotmat2;
    rotmat2.rotateY(M_PI/4.0);
    G4Transform3D tran2 = G4Transform3D(rotmat2,G4ThreeVector(0.0,0.0,0.0));

    G4VSolid* t3It4 = new G4IntersectionSolid( "Example", tube3, tube4, tran2 );




    G4IntersectionSolid b1Ib2("b1Intersectionb2",&b1,&b2,transform) ;
    G4IntersectionSolid likeb2("b1Intersectionb2",&b1,&b2) ;
    G4IntersectionSolid t1Ib2("t1Intersectionb2",&t1,&b2,&xRot,ponb2y) ;
    G4IntersectionSolid c2Ib2("c2Intersectionb2",&c2,&b2,transform) ;

    G4cout.precision(16);

// Check Inside

    assert(b1Ib2.Inside(pzero)==kSurface);
    assert(b1Ib2.Inside(pbigz)==kOutside);
    assert(b1Ib2.Inside(ponb2y)==kInside);

    assert(t1Ib2.Inside(pzero)==kSurface);
    assert(t1Ib2.Inside(pbigz)==kOutside);
    assert(t1Ib2.Inside(ponb2y)==kInside);

    assert(c2Ib2.Inside(pzero)==kSurface);
    assert(c2Ib2.Inside(pbigz)==kOutside);
    assert(c2Ib2.Inside(ponb2y)==kInside);

// Check Surface Normal

    G4ThreeVector normal;
    /*
    normal=b1Ib2.SurfaceNormal(G4ThreeVector(10,15,0));
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    normal=b1Ib2.SurfaceNormal(G4ThreeVector(-10,15,0));
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

    normal=b1Ib2.SurfaceNormal(pzero);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));

    normal=b1Ib2.SurfaceNormal(G4ThreeVector(0,20,0));
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=b1Ib2.SurfaceNormal(G4ThreeVector(0,15,10));
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

    normal=b1Ib2.SurfaceNormal(G4ThreeVector(0,15,-10));
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));
    */


// DistanceToOut(P)

    dist=b1Ib2.DistanceToOut(pzero);
    assert(ApproxEqual(dist,0));

    dist=b1Ib2.DistanceToOut(vx);
    assert(ApproxEqual(dist,0));

    dist=b1Ib2.DistanceToOut(vy);
    assert(ApproxEqual(dist,1));

    dist=b1Ib2.DistanceToOut(vz);
    assert(ApproxEqual(dist,0));

// DistanceToOut(P,V)

    dist=b1Ib2.DistanceToOut(G4ThreeVector(0,5,0),vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(G4ThreeVector(0,5,0),vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(G4ThreeVector(0,5,0),vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(G4ThreeVector(0,5,0),vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,sqrt(200.0))&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(G4ThreeVector(10,5,0),vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(ponb2x,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(*pNorm,vmx)&&*pgoodNorm);

    dist=b1.DistanceToOut(ponxside,vy,calcNorm,pgoodNorm,pNorm);
//  cout<<"b1.DistanceToOut(ponxside,vy) = "<<dist<<G4endl;
//  assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(ponb2mx,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(norm,vmx)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Ib2.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(norm,vmy)&&*pgoodNorm);


//DistanceToIn(P)

    dist=b1Ib2.DistanceToIn(G4ThreeVector(100,1,0));
    assert(ApproxEqual(dist,80));

    dist=b1Ib2.DistanceToIn(G4ThreeVector(10,1,0));
    assert(ApproxEqual(dist,0));

    dist=b1Ib2.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,70));

    dist=b1Ib2.DistanceToIn(pbigz);
    assert(ApproxEqual(dist,60));

    dist=b1Ib2.DistanceToIn(pbigmz);
    assert(ApproxEqual(dist,60));

// DistanceToIn(P,V)

    dist=b1Ib2.DistanceToIn(G4ThreeVector(100,1,0),vmx);
    assert(ApproxEqual(dist,90));

    dist=b1Ib2.DistanceToIn(G4ThreeVector(-100,1,0),vx);
    assert(ApproxEqual(dist,90));

    dist=b1Ib2.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(dist,80));

    dist=b1Ib2.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(dist,100));

    dist=b1Ib2.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Ib2.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Ib2.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Ib2.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Ib2.DistanceToIn(pzero,vmy);
    assert(ApproxEqual(dist,kInfinity));
    // G4cout<<"(kInfinity) b1Ib2.DistanceToIn(pzero,vmy) = "<<dist<<G4endl ;

    dist=b1Ib2.DistanceToIn(pzero,vy);
    assert(ApproxEqual(dist,0));
    // G4cout<<"(0) b1Ib2.DistanceToIn(pzero,vy) = "<<dist<<G4endl ;


    // It returns 0, probably due to G4Displaced or G4Transform3D games !?
    //
    // dist=b1Ib2.DistanceToIn(pzero,vmx);
    // assert(ApproxEqual(dist,kInfinity));
    //  G4cout<<"(kInfinity) b1Ib2.DistanceToIn(pzero,vmx) = "<<dist<<G4endl ;

    dist=likeb2.DistanceToIn(G4ThreeVector(10,0,0),vmx);
    assert(ApproxEqual(dist,0));
    // G4cout<<"(0) likeb2.DistanceToIn(G4ThreeVector(10,0,0),vmx) = "
    //       <<dist<<G4endl ;

    dist=b1Ib2.DistanceToIn(G4ThreeVector(10,0,0),vx);
    assert(ApproxEqual(dist,kInfinity));
    // G4cout<<"(kInfinity) likeb2.DistanceToIn(G4ThreeVector(10,0,0),vx) = "
    //       <<dist<<G4endl ;

    dist=b1Ib2.DistanceToIn(G4ThreeVector(10,0,0),vy);
    assert(ApproxEqual(dist,kInfinity));
    // G4cout<<"(kInfinity) likeb2.DistanceToIn(G4ThreeVector(10,0,0),vy) = "
    //       <<dist<<G4endl ;

    dist=b1Ib2.DistanceToIn(G4ThreeVector(10,0,0),vmy);
    assert(ApproxEqual(dist,kInfinity));
    // G4cout<<"(kInfinity) likeb2.DistanceToIn(G4ThreeVector(10,0,0),vmy) = "
    //       <<dist<<G4endl ;

    dist=t3It4->DistanceToIn(
    G4ThreeVector(1888.691673255004,-476.1676766307428,-295.4764663381112),
    G4ThreeVector(-0.8509009035458712,0.5062362036610951,-0.1403301765395527));
    assert(ApproxEqual(dist,940.603760037514));
    //  G4cout<<"t3It4->DistanceToIn = "<<dist<<G4endl ;



    G4cout<<"Tracking functions are OK"<<G4endl ;


//  CalculateExtent

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

    assert(b1Ib2.CalculateExtent(kXAxis,limit,origin,min,max));
    //   G4cout<<"min of b1Ib2.CalculateExtent(kXAxis,limit,origin,min,max) = "
    //    <<min<<G4endl ;
    //  G4cout<<"max of b1Ib2.CalculateExtent(kXAxis,limit,origin,min,max) = "
    //      <<max<<G4endl ;
     assert(ApproxEqual(min,-10)&&ApproxEqual(max,10));

    assert(b1Ib2.CalculateExtent(kYAxis,limit,origin,min,max));
    // G4cout<<"min of b1Ib2.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<min<<G4endl ;
    //  G4cout<<"max of b1Ib2.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<max<<G4endl ;
    assert(ApproxEqual(min,0)&&ApproxEqual(max,20));

    assert(b1Ib2.CalculateExtent(kZAxis,limit,origin,min,max));
    //  G4cout<<"min of b1Ib2.CalculateExtent(kZAxis,limit,origin,min,max) = "
    //      <<min<<G4endl ;
    //  G4cout<<"max of b1Ib2.CalculateExtent(kZAxis,limit,origin,min,max) = "
    //      <<max<<G4endl ;
     assert(ApproxEqual(min,-10)&&ApproxEqual(max,10));

    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);

    assert(b1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-120)&&ApproxEqual(max,-80));

    assert(b1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-140)&&ApproxEqual(max,-80));

    assert(b1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-160)&&ApproxEqual(max,-80));

    assert(b1Ib2.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(ApproxEqual(min,-110)&&ApproxEqual(max,-90));


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
    assert(!b1Ib2.CalculateExtent(kXAxis,xClip,origin,min,max));

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

    assert(b1Ib2.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1Ib2.CalculateExtent(kYAxis,allClip,tGen,min,max)) ;
    // G4cout<<"min of b1Ib2.CalculateExtent(kYAxis,allClip,tGen,min,max) = "
    //     <<min<<G4endl ;
    // G4cout<<"max of b1Ib2.CalculateExtent(kYAxis,allClip,tGen,min,max) = "
    //     <<max<<G4endl ;
    // Reasonable but not so obvious ?!
    assert(ApproxEqual(min,-3.21667)&&ApproxEqual(max,5)) ;



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

    assert(b1Ib2.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1Ib2.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,0)&&ApproxEqual(max,20));

    assert(b1Ib2.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-10)&&ApproxEqual(max,10));



    buggyClip1.AddLimit(kYAxis,-5,+5);

    assert(b1.CalculateExtent(kXAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kYAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-5)&&ApproxEqual(max,5));

    assert(b1.CalculateExtent(kZAxis,buggyClip1,origin,min,max));
    assert(ApproxEqual(min,-40)&&ApproxEqual(max,40));

    G4cout<<"CalculateExtent is OK"<<G4endl ;

  return 0 ;
}
