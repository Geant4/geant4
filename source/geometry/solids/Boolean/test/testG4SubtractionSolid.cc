

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

#include "G4SubtractionSolid.hh"
// #include "G4DisplacedSolid.hh"


int main()
{
    G4ThreeVector pzero(0,0,0);

    G4ThreeVector ponxside(20,0,0),ponyside(0,30,0),ponzside(0,0,40),

                    ponb2x(10,0,0),  ponb2y(0,10,0),  ponb2z(0,0,10),

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
    
// NOTE: xRot = rotation such that x axis->-x axis & y axis->-y axis

    xRot.rotateZ(-M_PI) ;

    G4Transform3D transform(xRot,G4ThreeVector(0,30,0)) ;

    G4Box b1("Test Box #1",20,30,40);
    G4Box b2("Test Box #2",10,10,10);
    G4Box b3("Test Box #3",10,20,50);
    G4Box b4("Test Box #4",20,20,40);

    G4Tubs t1("Solid Tube #1",0,50,50,0,360);
    
    // t2\t3 for DistanceToIn

    G4Tubs t2("Hole Tube #2",50,60,50,0,2*M_PI); 
 
    G4Tubs t3("Hole Tube #3",45,55,50,M_PI/4.,M_PI*3./2.);

    G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*M_PI),

	   c2("Full Cone",0,50,0,100,50,0,2*M_PI) ;

    G4SubtractionSolid b1Sb2("b1Sb2",&b1,&b2),

                       t1Sb2("t1Sb2",&t1,&b2),

                       t2St3("t2St3",&t2,&t3),

                       c2Sb2("c2Sb2",&c2,&b2) ;

    // Tube t1 \ box t3, which was rotated by -pi and translated +30y

    G4SubtractionSolid   t1Sb3("t1Subtractionb3",&t1,&b3,transform) ;
    G4SubtractionSolid   b1Sb4("t1Subtractionb3",&b1,&b4,transform) ;

// Check Inside

    assert(b1Sb2.Inside(pzero)==kOutside);
    assert(b1Sb2.Inside(pbigz)==kOutside);
    assert(b1Sb2.Inside(ponb2x)==kSurface);
    assert(b1Sb2.Inside(ponb2y)==kSurface);
    assert(b1Sb2.Inside(ponb2z)==kSurface);

    assert(t1Sb3.Inside(pzero)==kInside);
    assert(t1Sb3.Inside(pbigz)==kOutside);
    assert(t1Sb3.Inside(ponb2x)==kInside);
    assert(t1Sb3.Inside(ponb2y)==kSurface);
    assert(t1Sb3.Inside(ponb2z)==kInside);

    assert(c2Sb2.Inside(pzero)==kOutside);
    assert(c2Sb2.Inside(pbigz)==kOutside);
    assert(c2Sb2.Inside(ponb2x)==kSurface);
    assert(c2Sb2.Inside(ponb2y)==kSurface);
    assert(c2Sb2.Inside(ponb2z)==kSurface);

// Check Surface Normal

    G4ThreeVector normal;

    normal=b1Sb2.SurfaceNormal(ponb2x);
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

    normal=b1Sb2.SurfaceNormal(ponb2mx);
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    normal=b1Sb2.SurfaceNormal(ponb2y);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));

    normal=b1Sb2.SurfaceNormal(ponb2my);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=b1Sb2.SurfaceNormal(ponb2z);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    normal=b1Sb2.SurfaceNormal(ponb2mz);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

    normal=b1Sb2.SurfaceNormal(ponb2zy);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    normal=b1Sb2.SurfaceNormal(ponb2mzy);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));


    normal=t1Sb3.SurfaceNormal(ponb2y);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));


// DistanceToOut(P)

    dist=b1Sb2.DistanceToOut(ponb2x);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponb2y);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponb2z);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToOut(ponxside);
    assert(ApproxEqual(dist,0));

    dist=t1Sb3.DistanceToOut(ponb2y);
    assert(ApproxEqual(dist,0));

    dist=t1Sb3.DistanceToOut(pzero);
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToOut(ponb2x);
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToOut(ponb2z);
    assert(ApproxEqual(dist,10));


// DistanceToOut(P,V)

    dist=b1Sb2.DistanceToOut(ponb2x,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponxside,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponyside,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmy)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2z,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponzside,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmz)); // &&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,sqrt(200))&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2x,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vmx)); // &&*pgoodNorm);

    dist=b1.DistanceToOut(ponxside,vy,calcNorm,pgoodNorm,pNorm);
//  G4cout<<"b1.DistanceToOut(ponxside,vy) = "<<dist<<G4endl;
//  assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2mx,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2my,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2z,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vz)&&*pgoodNorm);

    dist=b1Sb2.DistanceToOut(ponb2mz,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmz)&&*pgoodNorm);

    // With placement

    dist=t1Sb3.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vy)); // &&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vx)&&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vmy)&&*pgoodNorm);

    dist=t1Sb3.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vz)&&*pgoodNorm);


//DistanceToIn(P)

    dist=b1Sb2.DistanceToIn(pbigx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(ponb2x);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToIn(ponxside);
    assert(ApproxEqual(dist,0));

    dist=b1Sb2.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigmz);
    assert(ApproxEqual(dist,60));

    // With placement

    dist=t1Sb3.DistanceToIn(pbigy);
    assert(ApproxEqual(dist,50));

    dist=t1Sb3.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,50));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,20,0));
    assert(ApproxEqual(dist,10));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,15,0));
    assert(ApproxEqual(dist,5));

// DistanceToIn(P,V)

    dist=b1Sb2.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(dist,80));

    dist=b1Sb2.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(dist,70));

    dist=b1Sb2.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(dist,60));

    dist=b1Sb2.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Sb2.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    // Cases with a number of iterations inside 'while'

    dist=t2St3.DistanceToIn(G4ThreeVector(2.5,-52.5,0),vy) ;
    assert(ApproxEqual(dist,107.443));
// G4cout<<"t2St3.DistanceToIn(G4ThreeVector(2.5,-52.5,0),vy) = "<<dist<<G4endl ;

    dist=t2St3.DistanceToIn(G4ThreeVector(2.5,-62.5,0),vy) ;
    assert(ApproxEqual(dist,2.55211));
//    G4cout<<"t2St3.DistanceToIn(G4ThreeVector(2.5,-62.5,0),vy) = "<<dist<<G4endl ;

    // With placement

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,100,0),vmy);
    assert(ApproxEqual(dist,50));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vmy);
    assert(ApproxEqual(dist,26));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vy);
    assert(ApproxEqual(dist,kInfinity));

    dist=t1Sb3.DistanceToIn(G4ThreeVector(0,36,0),vx);
    assert(ApproxEqual(dist,10));

    G4cout<<"Tracking functions are OK"<<G4endl ;

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

    assert(b1Sb4.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));

    assert(b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max));
    //  G4cout<<"min of b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<min<<G4endl ;
    // G4cout<<"max of b1Sb4.CalculateExtent(kYAxis,limit,origin,min,max) = "
    //      <<max<<G4endl ;
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,30));

    assert(b1Sb4.CalculateExtent(kZAxis,limit,origin,min,max));
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

    G4cout<<"CalculateExtent is OK "<<G4endl ;

  G4int out =0 ;
  return out ;
}








