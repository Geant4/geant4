


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

#include "G4UnionSolid.hh"
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

    G4ThreeVector vxz(1/sqrt(2.0),0,1/sqrt(2.0));
    G4ThreeVector vxmz(1/sqrt(2.0),0,-1/sqrt(2.0));
    G4ThreeVector vmxz(-1/sqrt(2.0),0,1/sqrt(2.0));
    G4ThreeVector vmxmz(-1/sqrt(2.0),0,-1/sqrt(2.0));

    G4double dist;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4RotationMatrix identity, xRot ;
    
// NOTE: xRot = rotation such that x axis->-x axis & y axis->-y axis

    xRot.rotateZ(-M_PI) ;

    G4Transform3D transform(xRot,G4ThreeVector(0,40,0)) ;

    G4Box b1("Test Box #1",20,30,40);
    G4Box b2("Test Box #2",10,10,10);
    G4Box b3("Test Box #3",10,20,50);
    G4Box* b4 = new G4Box("b4",50,50,50) ;
    G4Box* b5 = new G4Box("b5",10,10,60) ;

    G4Tubs t1("Solid Tube #1",0,50,50,0,360);
    G4Tubs t2("Hole Tube #2",45,50,50,0,360);
    G4Tubs t3("disk",0,50,5,0,360);
    G4Tubs t4("t4",0,45,50,0,360);

    G4Cons c1("Hollow Full Tube",50,100,50,100,50,0,2*M_PI),
	   c2("Full Cone",0,50,0,100,50,0,2*M_PI) ;

    G4UnionSolid b1Ub2("b1Unionb2",&b1,&b2),
                        t1Ub2("t1Unionb2",&t1,&b2),
                        c2Ub2("c2Unionb2",&c2,&b2) ;

    G4UnionSolid b4Ub5("b4Ub5",b4,b5) ;

    G4UnionSolid t2Ut4("t2Ut4",&t2,&t4) ;

    // With placement

    G4UnionSolid b1Ub3("b1Unionb3",&b1,&b3,transform) ;

    G4UnionSolid t1Ub3("t1Unionb3",&t1,&b3,transform) ;

    G4UnionSolid t2Ut1Trans("t2Ut1Trans",&t2,&t1,
                                         &identity,G4ThreeVector(0,0,110)) ;

    G4UnionSolid t2Ut1Ut1("t2Ut1Ut1",&t2Ut1Trans,&t1,
                                         &identity,G4ThreeVector(0,0,-110)) ;

    G4UnionSolid t2Ut3("t2Ut3",&t2,&t3,
                                         &identity,G4ThreeVector(0,0,50)) ;

    G4UnionSolid t2Ut3Ut3("t2Ut3Ut3",&t2Ut3,&t3,
                                         &identity,G4ThreeVector(0,0,-50)) ;

    G4UnionSolid t3Ut3("t3Ut3",&t3,&t3,
                                         &identity,G4ThreeVector(0,0,10)) ;

    G4UnionSolid t3Ut3Ut3("t3Ut3Ut3",&t3Ut3,&t3,
                                         &identity,G4ThreeVector(0,0,-10)) ;




// Check Inside

    assert(b1Ub2.Inside(pzero)==kInside);
    assert(b1Ub2.Inside(pbigz)==kOutside);
    assert(b1Ub2.Inside(ponxside)==kSurface);
    assert(b1Ub2.Inside(ponyside)==kSurface);
    assert(b1Ub2.Inside(ponzside)==kSurface);

    assert(t1Ub2.Inside(pzero)==kInside);
    assert(t1Ub2.Inside(pbigz)==kOutside);
    assert(t1Ub2.Inside(ponb2x)==kInside);
    assert(t1Ub2.Inside(ponb2y)==kInside);
    assert(t1Ub2.Inside(ponb2z)==kInside);

    assert(c2Ub2.Inside(pzero)==kInside);
    assert(c2Ub2.Inside(pbigz)==kOutside);
    assert(c2Ub2.Inside(ponb2x)==kInside);
    assert(c2Ub2.Inside(ponb2y)==kInside);
    assert(c2Ub2.Inside(ponb2z)==kInside);

    // With placement

    assert(t1Ub3.Inside(pzero)==kInside);
    assert(t1Ub3.Inside(G4ThreeVector(0,60,0))==kSurface);
    assert(t1Ub3.Inside(G4ThreeVector(50,0,0))==kSurface);
    assert(t1Ub3.Inside(G4ThreeVector(0,65,0))==kOutside);
    assert(t1Ub3.Inside(G4ThreeVector(0,50,0))==kInside);


// Check Surface Normal

    G4ThreeVector normal;

    normal=b1Ub2.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    normal=b1Ub2.SurfaceNormal(ponmxside);
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

    normal=b1Ub2.SurfaceNormal(ponyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=b1Ub2.SurfaceNormal(ponmyside);
    assert(ApproxEqual(normal,G4ThreeVector(0,-1,0)));

    normal=b1Ub2.SurfaceNormal(ponzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,1)));

    normal=b1Ub2.SurfaceNormal(ponmzside);
    assert(ApproxEqual(normal,G4ThreeVector(0,0,-1)));

    // Wirh placemenet

    normal=t1Ub3.SurfaceNormal(G4ThreeVector(0,60,0));
    assert(ApproxEqual(normal,G4ThreeVector(0,1,0)));

    normal=t1Ub3.SurfaceNormal(G4ThreeVector(10,55,0));
    assert(ApproxEqual(normal,G4ThreeVector(1,0,0)));

    normal=t1Ub3.SurfaceNormal(G4ThreeVector(-10,55,0));
    assert(ApproxEqual(normal,G4ThreeVector(-1,0,0)));

    normal=t1Ub3.SurfaceNormal(G4ThreeVector(50.0/sqrt(2.0),50.0/sqrt(2.0),0));
    assert(ApproxEqual(normal,vxy));


// DistanceToOut(P)

    dist=b1Ub2.DistanceToOut(pzero);
    assert(ApproxEqual(dist,20));

    dist=b1Ub2.DistanceToOut(vx);
    assert(ApproxEqual(dist,19));

    // With placement

    dist=t1Ub3.DistanceToOut(pzero);
    assert(ApproxEqual(dist,50));
 
    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0));
    assert(ApproxEqual(dist,5));
 
    dist=t1Ub3.DistanceToOut(G4ThreeVector(45,0,0));
    assert(ApproxEqual(dist,5));
 
    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,45,0));
    assert(ApproxEqual(dist,10));
 
    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,-45,0));
    assert(ApproxEqual(dist,5));

 
// DistanceToOut(P,V)

    dist=b4Ub5.DistanceToOut(G4ThreeVector(50,0,0),vmx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,100));
    //    G4cout<<"(100) b4Ub5.DistanceToOut(G4ThreeVector(50,0,0),vmx) = "
    //          <<dist<<G4endl ;

    dist=b4Ub5.DistanceToOut(G4ThreeVector(-10,0,-60),vxz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,84.8528));
    //  G4cout<<"(84.85) b4Ub5.DistanceToOut(G4ThreeVector(-10,0,-60),vxz) = "
    //        <<dist<<G4endl ;


    dist=b1Ub2.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(*pNorm,vx)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmx)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vy)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmy)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,40)&&ApproxEqual(norm,vz)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,40)&&ApproxEqual(norm,vmz)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,sqrt(800.))); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2x,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(*pNorm,vx)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2x,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(*pNorm,vmx)); // &&*pgoodNorm);

    dist=b1.DistanceToOut(ponxside,vy,calcNorm,pgoodNorm,pNorm);
//  G4cout<<"b1.DistanceToOut(ponxside,vy) = "<<dist<<G4endl;
//  assert(ApproxEqual(dist,0)&&ApproxEqual(*pNorm,vy)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2mx,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)); //&&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2y,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vy)); // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2my,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,20)&&ApproxEqual(norm,vmy)) ; // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2z,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vz)) ; // &&*pgoodNorm);

    dist=b1Ub2.DistanceToOut(ponb2mz,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30)&&ApproxEqual(norm,vmz)) ; // &&*pgoodNorm);

    // With placement

    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0),vmy,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,105)&&ApproxEqual(norm,vmy)); // &&*pgoodNorm);

    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0),vy,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,5)&&ApproxEqual(norm,vy)); // &&*pgoodNorm);

    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vx)); // &&*pgoodNorm);

    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0),vmx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,10)&&ApproxEqual(norm,vmx)); // &&*pgoodNorm);

    dist=t1Ub3.DistanceToOut(G4ThreeVector(0,55,0),vz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,50)&&ApproxEqual(norm,vz)); // &&*pgoodNorm);

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,0),vmx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,5));
    // G4cout<<"(5) t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,0),vmx) = "
    //       <<dist<<G4endl ;

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,0),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0));
    //   G4cout<<"(0) t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,0),vx) = "
    //         <<dist<<G4endl ;

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,70),vmx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,100));
    //   G4cout<<"(100) t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,70),vmx) = "
    //         <<dist<<G4endl ;

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,70),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0));
    //   G4cout<<"(0) t2Ut1Ut1.DistanceToOut(G4ThreeVector(50,0,70),vx) = "
    //         <<dist<<G4endl ;

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(0,0,60),vmz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,0));
    //   G4cout<<"(0) t2Ut1Ut1.DistanceToOut(G4ThreeVector(0,0,60),vmz) = "
    //         <<dist<<G4endl ;

    dist=t2Ut1Ut1.DistanceToOut(G4ThreeVector(0,0,60),vz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,100));
    //  G4cout<<"(100) t2Ut1Ut1.DistanceToOut(G4ThreeVector(0,0,60),vz) = "
    //        <<dist<<G4endl ;

    dist=t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-48),vxmz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,9.89949));
    //  G4cout<<"(9.89) t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-48),vxmz) = "
    //        <<dist<<G4endl ;

    dist=t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-50),vxz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,7.07107));
    //  G4cout<<"(7.07) t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-50),vxz) = "
    //        <<dist<<G4endl ;

    dist=t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-50),vxmz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,7.07107));
    // G4cout<<"(7.07) t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-50),vxmz) = "
    //       <<dist<<G4endl ;

    dist=t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-52),vxmz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,4.24264));
    // G4cout<<"(4.24) t2Ut3Ut3.DistanceToOut(G4ThreeVector(-50,0,-52),vxmz) = "
    //       <<dist<<G4endl ;

    dist=t3Ut3Ut3.DistanceToOut(G4ThreeVector(0,0,-15),vz,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,30));
    //  G4cout<<"(30) t3Ut3Ut3.DistanceToOut(G4ThreeVector(0,0,-15),vz) = "
    //        <<dist<<G4endl ;

    dist=t2Ut4.DistanceToOut(G4ThreeVector(-50,0,0),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,100));
    // G4cout<<"(100) t2Ut4.DistanceToOut(G4ThreeVector(-50,0,0),vx) = "
    //       <<dist<<G4endl ;

    dist=t2Ut4.DistanceToOut(G4ThreeVector(-15,-45,0),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,36.7945));
    // G4cout<<"(36.8) t2Ut4.DistanceToOut(G4ThreeVector(-15,-45,0),vx) = "
    //    <<dist<<G4endl ;

    dist=t2Ut4.DistanceToOut(G4ThreeVector(0,-45,0),vx,
                             calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(dist,21.7945));
    // G4cout<<"(21.8) t2Ut4.DistanceToOut(G4ThreeVector(0,-45,0),vx) = "
    //      <<dist<<G4endl ;

// DistanceToIn(P)

    dist=b1Ub2.DistanceToIn(pbigx);
    assert(ApproxEqual(dist,80));
    //   G4cout<<"b1Ub2.DistanceToIn(pbigx) = "<<dist<<G4endl ;
    //   G4cout<<"b1.DistanceToIn(pbigx) = "<<b1.DistanceToIn(pbigx)<<G4endl ;
    //   G4cout<<"b2.DistanceToIn(pbigx) = "<<b2.DistanceToIn(pbigx)<<G4endl ;

    dist=b1Ub2.DistanceToIn(ponxside);
    assert(ApproxEqual(dist,0));

    dist=b1Ub2.DistanceToIn(ponxside);
    assert(ApproxEqual(dist,0));

    dist=b1Ub2.DistanceToIn(pbigmy);
    assert(ApproxEqual(dist,70));

    dist=b1Ub2.DistanceToIn(pbigz);
    assert(ApproxEqual(dist,60));

    dist=b1Ub2.DistanceToIn(pbigmz);
    assert(ApproxEqual(dist,60));

    // With placement

    dist=t1Ub3.DistanceToIn(pbigmz);
    assert(ApproxEqual(dist,50));

    dist=t1Ub3.DistanceToIn(pbigy);
    assert(ApproxEqual(dist,40));

    dist=t1Ub3.DistanceToIn(G4ThreeVector(0,65,0));
    assert(ApproxEqual(dist,5));

    dist=t1Ub3.DistanceToIn(G4ThreeVector(0,60,0));
    assert(ApproxEqual(dist,0));


// DistanceToIn(P,V)

    dist=b1Ub2.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(dist,80));

    dist=b1Ub2.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(dist,80));

    dist=b1Ub2.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(dist,70));

    dist=b1Ub2.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(dist,70));

    dist=b1Ub2.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(dist,60));

    dist=b1Ub2.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(dist,60));

    dist=b1Ub2.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(dist,kInfinity));

    dist=b1Ub2.DistanceToIn(pbigmx,vxy);
    assert(ApproxEqual(dist,kInfinity));


    // With placement

    dist=t1Ub3.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(dist,50));

    dist=t1Ub3.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(dist,40));

    dist=t1Ub3.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(dist,50));

    dist=t1Ub3.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(dist,50));

    dist=t2Ut1Trans.DistanceToIn(pzero,vx);
    //  G4cout<<"t2Ut1Trans.DistanceToIn(pzero,vx) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,45));

    dist=t2Ut1Trans.DistanceToIn(pzero,vz);
    //   G4cout<<"t2Ut1Trans.DistanceToIn(pzero,vz) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,60));

    dist=t2Ut1Trans.DistanceToIn(pzero,vmz);
    //   G4cout<<"t2Ut1Trans.DistanceToIn(pzero,vmz) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,kInfinity));

    dist=t2Ut1Trans.DistanceToIn(G4ThreeVector(0,0,55),vx);
    //   G4cout<<"t2Ut1Trans.DistanceToIn(G4ThreeVector(0,0,55),vx) = "
    //         <<dist<<G4endl ;
    assert(ApproxEqual(dist,kInfinity));

    dist=t2Ut1Ut1.DistanceToIn(pzero,vmz);
    //  G4cout<<"t2Ut1Ut1.DistanceToIn(pzero,vmz) = "<<dist<<G4endl ;
    assert(ApproxEqual(dist,60));


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

    assert(b1Ub3.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-20)&&ApproxEqual(max,20));

    assert(b1Ub3.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-30)&&ApproxEqual(max,60));

    assert(b1Ub3.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(ApproxEqual(min,-50)&&ApproxEqual(max,50));

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

    G4cout<<"CalculateExtent is OK"<<G4endl ;


  return 0 ;
}
