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
// $Id: testG4Tubs.cc,v 1.12 2002-04-16 08:24:19 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// testG4Tubs
//
//  Test file for class G4Tubs [NOT thorough]
//
//             Ensure asserts are compiled in
//
// History:
//
// ~95-96 P. Kent R&D
// 21.5.99 V. Grichine tests of phi-intersections (t3 and t4)

#include <assert.h>
#include <math.h>

#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4Tubs.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

///////////////////////////////////////////////////////////////////
//
// Dave's auxiliary function

const G4String OutputInside(const EInside a)
{
	switch(a) 
        {
		case kInside:  return "Inside"; 
		case kOutside: return "Outside";
		case kSurface: return "Surface";
	}
	return "????";
}


G4bool testG4Tubs()
{
    G4cout.precision(16) ;

    G4ThreeVector pzero(0,0,0);

    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector ponxside(50,0,0);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxy(-1/sqrt(2.0),1/sqrt(2.0),0);
    G4ThreeVector vmxmy(-1/sqrt(2.0),-1/sqrt(2.0),0);
    G4ThreeVector vxmy(1/sqrt(2.0),-1/sqrt(2.0),0);

    G4double Dist;
    G4ThreeVector *pNorm,norm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

    G4Tubs t1("Solid Tube #1",0,50*mm,50*mm,0,2*pi);

    G4Tubs t2("Hole Tube #2",45*mm,50*mm,50*mm,0,2*pi);

    G4Tubs t2a("Hole Tube #2",5*mm,50*mm,50*mm,0,2*pi);

    G4Tubs t2b("Hole Tube #2",15*mm,50*mm,50*mm,0,2*pi);

    G4Tubs t2c("Hole Tube #2",25*mm,50*mm,50*mm,0,2*pi);

    G4Tubs t2d("Hole Tube #2",35*mm,50*mm,50*mm,0,2*pi);

    G4Tubs t3("Solid Sector #3",0,50*mm,50*mm,halfpi,halfpi);

    G4Tubs t4("Hole Sector #4",45*mm,50*mm,50*mm,halfpi,halfpi);

  G4Tubs t5("Hole Sector #5",50*mm,100*mm,50*mm,0.0,270.0*deg);

  G4Tubs tube6("tube6",750,760,350,0.31415926535897931,5.6548667764616276);

  G4Tubs tube7("tube7",2200,3200,2500,-0.68977164349384879,3.831364227270472);

  G4Tubs tube8("tube8",2550,2580,2000,0,2*pi);

  G4Tubs tube9("tube9",1150,1180,2000,0,2*pi);

  G4Tubs tube10("tube10",400*mm,405*mm,400*mm,0*degree,360*degree) ;


// Check name
    assert(t1.GetName()=="Solid Tube #1");

// Check Inside

	//
	// Make a tub
	//
	G4Tubs *arc = new G4Tubs( "outer", 1*m, 1.1*m, 0.01*m, -15*deg, 30*deg );
	
	//
	// First issue: 
	//   A point on the start phi surface just beyond the
	//   start angle but still well within tolerance 
	//   is found to be "outside" by G4Tubs::Inside
	//
	//   pt1 = exactly on phi surface (within precision)
	//   pt2 = t1 but slightly higher, and still on tolerant surface
	//   pt3 = t1 but slightly lower, and still on tolerant surface
	//
	G4ThreeVector pt1( 1.05*m*cos(-15*deg),
	                   1.05*m*sin(-15*deg),
			      0*m );
 			  
        G4ThreeVector pt2 = pt1 + G4ThreeVector(0,0.001*kCarTolerance,0) ;
        G4ThreeVector pt3 = pt1 - G4ThreeVector(0,0.001*kCarTolerance,0) ;
	
	EInside a1 = arc->Inside(pt1);
	EInside a2 = arc->Inside(pt2);
	EInside a3 = arc->Inside(pt3);
	
	// G4cout << "Point pt1 is " << OutputInside(a1) << G4endl;

        assert(a1==kSurface);
	// G4cout << "Point pt2 is " << OutputInside(a2) << G4endl;
        assert(a2==kSurface);
	// G4cout << "Point pt3 is " << OutputInside(a3) << G4endl;
	assert(a3==kSurface);


    assert(t1.Inside(pzero)==kInside);
    assert(t1.Inside(pbigx)==kOutside);

    EInside in = t5.Inside(G4ThreeVector(60,-0.001*kCarTolerance,0)) ;
    assert(in == kSurface);
    //    G4cout<<"t5.Inside(G4ThreeVector(60,-0.001*kCarTolerance,0)) = "
    //     <<OutputInside(in)<<G4endl;
    in = tube10.Inside(G4ThreeVector(-114.8213313833317*mm,
					   382.7843220719649*mm,
                                           -32.20788536438663*mm)) ;
    //  assert(in == kSurface);
    G4cout<<"tube10.Inside(G4ThreeVector(-114.821...)) = "
         <<OutputInside(in)<<G4endl;

// Check Surface Normal

    G4ThreeVector normal;

    normal=t1.SurfaceNormal(ponxside);
    assert(ApproxEqual(normal,vx));

// DistanceToOut(P)
    Dist=t1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,50));

// DistanceToOut(P,V)
    Dist=t1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vx)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmx)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vy)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmy)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
    Dist=t1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);

    Dist=t2.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut(pzero,vxy) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(ponxside,vmx,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut(ponxside,vmx) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(ponxside,vmxmy,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut(ponxside,vmxmy) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(ponxside,vz,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut(ponxside,vz) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(pbigx,vx,calcNorm,pgoodNorm,pNorm);
    //   G4cout<<"Dist=t2.DistanceToOut(pbigx,vx) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(pbigx,vxy,calcNorm,pgoodNorm,pNorm);
    //   G4cout<<"Dist=t2.DistanceToOut(pbigx,vxy) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(pbigx,vz,calcNorm,pgoodNorm,pNorm);
    //   G4cout<<"Dist=t2.DistanceToOut(pbigx,vz) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(G4ThreeVector(45.5,0,0),vx,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut((45.5,0,0),vx) = "<<Dist<<G4endl;

    Dist=t2.DistanceToOut(G4ThreeVector(49.5,0,0),vx,calcNorm,pgoodNorm,pNorm);
    //  G4cout<<"Dist=t2.DistanceToOut((49.5,0,0),vx) = "<<Dist<<G4endl;


    Dist=t3.DistanceToOut(G4ThreeVector(0,10,0),vx,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((0,10,0),vx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0));

    Dist=t3.DistanceToOut(G4ThreeVector(0.5,10,0),vx,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((0.5,10,0),vx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,48.489795));

    Dist=t3.DistanceToOut(G4ThreeVector(-0.5,9,0),vx,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((-0.5,9,0),vx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,0.5));

    Dist=t3.DistanceToOut(G4ThreeVector(-5,9.5,0),vx,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((-5,9.5,0),vx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,5));

    Dist=t3.DistanceToOut(G4ThreeVector(-5,9.5,0),vmy,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((-5,9.5,0),vmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,9.5));

    Dist=t3.DistanceToOut(G4ThreeVector(-5,9,0),vxmy,calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=t3.DistanceToOut((-5,9,0),vxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,7.0710678));

    // bug #76
    Dist=tube6.DistanceToOut(
    G4ThreeVector(-388.20504321896431,-641.71398957741451,332.85995254027955),
    G4ThreeVector(-0.47312863350457468,-0.782046391443315, 0.40565100491504164),
    calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=tube6.DistanceToOut(p,v) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,10.940583));

    // bug #91
    Dist=tube7.DistanceToOut(
    G4ThreeVector(-2460,1030,-2500),
    G4ThreeVector(-0.086580540180167642,0.070084247882560638,0.9937766390194761),
    calcNorm,pgoodNorm,pNorm);
    // G4cout<<"Dist=tube7.DistanceToOut(p,v) = "<<Dist<<G4endl;
    // assert(ApproxEqual(Dist,4950.348576972614));

    Dist=tube8.DistanceToOut(
 G4ThreeVector(6.71645645882942,2579.415860329989,-1.519530725281157),
 G4ThreeVector(-0.6305220496340839,-0.07780451841562354,0.7722618738739774),
    calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=tube8.DistanceToOut(p,v) = "<<Dist<<G4endl;
    // assert(ApproxEqual(Dist,4950.348576972614));

    Dist=tube9.DistanceToOut(
 G4ThreeVector(2.267347771505638,1170.164934028592,4.820317321984064),
 G4ThreeVector(-0.1443054266272111,-0.01508874701037938,0.9894181489944458),
    calcNorm,pgoodNorm,pNorm);
    G4cout<<"Dist=tube9.DistanceToOut(p,v) = "<<Dist<<G4endl;
    // assert(ApproxEqual(Dist,4950.348576972614));


    G4cout<<G4endl ;


//DistanceToIn(P)

    Dist=t1.DistanceToIn(pbigx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigz);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmz);
    assert(ApproxEqual(Dist,50));

// DistanceToIn(P,V)

    Dist=t1.DistanceToIn(pbigx,vmx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmx,vx);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmy,vy);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigz,vmz);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigmz,vz);
    assert(ApproxEqual(Dist,50));
    Dist=t1.DistanceToIn(pbigx,vxy);
    assert(ApproxEqual(Dist,kInfinity));

    Dist=t2.DistanceToIn(G4ThreeVector(45.5,0,0),vx);
    //  G4cout<<"Dist=t2.DistanceToIn((45.5,0,0),vx) = "<<Dist<<G4endl;
   
    Dist=t2.DistanceToIn(G4ThreeVector(45.5,0,0),vmx);
    //  G4cout<<"Dist=t2.DistanceToIn((45.5,0,0),vmx) = "<<Dist<<G4endl;
   
    Dist=t2.DistanceToIn(G4ThreeVector(49.5,0,0),vmx);
    //  G4cout<<"Dist=t2.DistanceToIn((49.5,0,0),vmx) = "<<Dist<<G4endl;
   
    Dist=t2.DistanceToIn(G4ThreeVector(49.5,0,0),vx);
    //   G4cout<<"Dist=t2.DistanceToIn((49.5,0,0),vx) = "<<Dist<<G4endl;
   
    Dist=t3.DistanceToIn(G4ThreeVector(49.5,0,0),vmx);
    //  G4cout<<"Dist=t2.DistanceToIn((49.5,0,0),vmx) = "<<Dist<<G4endl;
   
    Dist=t3.DistanceToIn(G4ThreeVector(49.5,5,0),vmx);
    //  G4cout<<"Dist=t2.DistanceToIn((49.5,5,0),vmx) = "<<Dist<<G4endl;
   
    Dist=t3.DistanceToIn(G4ThreeVector(49.5,-0.5,0),vmx);
    //  G4cout<<"Dist=t2.DistanceToIn((49.5,-0.5,0),vmx) = "<<Dist<<G4endl;
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-20.0,0),vxy);
    // G4cout<<"Dist=t5.DistanceToIn((30.0,-20.0,0),vxy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,28.284271));
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-70.0,0),vxy);
    // G4cout<<"Dist=t5.DistanceToIn((30.0,-70.0,0),vxy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,kInfinity));
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-20.0,0),vmxmy);
    //  G4cout<<"Dist=t5.DistanceToIn((30.0,-20.0,0),vmxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,42.426407));
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-70.0,0),vmxmy);
    // G4cout<<"Dist=t5.DistanceToIn((30.0,-70.0,0),vmxmy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,kInfinity));
   
    Dist=t5.DistanceToIn(G4ThreeVector(50.0,-20.0,0),vy);
    // G4cout<<"Dist=t5.DistanceToIn((50.0,-20.0,0),vy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,20));

    Dist=t5.DistanceToIn(G4ThreeVector(100.0,-20.0,0),vy);
    // G4cout<<"Dist=t5.DistanceToIn((100.0,-20.0,0),vy) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,kInfinity));
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-50.0,0),vmx);
    //  G4cout<<"Dist=t5.DistanceToIn((30.0,-50.0,0),vmx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,30));
   
    Dist=t5.DistanceToIn(G4ThreeVector(30.0,-100.0,0),vmx);
    //  G4cout<<"Dist=t5.DistanceToIn((30.0,-100.0,0),vmx) = "<<Dist<<G4endl;
    assert(ApproxEqual(Dist,kInfinity));
   


// CalculateExtent

    G4VoxelLimits unlimit;		// Unlimited

    G4VoxelLimits limitX, limitY, limitZ, limitXYZ,limitXsYZ, limitXYsZ;

    limitX.AddLimit(kXAxis,-20,-10);
    limitY.AddLimit(kYAxis,30,40);
    limitZ.AddLimit(kZAxis,-40,-30);

    limitXYZ.AddLimit(kXAxis,-20,-10);
    limitXYZ.AddLimit(kYAxis,30,40);
    limitXYZ.AddLimit(kZAxis,-40,-30);

    limitXsYZ.AddLimit(kXAxis,-60,-10);
    limitXsYZ.AddLimit(kYAxis,30,40);
    limitXsYZ.AddLimit(kZAxis,-40,-30);

    limitXYsZ.AddLimit(kXAxis,-20,-10);
    limitXYsZ.AddLimit(kYAxis,30,60);
    limitXYsZ.AddLimit(kZAxis,-40,-30);

    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    G4bool clipped;

    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);

    G4RotationMatrix r90Z;
    r90Z.rotateZ(halfpi);
    G4AffineTransform tRotZ(r90Z,pzero);
    G4AffineTransform tRotZpos(r90Z,pmxmymz);


    assert(t1.CalculateExtent(kXAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kYAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kZAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=50);

    assert(t3.CalculateExtent(kXAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=0);
    assert(t3.CalculateExtent(kYAxis,unlimit,origin,min,max));
    assert(min<=0&&max>=50);
    assert(t3.CalculateExtent(kZAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=50);

    /////////////////////////////////////////////////////

    assert(t1.CalculateExtent(kXAxis,limitXYZ,origin,min,max));
    G4cout<<"t1.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    assert(t1.CalculateExtent(kYAxis,limitXYZ,origin,min,max));
    G4cout<<"t1.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    assert(t1.CalculateExtent(kZAxis,limitXYZ,origin,min,max));
    G4cout<<"t1.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t2.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t2.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t2.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2a.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t2a.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2a.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t2a.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2a.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t2a.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2b.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t2b.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2b.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t2b.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2b.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t2b.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2c.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t2c.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2c.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t2c.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2c.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t2c.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2d.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t2d.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2d.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t2d.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2d.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t2d.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);


    clipped=t3.CalculateExtent(kXAxis,limitXYZ,origin,min,max);
    G4cout<<"t3.CE(kXAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);


    clipped=t3.CalculateExtent(kYAxis,limitXYZ,origin,min,max);
    G4cout<<"t3.CE(kYAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t3.CalculateExtent(kZAxis,limitXYZ,origin,min,max);
    G4cout<<"t3.CE(kZAxis,limitXYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t1.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t1.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t1.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t1.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t1.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t1.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2a.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2a.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2a.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2a.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2a.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2a.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2b.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2b.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2b.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2b.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2b.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2b.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2c.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2c.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2c.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2c.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2c.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2c.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2d.CalculateExtent(kXAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2d.CE(kXAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2d.CalculateExtent(kYAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2.CE(kYAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2.CalculateExtent(kZAxis,limitXsYZ,origin,min,max);
    G4cout<<"t2d.CE(kZAxis,limitXsYZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);


     ///////////////////////////////////////////////////////////
    

    assert(t3.CalculateExtent(kXAxis,limitX,origin,min,max));
    G4cout<<"t3.CE(kXAxis,limitX,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=-20&&max>=-30);



    assert(t3.CalculateExtent(kYAxis,limitY,origin,min,max));
    G4cout<<"t3.CE(kYAxis,limitY,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=30&&max>=40);

    assert(t3.CalculateExtent(kZAxis,limitZ,origin,min,max));
    G4cout<<"t3.CE(kZAxis,limitZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    assert(min<=-40&&max>=-30);
    
    assert(t2.CalculateExtent(kXAxis,limitX,origin,min,max));
    G4cout<<"t2.CE(kXAxis,limitX,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=-20&&max>=-30);

    assert(t2.CalculateExtent(kYAxis,limitY,origin,min,max));
    G4cout<<"t2.CE(kYAxis,limitY,origin,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=30&&max>=40);

    assert(t2.CalculateExtent(kZAxis,limitZ,origin,min,max));
    G4cout<<"t2.CE(kZAxis,limitZ,origin,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    assert(min<=-40&&max>=-30);
    
    assert(t4.CalculateExtent(kXAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=0);
    assert(t4.CalculateExtent(kYAxis,unlimit,origin,min,max));
    assert(min<=0&&max>=50);
    assert(t4.CalculateExtent(kZAxis,unlimit,origin,min,max));
    assert(min<=-50&&max>=50);

    //////////////////////////////////////////////////////////
    

    assert(t1.CalculateExtent(kXAxis,unlimit,tPosOnly,min,max));
    assert(min<=-150&&max>=-50);
    assert(t1.CalculateExtent(kYAxis,unlimit,tPosOnly,min,max));
    assert(min<=-160&&max>=-60);
    assert(t1.CalculateExtent(kZAxis,unlimit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);

    assert(t3.CalculateExtent(kXAxis,unlimit,tPosOnly,min,max));
    assert(min<=-150&&max>=-100);
    assert(t3.CalculateExtent(kYAxis,unlimit,tPosOnly,min,max));
    assert(min<=-110&&max>=-60);
    assert(t3.CalculateExtent(kZAxis,unlimit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);

    assert(t4.CalculateExtent(kXAxis,unlimit,tPosOnly,min,max));
    assert(min<=-150&&max>=-100);
    assert(t4.CalculateExtent(kYAxis,unlimit,tPosOnly,min,max));
    assert(min<=-110&&max>=-60);
    assert(t4.CalculateExtent(kZAxis,unlimit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);


    /////////////////////////////////////////////////////////////

    assert(t1.CalculateExtent(kXAxis,unlimit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kYAxis,unlimit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(t1.CalculateExtent(kZAxis,unlimit,tRotZ,min,max));
    assert(min<=-50&&max>=50);

    assert(t3.CalculateExtent(kXAxis,unlimit,tRotZ,min,max));
    G4cout<<"t3.CE(kXAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=0&&max>=50);

    assert(t3.CalculateExtent(kYAxis,unlimit,tRotZ,min,max));
    G4cout<<"t3.CE(kYAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=0&&max>=50);

    assert(t3.CalculateExtent(kZAxis,unlimit,tRotZ,min,max));
    G4cout<<"t3.CE(kZAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    assert(min<=-50&&max>=50);

    clipped=t1.CalculateExtent(kXAxis,limitX,tRotZ,min,max);
    G4cout<<"t1.CE(kXAxis,limitX,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t1.CalculateExtent(kYAxis,limitY,tRotZ,min,max);
    G4cout<<"t1.CE(kYAxis,limitY,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t1.CalculateExtent(kZAxis,limitZ,tRotZ,min,max);
    G4cout<<"t1.CE(kZAxis,limitZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t3.CalculateExtent(kXAxis,limitX,tRotZ,min,max);
    G4cout<<"t3.CE(kXAxis,limitX,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t3.CalculateExtent(kYAxis,limitY,tRotZ,min,max);
    G4cout<<"t3.CE(kYAxis,limitY,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t3.CalculateExtent(kZAxis,limitZ,tRotZ,min,max);
    G4cout<<"t3.CE(kZAxis,limitZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2.CalculateExtent(kXAxis,limitX,tRotZ,min,max);
    G4cout<<"t2.CE(kXAxis,limitX,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2.CalculateExtent(kYAxis,limitY,tRotZ,min,max);
    G4cout<<"t2.CE(kYAxis,limitY,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2.CalculateExtent(kZAxis,limitZ,tRotZ,min,max);
    G4cout<<"t2.CE(kZAxis,limitZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    /* *******************************
    ******************************** */

    clipped=t1.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t1.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);


    clipped=t1.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t1.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t1.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t1.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2a.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2a.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2a.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2a.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2a.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2a.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2b.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2b.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2b.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2b.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2b.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2b.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2c.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2c.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2c.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2c.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2c.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2c.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t2d.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2d.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t2d.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2d.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t2d.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t2d.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);

    clipped=t3.CalculateExtent(kXAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t3.CE(kXAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=-20&&max>=-30);

    clipped=t3.CalculateExtent(kYAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t3.CE(kYAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=30&&max>=40);

    clipped=t3.CalculateExtent(kZAxis,limitXYZ,tRotZ,min,max);
    G4cout<<"t3.CE(kZAxis,limitXYZ,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-40&&max>=-30);


    assert(t4.CalculateExtent(kXAxis,unlimit,tRotZ,min,max));
    G4cout<<"t4.CE(kXAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=0&&max>=50);

    assert(t4.CalculateExtent(kYAxis,unlimit,tRotZ,min,max));
    G4cout<<"t4.CE(kYAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl;
    // assert(min<=0&&max>=50);

    assert(t4.CalculateExtent(kZAxis,unlimit,tRotZ,min,max));
    G4cout<<"t4.CE(kZAxis,unlimit,tRotZ,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    // assert(min<=-50&&max>=50);

    assert(t3.CalculateExtent(kXAxis,unlimit,tRotZpos,min,max));
    G4cout<<"t3.CE(kXAxis,unlimit,tRotZpos,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=-100&&max>=-50);

    assert(t3.CalculateExtent(kYAxis,unlimit,tRotZpos,min,max));
    G4cout<<"t3.CE(kYAxis,unlimit,tRotZpos,min = "
          <<min<<"; max = "<<max<<G4endl;
    assert(min<=-110&&max>=-60);

    assert(t3.CalculateExtent(kZAxis,unlimit,tRotZpos,min,max));
    G4cout<<"t3.CE(kZAxis,unlimit,tRotZpos,min = "
          <<min<<"; max = "<<max<<G4endl<<G4endl;
    assert(min<=-170&&max>=-70);

// Check that clipped away
    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-60);
    assert(!t1.CalculateExtent(kXAxis,xClip,origin,min,max));

// Assert clipped to volume
    G4VoxelLimits allClip;
    allClip.AddLimit(kXAxis,-5,+5);
    allClip.AddLimit(kYAxis,-5,+5);
    allClip.AddLimit(kZAxis,-5,+5);

    G4RotationMatrix genRot;
    genRot.rotateX(pi/6);
    genRot.rotateY(pi/6);
    genRot.rotateZ(pi/6);
    G4AffineTransform tGen(genRot,vx);

    assert(t1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(t1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(t1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);


// Test t1 z clipping ok

    for (G4double zTest=-100;zTest<100;zTest+=9)
    {
      G4VoxelLimits zTestClip;
      zTestClip.AddLimit(kZAxis,-kInfinity,zTest);
      if (zTest<-50)
      {
	assert(!t1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
      }
      else
      {
	assert(t1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
	G4double testMin=-50;
	G4double testMax=(zTest<50) ? zTest : 50;
	assert ( ApproxEqual(min,testMin) && ApproxEqual(max,testMax) );
      }
    }
    G4cout<<"Test t1 z clipping ok"<<G4endl;

// Test t1 y clipping ok

    for (G4double xTest=-100;xTest<100;xTest+=9)
    {
      G4VoxelLimits xTestClip;
      xTestClip.AddLimit(kXAxis,-kInfinity,xTest);
      if (xTest<-50)
      {
        assert(!t1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
      }
      else
      {
	assert(t1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
// Calc max y coordinate
// G4double testMax=(xTest<0) ? sqrt(50*50-xTest*xTest) : 50;
// assert (ApproxEqual(min,-testMax)&&ApproxEqual(max,testMax));
      }
    }
    G4cout<<"Test t1 y clipping ok"<<G4endl;

// Test t1 x clipping ok

    for (G4double yTest=-100;yTest<100;yTest+=9)
    {
      G4VoxelLimits yTestClip;
      yTestClip.AddLimit(kYAxis,-kInfinity,yTest);
      if (yTest<-50)
      {
	assert(!t1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
      }
      else
      {
	assert(t1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
// Calc max y coordinate
	      //  G4double testMax=(yTest<0) ? sqrt(50*50-yTest*yTest) : 50;
	      //  assert (ApproxEqual(min,-testMax)&&ApproxEqual(max,testMax));
      }
    }
    G4cout<<"Test t1 x clipping ok"<<G4endl;


    /* ********************************
    ************************************ */

    return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4Tubs());
    return 0;
}







