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

// $Id: testG4Orb.cc,v 1.3 2004-09-10 09:38:29 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// G4Orb Test File
//
// o Basic asserts on each function +
//   awkward cases for tracking / geom algorithms
//
// o Add tests on dicovering bugs in G4Orb.cc...
//
// History:
// 28.03.95 P.Kent Initial version
// 20.10.96 V.Grichine Final modifications to commit

#include "G4ios.hh"
#include <assert.h>
#include <math.h>
#include "globals.hh"
#include "geomdefs.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4Orb.hh"

//const G4double kApproxEqualTolerance = kCarTolerance;

// Return true if the double check is approximately equal to target
//
// Process:
//
// Return true is difference < kApproxEqualTolerance

//G4bool ApproxEqual(const G4double check,const G4double target)
//{
//    return (fabs(check-target)<kApproxEqualTolerance) ? true : false ;
//}

// Return true if the 3vector check is approximately equal to target
//G4bool ApproxEqual(const G4ThreeVector& check, const G4ThreeVector& target)
//{
//    return (ApproxEqual(check.x(),target.x())&&
//	   ApproxEqual(check.y(),target.y())&&
//	    ApproxEqual(check.z(),target.z()))? true : false;
//}

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



int main(void)
{
    G4double Dist, vol, volCheck;
    G4ThreeVector pzero(0,0,0),px(30,0,0),py(0,30,0),pz(0,0,30);
    G4ThreeVector Pmx(-30,0,0),pmy(0,-30,0),pmz(0,0,-30);
    G4ThreeVector pbigx(100,0,0),pbigy(0,100,0),pbigz(0,0,100);
    G4ThreeVector pbigmx(-100,0,0),pbigmy(0,-100,0),pbigmz(0,0,-100);

    G4ThreeVector ponrmin1(45,0,0),ponrmax1(50,0,0),ponzmax(0,0,50),
	    ponrmin2(45/sqrt(2.),45/sqrt(2.),0),
            ponrmin3(0,0,-45),ponrminJ(0,0,-300),ponrmaxJ(0,0,-500),
	    ponrmax2(50/sqrt(2.),50/sqrt(2.),0);
    G4ThreeVector ponphi1(48/sqrt(2.),-48/sqrt(2.),0),
	          ponphi2(48/sqrt(2.),48/sqrt(2.),0),
	          pInPhi(48*0.866,-24,0),
	          pOverPhi(-48/sqrt(2.),48/sqrt(2.),0);
    G4ThreeVector pontheta1(0,48*sin(pi/4),48*cos(pi/4)),
	    pontheta2(0,48*sin(pi/4),-48*cos(pi/4));

    G4ThreeVector ptestphi1(-100,-45/sqrt(2.),0),
	    ptestphi2(-100,45/sqrt(2.),0);

    G4ThreeVector ptesttheta1(0,48/sqrt(2.),100),
	    ptesttheta2(0,48/sqrt(2.),-100);

    G4ThreeVector vx(1,0,0),vy(0,1,0),vz(0,0,1);
    G4ThreeVector vmx(-1,0,0),vmy(0,-1,0),vmz(0,0,-1);
    G4ThreeVector vxy(1/sqrt(2.),1/sqrt(2.),0),vmxmy(-1/sqrt(2.),-1/sqrt(2.),0);
    G4ThreeVector vxmy(1/sqrt(2.),-1/sqrt(2.),0),vmxy(-1/sqrt(2.),1/sqrt(2.),0);
    G4ThreeVector v345exit1(-0.8,0.6,0),v345exit2(0.8,0.6,0),
	          v345exit3(0.6,0.8,0);
    G4ThreeVector norm,*pNorm;
    G4bool *pgoodNorm,goodNorm,calcNorm=true;

    pNorm=&norm;
    pgoodNorm=&goodNorm;

  G4Orb s1("Solid G4Orb",50);

  G4Orb s10("s10",0.018*mm);

  G4Orb* solidD1= new G4Orb("D1", 2.7*cm);

  G4ThreeVector s10p(0.01160957408065766*mm,
                     0.01308205826682229*mm,
                     0.004293345210644617*mm);

  G4ThreeVector positionKip(-5.1427112977805294,
                            25.37671986392073,
                            -7.6534050889634502);
  G4ThreeVector directionKip(-0.90020960513809678,
                             -0.052042885743481759,
                              0.43233575477931813);


  G4cout.precision(16);


#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif

    //////////////// Check name /////////////////////////

    assert(s1.GetName()=="Solid G4Orb");


    // check cubic volume cases

    vol = s1.GetCubicVolume();
    volCheck = 4*pi*50*50*50/3;
    assert(ApproxEqual(vol,volCheck));
    
// Check G4Orb::Inside


    EInside   inside = s10.Inside(s10p);
    G4cout<<"s10.Inside(s10p) = "
          <<OutputInside(inside)<<G4endl ;
    G4cout<<"p radius = "
          <<s10p.mag()<<G4endl ;
    

    assert(s1.Inside(pzero)==kInside);
    // assert(s1.Inside(pz)==kInside);

// Checking G4Orb::SurfaceNormal
    norm=s1.SurfaceNormal(ponrmax1);
    assert(ApproxEqual(norm,vx));

// Checking G4Orb::DistanceToOut(P)
    Dist=s1.DistanceToOut(pzero);
    assert(ApproxEqual(Dist,50));
    Dist=s1.DistanceToOut(ponrmax1);
    assert(ApproxEqual(Dist,0));

// Checking G4Orb::DistanceToOut(p,v)


        Dist=s1.DistanceToOut(pz,vz,calcNorm,pgoodNorm,pNorm);
	G4cout<<"Dist=s1.DistanceToOut(pz,vz) = "<<Dist<<G4endl;

     Dist=s1.DistanceToOut(ponrmax1,vx,calcNorm,pgoodNorm,pNorm);
     *pNorm=pNorm->unit();
     assert(ApproxEqual(Dist,0)&&*pgoodNorm&&ApproxEqual(*pNorm,vx));


    Dist=s1.DistanceToOut(pzero,vx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vx)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmx,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmx)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vy)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmy)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);

    Dist=s1.DistanceToOut(pzero,vmz,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vmz)&&*pgoodNorm);
    Dist=s1.DistanceToOut(pzero,vxy,calcNorm,pgoodNorm,pNorm);
    assert(ApproxEqual(Dist,50)&&ApproxEqual(pNorm->unit(),vxy)&&*pgoodNorm);
    
    

     
// Checking G4Orb::DistanceToIn(P)

    Dist=s1.DistanceToIn(ponrmax1);
    assert(ApproxEqual(Dist,0));

// Checking G4Orb::DistanceToIn(P,V)

    Dist=s1.DistanceToIn(ponzmax,vz);
    G4cout<<"s1.DistanceToIn(ponzmax,vz) = "<<Dist<<G4endl;
    Dist=s1.DistanceToIn(pbigy,vy);
    assert(Dist==kInfinity);
    Dist=s1.DistanceToIn(pbigy,vmy);
    assert(ApproxEqual(Dist,50));

  // Checking In/Out both return zeros

  Dist = solidD1->DistanceToIn(positionKip, directionKip);
  G4cout << " Dist In  (pos, uKipp) = " << Dist << G4endl;

  Dist = solidD1->DistanceToOut(positionKip, directionKip);
  G4cout << " Dist Out (pos, uKipp) = " << Dist << G4endl;
     
  G4cout << "   Dot( pos, uKipp) = " << positionKip.dot(directionKip)<< G4endl;

  inside = solidD1->Inside(positionKip);
  G4cout<<"solidD1->Inside(positionKip) = " << OutputInside(inside)<<G4endl ;



     ///////////////////////////////////////////////////////////////////////////


// CalculateExtent
    G4VoxelLimits limit;		// Unlimited
    G4RotationMatrix noRot;
    G4AffineTransform origin;
    G4double min,max;
    assert(s1.CalculateExtent(kXAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kYAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kZAxis,limit,origin,min,max));
    assert(min<=-50&&max>=50);

    G4ThreeVector pmxmymz(-100,-110,-120);
    G4AffineTransform tPosOnly(pmxmymz);
    assert(s1.CalculateExtent(kXAxis,limit,tPosOnly,min,max));
    assert(min<=-150&&max>=-50);
    assert(s1.CalculateExtent(kYAxis,limit,tPosOnly,min,max));
    assert(min<=-160&&max>=-60);
    assert(s1.CalculateExtent(kZAxis,limit,tPosOnly,min,max));
    assert(min<=-170&&max>=-70);

    
    G4RotationMatrix r90Z;
    r90Z.rotateZ(halfpi);
    G4AffineTransform tRotZ(r90Z,pzero);
    assert(s1.CalculateExtent(kXAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kYAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);
    assert(s1.CalculateExtent(kZAxis,limit,tRotZ,min,max));
    assert(min<=-50&&max>=50);

    
// Check that clipped away
    G4VoxelLimits xClip;
    xClip.AddLimit(kXAxis,-100,-60);
    assert(!s1.CalculateExtent(kXAxis,xClip,origin,min,max));

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
    assert(s1.CalculateExtent(kXAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s1.CalculateExtent(kYAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);
    assert(s1.CalculateExtent(kZAxis,allClip,tGen,min,max));
    assert(min<=-5&&max>=5);


// Test z clipping ok
    for (G4double zTest=-100;zTest<100;zTest+=9)
    {
      G4VoxelLimits zTestClip;
      zTestClip.AddLimit(kZAxis,-kInfinity,zTest);
      if (zTest<-50)
      {
	assert(!s1.CalculateExtent(kZAxis,zTestClip,origin,min,max));	   
      }
      else
      {
	assert(s1.CalculateExtent(kZAxis,zTestClip,origin,min,max));
		    
	G4double testMin=-50;
        G4double testMax=(zTest<50) ? zTest : 50;
	assert (ApproxEqual(min,testMin) && ApproxEqual(max,testMax));
      }
    }

// Test y clipping ok
    for (G4double xTest=-100;xTest<100;xTest+=9)
    {
      G4VoxelLimits xTestClip;
      xTestClip.AddLimit(kXAxis,-kInfinity,xTest);
      if (xTest<-50)
      {
        assert(!s1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
      }
      else
      {
        assert(s1.CalculateExtent(kYAxis,xTestClip,origin,min,max));
// Calc max y coordinate
	G4double testMax=(xTest<0) ? sqrt(50*50-xTest*xTest) : 50;
	assert (ApproxEqual(min,-testMax) && ApproxEqual(max,testMax));
      }
    }

// Test x clipping ok
    for (G4double yTest=-100;yTest<100;yTest+=9)
    {
      G4VoxelLimits yTestClip;
      yTestClip.AddLimit(kYAxis,-kInfinity,yTest);
      if (yTest<-50)
      {
        assert(!s1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
      }
      else
      {
        assert(s1.CalculateExtent(kXAxis,yTestClip,origin,min,max));
// Calc max y coordinate
	G4double testMax=(yTest<0) ? sqrt(50*50-yTest*yTest) : 50;
	assert (ApproxEqual(min,-testMax) && ApproxEqual(max,testMax));
      }
    }
    return 0;
}











