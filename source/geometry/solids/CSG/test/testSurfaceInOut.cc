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
// Test File for tracking functions on a solid surface
//
// o Basic asserts on each function +
//   awkward cases for tracking / geom algorithms
//
// o Add tests on dicovering bugs in G4Sphere.cc...
//
// History:
//
// 14.07.04 V.Grichine creation for sphere

#include "G4ios.hh"
#include <assert.h>
#include <math.h>
#include "globals.hh"
#include "geomdefs.hh"
#include "Randomize.hh"

#include "ApproxEqual.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"

#include "G4Box.hh"
#include "G4Orb.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4Cons.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4Torus.hh"
#include "G4Trd.hh"



///////////////////////////////////////////////////////////////////////////
//
//







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


/////////////////////////////////////////////////////////////////////////////
//
// Random unit direction vector

G4ThreeVector GetRandomUnitVector()
{
  G4double cosTheta, sinTheta, phi, vx, vy, vz;

  cosTheta = -1. + 2.*G4UniformRand();
  if( cosTheta > 1.)  cosTheta = 1.;
  if( cosTheta < -1.) cosTheta = -1.;
  sinTheta = sqrt( 1. - cosTheta*cosTheta );
  
  phi = 2*pi*G4UniformRand();

  vx = sinTheta*cos(phi);
  vy = sinTheta*sin(phi);
  vz = cosTheta;

  return G4ThreeVector(vx,vy,vz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random  vector on box surface 

G4ThreeVector GetVectorOnBox( G4Box& box )
{
  G4double rand, a, b, c, px, py, pz;
  G4double part = 1./6.;

  a = box.GetXHalfLength();
  b = box.GetYHalfLength();
  c = box.GetZHalfLength();
  
  rand = G4UniformRand();

  if      ( rand < part )
  {
    px = -a - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 2*part )
  {
    px =  a - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 3*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 4*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py =  b - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (2.*c + kCarTolerance)*G4UniformRand();  
  }
  else if ( rand < 5*part )
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz = -c - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();  
  }
  else 
  {
    px = -a - 0.5*kCarTolerance + (2.*a + kCarTolerance)*G4UniformRand(); 
    py = -b - 0.5*kCarTolerance + (2.*b + kCarTolerance)*G4UniformRand(); 
    pz =  c - 0.5*kCarTolerance + (kCarTolerance)*G4UniformRand();  
  }

  return G4ThreeVector(px,py,pz);
}

/////////////////////////////////////////////////////////////////////////////
//
// Random vector on orb surface

G4ThreeVector GetVectorOnOrb(G4Orb& orb)
{
  G4double cosTheta, sinTheta, phi, radius, px, py, pz;

  radius = orb.GetRadius();
  radius += -0.5*kCarTolerance + (kCarTolerance)*G4UniformRand(); 

  cosTheta = -1. + 2.*G4UniformRand();
  if( cosTheta > 1.)  cosTheta = 1.;
  if( cosTheta < -1.) cosTheta = -1.;
  sinTheta = sqrt( 1. - cosTheta*cosTheta );
  
  phi = 2*pi*G4UniformRand();

  px = radius*sinTheta*cos(phi);
  py = radius*sinTheta*sin(phi);
  pz = radius*cosTheta;

  return G4ThreeVector(px,py,pz);
}


//////////////////////////////////////////////////////////////////////
//
// Main executable function

int main(void)
{
  G4double distIn, distOut;
  EInside surfaceP;
  G4ThreeVector norm,*pNorm;
  G4bool *pgoodNorm, goodNorm, calcNorm=true;

  enum Esolid {kBox, kOrb, kSphere, kCons, kTubs, kTorus, kPara, kTrapezoid, kTrd};

  Esolid useCase = kOrb;

  pNorm=&norm;
  pgoodNorm=&goodNorm;

  // Boxes for test

  G4Box b1("Test Box #1",20,30,40);
  G4Box b2("Test Box #2",10,10,10);
  G4Box box3("BABAR Box",0.14999999999999999, 
                           24.707000000000001,  
                           22.699999999999999) ;

  // orbs for test

  G4Orb  o1("Solid G4Orb",50);
  G4Orb  o10("s10",0.018*mm);
  G4Orb* solidO1= new G4Orb("O1", 2.7*cm);


  // spheres for test

  G4Sphere s1("Solid G4Sphere",0,50,0,twopi,0,pi);
  G4Sphere s2("Spherical Shell",45,50,0,twopi,0,pi);
  G4Sphere s3("Band (theta segment)",45,50,0,twopi,pi/4,halfpi);
  G4Sphere s32("Band (theta segment2)",45,50,0,twopi,0,pi/4);
  G4Sphere s33("Band (theta segment1)",45,50,0,twopi,pi*3/4,pi/4);
  G4Sphere s34("Band (theta segment)",4,50,0,twopi,pi/4,halfpi);
  G4Sphere s4("Band (phi segment)",45,50,-pi/4,halfpi,0,twopi);
  //    G4cout<<"s4.fSPhi = "<<s4.GetSPhi()<<G4endl;
  G4Sphere s41("Band (phi segment)",5,50,-pi,3.*pi/2.,0,twopi);
  G4Sphere s42("Band (phi segment)",5,50,-pi/2,3.*pi/2.,0,twopi);
  G4Sphere s5("Patch (phi/theta seg)",45,50,-pi/4,halfpi,pi/4,halfpi);
  G4Sphere s6("John example",300,500,0,5.76,0,pi) ; 
  G4Sphere s7("sphere7",1400.,1550.,0.022321428571428572,0.014642857142857141,
	                  1.5631177553663251,0.014642857142857141    );
  G4Sphere s8("sphere",278.746*mm, 280.0*mm, 0.0*degree, 360.0*degree,
		                               0.0*degree, 90.0*degree);

  G4Sphere b216("b216", 1400.0, 1550.0, 
                  0.022321428571428572, 
		  0.014642857142857141,
                  1.578117755366325,
                  0.014642857142867141);

  G4Sphere s9("s9",0*mm,410*mm,0*degree,360*degree,90*degree,90*degree);

  G4Sphere b402("b402", 475*mm, 480*mm, 
               0*degree,360*degree,17.8*degree,144.4*degree);
   

  G4Sphere s10("s10",0*mm,0.018*mm,0*degree,360*degree,0*degree,180*degree);


  G4Sphere s11("s11",5000000.*mm,
                    3700000000.*mm,
                   0*degree,360*degree,0*degree,180*degree);


  G4Sphere sAlex("sAlex",500.*mm,
                    501.*mm,
                   0*degree,360*degree,0*degree,180*degree);

  G4Sphere sLHCB("sLHCB",8600*mm, 8606*mm, 
    -1.699135525184141*degree,
    3.398271050368263*degree,88.52855940538514*degree,2.942881189229715*degree );

  G4Sphere spAroundX("SpAroundX",  10.*mm, 1000.*mm, -1.0*degree, 
                                                     2.0*degree, 
                                        0.*degree, 180.0*degree );


#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif

  G4cout.precision(20);

  G4int i,j, iMax=10000, jMax=1000;
  
  // Check box tracking function 

  switch (useCase)
  {

    case kBox:
    G4cout<<"Testing G4Box:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      G4ThreeVector p = GetVectorOnBox(b1);
    
      surfaceP = b1.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = b1.DistanceToIn(p,v);
          distOut = b1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

          if( distIn < kCarTolerance && distOut < kCarTolerance )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }     
        }
      }
    }
    break;
    
    case kOrb:
      G4cout<<"Testing G4Orb:"<<G4endl<<G4endl;
    for(i=0;i<iMax;i++)
    {
      G4ThreeVector p = GetVectorOnOrb(o1);
    
      surfaceP = o1.Inside(p);
      if(surfaceP != kSurface)
      {
        G4cout<<"p is out of surface: "<<G4endl;
        G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl<<G4endl;

      }
      else
      {
        for(j=0;j<jMax;j++)
        {
          G4ThreeVector v = GetRandomUnitVector();

          distIn  = o1.DistanceToIn(p,v);
          distOut = o1.DistanceToOut(p,v,calcNorm,pgoodNorm,pNorm); 

          if( distIn < kCarTolerance && distOut < kCarTolerance )
	  {
	    G4cout<<" distIn < kCarTolerance && distOut < kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }
          else if(distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance)
	  {
	    G4cout<<" distIn > 100000*kCarTolerance && distOut > 100*kCarTolerance"<<G4endl;
            G4cout<<"distIn = "<<distIn<<";  distOut = "<<distOut<<G4endl;
            G4cout<<"location p: "<<G4endl;
            G4cout<<"( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "<<G4endl;
            G4cout<<" direction v: "<<G4endl;
            G4cout<<"( "<<v.x()<<", "<<v.y()<<", "<<v.z()<<" ); "<<G4endl<<G4endl;
	  }     
        }
      }
    }
    break;

    default:
    break;   
  }
  return 0;
}

