// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: testG4ReplicaNavigation.cc,v 1.3 1999-12-15 14:50:29 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Test private location & distance computation functions of
// G4ReplicaNavigation                      Paul Kent Aug 96


#include <assert.h>
#include "ApproxEqual.hh"
#include "globals.hh"
#include "G4ReplicaNavigation.hh"
#include "G4PVReplica.hh"

class G4ReplicaNavigationTester
{
public:
  EInside Inside(const G4VPhysicalVolume *pVol,
		 const G4int replicaNo,
		 const G4ThreeVector &localPoint) const
  {
    return nav.Inside(pVol,replicaNo,localPoint);
  }

  G4double DistanceToOut(const G4VPhysicalVolume *pVol,
			 const G4int replicaNo,
			 const G4ThreeVector &localPoint) const
  {
    return nav.DistanceToOut(pVol,replicaNo,localPoint);
  }

  G4double DistanceToOut(const G4VPhysicalVolume *pVol,
			 const G4int replicaNo,
			 const G4ThreeVector &localPoint,
			 const G4ThreeVector &localDirection) const
  {
    return nav.DistanceToOut(pVol,replicaNo,localPoint,localDirection);
  }

private:
  G4ReplicaNavigation nav;
};

G4bool testG4ReplicaNavigation()
{
  EInside in;
  G4double Dist;

  G4LogicalVolume *pMotherVol= 0;

  G4ReplicaNavigationTester repNav;
  G4PVReplica xRep("Test",0,pMotherVol,kXAxis,3,40);
  G4PVReplica yRep("Test",0,pMotherVol,kYAxis,3,40);
  G4PVReplica zRep("Test",0,pMotherVol,kZAxis,3,40);
  G4PVReplica phiRep("Test",0,pMotherVol,kPhi,4,M_PI*0.5);
  G4PVReplica radRep("Test",0,pMotherVol,kRho,4,20);

  in=repNav.Inside(&xRep,0,G4ThreeVector(21,0,0));
  assert(in==kOutside);
  in=repNav.Inside(&xRep,0,G4ThreeVector(20,0,0));
  assert(in==kSurface);
  in=repNav.Inside(&xRep,0,G4ThreeVector(19,0,0));
  assert(in==kInside);
  in=repNav.Inside(&xRep,0,G4ThreeVector(-20,0,0));
  assert(in==kSurface);
  in=repNav.Inside(&xRep,0,G4ThreeVector(-21,0,0));
  assert(in==kOutside);

  in=repNav.Inside(&yRep,0,G4ThreeVector(0,21,0));
  assert(in==kOutside);
  in=repNav.Inside(&yRep,0,G4ThreeVector(0,20,0));
  assert(in==kSurface);
  in=repNav.Inside(&yRep,0,G4ThreeVector(0,19,0));
  assert(in==kInside);
  in=repNav.Inside(&yRep,0,G4ThreeVector(0,-20,0));
  assert(in==kSurface);
  in=repNav.Inside(&yRep,0,G4ThreeVector(0,-21,0));
  assert(in==kOutside);

  in=repNav.Inside(&zRep,0,G4ThreeVector(0,0,21));
  assert(in==kOutside);
  in=repNav.Inside(&zRep,0,G4ThreeVector(0,0,20));
  assert(in==kSurface);
  in=repNav.Inside(&zRep,0,G4ThreeVector(0,0,19));
  assert(in==kInside);
  in=repNav.Inside(&zRep,0,G4ThreeVector(0,0,-20));
  assert(in==kSurface);
  in=repNav.Inside(&zRep,0,G4ThreeVector(0,0,-21));
  assert(in==kOutside);

  in=repNav.Inside(&phiRep,0,G4ThreeVector(0,0,0));
  assert(in==kSurface);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(10,0,0));
  assert(in==kInside);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(-10,0,0));
  assert(in==kOutside);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(10,10,0));
  assert(in==kSurface);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(10,10.1,0));
  assert(in==kOutside);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(10,-10,0));
  assert(in==kSurface);
  in=repNav.Inside(&phiRep,0,G4ThreeVector(10,-10.1,0));
  assert(in==kOutside);

  in=repNav.Inside(&radRep,0,G4ThreeVector(0,0,0));
  assert(in==kInside);
  in=repNav.Inside(&radRep,0,G4ThreeVector(0,20,0));
  assert(in==kSurface);
  in=repNav.Inside(&radRep,0,G4ThreeVector(0,21,0));
  assert(in==kOutside);

  in=repNav.Inside(&radRep,1,G4ThreeVector(0,0,0));
  assert(in==kOutside);
  in=repNav.Inside(&radRep,1,G4ThreeVector(0,20,0));
  assert(in==kSurface);
  in=repNav.Inside(&radRep,1,G4ThreeVector(0,30,0));
  assert(in==kInside);

  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(0,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(20,20,20));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(-21,-21,-21));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&yRep,0,G4ThreeVector(0,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&yRep,0,G4ThreeVector(20,20,20));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&yRep,0,G4ThreeVector(-21,-21,-21));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&zRep,0,G4ThreeVector(0,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&zRep,0,G4ThreeVector(20,20,20));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&zRep,0,G4ThreeVector(-21,-21,-21));
  assert(Dist==0);


  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,0,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0));
  assert(ApproxEqual(Dist,10*sin(M_PI*0.25)));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(-10,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,10,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,-10,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,5,0));
  assert(ApproxEqual(Dist,sqrt(125.)*sin(M_PI*0.25-atan(0.5))));

  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,20,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,21,0));
  assert(Dist==0);

  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(0,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(0,20,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(0,21,0));
  assert(Dist==1);
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(21,21,0));
  cout.precision(8);
  // cout << " Dist is " << Dist << " and expected= " << sqrt(2.*441.)-20. << G4endl;
  // cout << "   a difference of " << Dist-(sqrt(2.*441.)-20.) << G4endl;
  assert( Dist - (sqrt(2.*441.)-20.) < 1.e-14 );
  // assert(ApproxEqual(Dist, sqrt(2.*441.)-20.));


  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(1,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(20,0,0),
			    G4ThreeVector(1,0,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(20,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,40));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(21,0,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(20,0,0),
			    G4ThreeVector(-1/sqrt(2.),-1/sqrt(2.),0));
  assert(ApproxEqual(Dist,40*sqrt(2.)));
  Dist=repNav.DistanceToOut(&xRep,0,G4ThreeVector(20,0,0),
			    G4ThreeVector(0,1,0));
  assert(Dist==kInfinity);


  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==kInfinity);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(-1,0,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==kInfinity);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,-1,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==kInfinity);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,1,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==kInfinity);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(-1,0,0),
			    G4ThreeVector(-1,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,-1,0),
			    G4ThreeVector(-1,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,1,0),
			    G4ThreeVector(-1,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,10));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0),
			    G4ThreeVector(0,1,0));
  assert(ApproxEqual(Dist,10));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0),
			    G4ThreeVector(0,-1,0));
  assert(ApproxEqual(Dist,10));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0),
			    G4ThreeVector(-1/sqrt(2.),1/sqrt(2.),0));
  assert(ApproxEqual(Dist,10*sin(M_PI*0.25)));
  Dist=repNav.DistanceToOut(&phiRep,0,G4ThreeVector(10,0,0),
			    G4ThreeVector(-1/sqrt(2.),-1/sqrt(2.),0));
  assert(ApproxEqual(Dist,10*sin(M_PI*0.25)));


  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(1,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(0,0,0),
			    G4ThreeVector(-1/sqrt(2.),-1/sqrt(2.),0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(sqrt(200.),sqrt(200.),0),
			    G4ThreeVector(-1/sqrt(2.),-1/sqrt(2.),0));
  assert(ApproxEqual(Dist,40));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(sqrt(200.),sqrt(200.),0),
			    G4ThreeVector(1/sqrt(2.),1/sqrt(2.),0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(sqrt(200.),sqrt(200.),0),
			    G4ThreeVector(0,0,1));
  assert(Dist==kInfinity);
  Dist=repNav.DistanceToOut(&radRep,0,G4ThreeVector(21,0,0),
			    G4ThreeVector(1,0,0));
  assert(Dist==0);
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(20,0,0),
			    G4ThreeVector(1,0,0));
  assert(ApproxEqual(Dist,20));
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(20,0,0),
			    G4ThreeVector(-1,0,0));
  assert(ApproxEqual(Dist,0));
  Dist=repNav.DistanceToOut(&radRep,1,G4ThreeVector(20,0,0),
			    G4ThreeVector(0,-1,0));
  assert(ApproxEqual(Dist,sqrt(40.*40.-20.*20.)));


  return true;
}

int main()
{
#ifdef NDEBUG
    G4Exception("FAIL: *** Assertions must be compiled in! ***");
#endif
    assert(testG4ReplicaNavigation());
    return 0;
}

