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
// $Id: testG4ReplicaNavigation.cc,v 1.7 2002-11-18 16:43:06 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// Test private location & distance computation functions of
// G4ReplicaNavigation                      Paul Kent Aug 96


#include <assert.h>
#include "ApproxEqual.hh"
#include "globals.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ReplicaNavigation.hh"
#include "G4PVReplica.hh"
#include "G4PVPlacement.hh"

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

  // Define two worlds
  //
  G4Box* hall_box =
         new G4Box("expHall_box",3000,3000,3000);
  G4LogicalVolume* hall_log1 =
         new G4LogicalVolume(hall_box,0,"expHall_log",0,0,0);
  G4VPhysicalVolume* hall_phys1 =
         new G4PVPlacement(0,G4ThreeVector(),"expHall1",hall_log1,0,false,0);
  G4LogicalVolume* hall_log2 =
         new G4LogicalVolume(hall_box,0,"expHall_log",0,0,0);
  G4VPhysicalVolume* hall_phys2 =
         new G4PVPlacement(0,G4ThreeVector(),"expHall2",hall_log2,0,false,0);

  // Define volumes to be sliced
  //
  G4Box* fBox = new G4Box("Test Box",120.,120.,120.);
  G4LogicalVolume *pMotherVol1X= new G4LogicalVolume(fBox, 0, "lmoth1X", 0, 0, 0);
    new G4PVPlacement(0,G4ThreeVector(),"pmoth1",pMotherVol1X,hall_phys1,false,0);
  G4LogicalVolume *pMotherVol1Y= new G4LogicalVolume(fBox, 0, "lmoth1Y", 0, 0, 0);
    new G4PVPlacement(0,G4ThreeVector(),"pmoth1",pMotherVol1Y,hall_phys1,false,0);
  G4LogicalVolume *pMotherVol1Z= new G4LogicalVolume(fBox, 0, "lmoth1Z", 0, 0, 0);
    new G4PVPlacement(0,G4ThreeVector(),"pmoth1",pMotherVol1Z,hall_phys1,false,0);

  G4Sphere* fSphere = new G4Sphere("Test Sphere",0.,80.,0*deg,360*deg,0*deg,360*deg);
  G4LogicalVolume *pMotherVol2P= new G4LogicalVolume(fSphere, 0, "lmoth2P", 0, 0, 0);
    new G4PVPlacement(0,G4ThreeVector(),"pmoth2",pMotherVol2P,hall_phys2,false,0);
  G4LogicalVolume *pMotherVol2R= new G4LogicalVolume(fSphere, 0, "lmoth2R", 0, 0, 0);
    new G4PVPlacement(0,G4ThreeVector(),"pmoth2",pMotherVol2R,hall_phys2,false,0);

  G4ReplicaNavigationTester repNav;

  // Define the actual slices (cartesian axis)
  //
  G4Box* xBoxSlice = new G4Box("Sliced Box X",40.,120.,120.);
  G4LogicalVolume* xBoxLog= new G4LogicalVolume(xBoxSlice, 0, "xBoxSlice", 0, 0, 0);
  G4PVReplica xRep("TestX",xBoxLog,pMotherVol1X,kXAxis,3,40);

  G4Box* yBoxSlice = new G4Box("Sliced Box Y",120.,40.,120.);
  G4LogicalVolume* yBoxLog= new G4LogicalVolume(yBoxSlice, 0, "yBoxSlice", 0, 0, 0);
  G4PVReplica yRep("TestY",yBoxLog,pMotherVol1Y,kYAxis,3,40);

  G4Box* zBoxSlice = new G4Box("Sliced Box Z",120.,120.,40.);
  G4LogicalVolume* zBoxLog= new G4LogicalVolume(zBoxSlice, 0, "zBoxSlice", 0, 0, 0);
  G4PVReplica zRep("TestZ",zBoxLog,pMotherVol1Z,kZAxis,3,40);

  // Define the actual slices (Phi and Rho)
  //
  G4Sphere* fSphereP = new G4Sphere("Sliced Sphere Phi",0.,80.,0*deg,90*deg,0*deg,360*deg);
  G4LogicalVolume* phiSphereLog= new G4LogicalVolume(fSphereP, 0, "PhiSlice", 0, 0, 0);
  G4PVReplica phiRep("TestPhi",phiSphereLog,pMotherVol2P,kPhi,4,M_PI*0.5);
  G4Sphere* fSphereR = new G4Sphere("Sliced Sphere Rho",0.,20.,0*deg,360*deg,0*deg,360*deg);
  G4LogicalVolume* rhoSphereLog= new G4LogicalVolume(fSphereR, 0, "RhoSlice", 0, 0, 0);
  G4PVReplica radRep("TestRho",rhoSphereLog,pMotherVol2R,kRho,4,20);

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
  G4std::cout.precision(8);
  // G4cout << " Dist is " << Dist << " and expected= " << sqrt(2.*441.)-20. << G4endl;
  // G4cout << "   a difference of " << Dist-(sqrt(2.*441.)-20.) << G4endl;
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

