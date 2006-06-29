//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: testG4ParameterisedSolid1.cc,v 1.13 2006-06-29 18:58:43 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
//   Define geometry with parameterised volumes that parameterise solid type
//
//   Test the Navigation in this geometry 
//                  (which also include rotations as well as translations).
//

#include <assert.h>
#include "G4ios.hh"
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"
#include "G4Sphere.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Sample Parameterisation
class BoxesAndSpheres : public G4VPVParameterisation
{
 public:
  BoxesAndSpheres(G4double twistAngle, G4int noBoxes, G4int noSpheres)
  { 
    fRotationVec= new G4RotationMatrix();
    fTwistAngle= twistAngle;
    fNumBoxes= noBoxes;
    fNumSpheres= noSpheres;
    fBox= new G4Box("Test Box",10.,10.,10.);
    fSphere= new G4Sphere("Test Sphere",0.,1.,0*deg,180*deg,0*deg,90*deg);
  }

  virtual ~BoxesAndSpheres()
  {
//    delete fRotationVec;
//    delete fBox;
//    delete fSphere;
  }

  G4double GetTwistAngle() { return fTwistAngle; }
  void     SetTwistAngle(G4double newAngle ) { fTwistAngle= newAngle; }

  virtual G4VSolid* ComputeSolid(const G4int n,
				 G4VPhysicalVolume*) 
  {
    G4VSolid* mySolid=0;
    if( n < fNumBoxes ) {
       if( n >= 0 ) {
          mySolid = fBox;
       }else{
	  G4Exception(" Your Boxes+Spheres replica number was out of range");
       }
    }else{
       if( n < fNumBoxes + fNumSpheres ) {
          mySolid = fSphere;
       }else{
	  G4Exception(" Your Boxes+Spheres replica number was out of range");
       }
    }
    return mySolid;
  }				 
				     
  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(n*100*mm,0.,0.));
    *fRotationVec = G4RotationMatrix();             // Unit matrix
    fRotationVec->rotateZ( n * fTwistAngle );
    pRep->SetRotation( fRotationVec );
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int,
				 const G4VPhysicalVolume*) const
  {
    if( &pBox != fBox ){
      G4cerr << " Got another Box in ComputeDimensions(G4Box, , )" << G4endl;
    }
    pBox.SetXHalfLength(10*mm);
    pBox.SetYHalfLength(10*mm);
    pBox.SetZHalfLength(10*mm);
  }
  
  virtual void ComputeDimensions(G4Sphere &pSphere,
				 const G4int n,
				 const G4VPhysicalVolume*) const
  {
    if( &pSphere != fSphere )
    {
      G4cerr << " Got another sphere in ComputeDimensions(G4Sphere, , )"
             << G4endl;
    }
    G4int nrad= std::min(5, n-fNumBoxes+1);
    pSphere.SetInsideRadius( nrad *  5. * mm);
    pSphere.SetOuterRadius ( nrad * 10. * mm);
    pSphere.SetStartPhiAngle  (0.);
    pSphere.SetDeltaPhiAngle  (2*pi);
    pSphere.SetStartThetaAngle(0);
    pSphere.SetDeltaThetaAngle(pi);
  }
 
  virtual void ComputeDimensions(G4Tubs &,
				 const G4int ,
                                 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trd &, 
				 const G4int,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Cons &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trap &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Hype &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Orb &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Torus &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Para &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polycone &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polyhedra &,
				 const G4int ,
				 const G4VPhysicalVolume*) const {}
 private:
    G4RotationMatrix *fRotationVec;
    G4double  fTwistAngle;
    G4int     fNumBoxes;
    G4int     fNumSpheres;
    G4Box*    fBox;
    G4Sphere* fSphere;
};

G4double    angle1= 15.0*deg;          // pi/180. ;
BoxesAndSpheres myParam(angle1,3,4);

// Build simple geometry:
// 4 small cubes (G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    // The world volume
    //
    G4Box *myBigBox= new G4Box ("Big Cube", 1000*mm, 1000*mm, 1000*mm);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myBigBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set


    // A set of boxes
    G4Box *myBox=new G4Box("cube",10,10,10);
    G4LogicalVolume *boxLog=new G4LogicalVolume(myBox,0,
						"Rotating Box",0,0,0);

//  G4PVParameterised *paramP=
                              new G4PVParameterised("Rotating Block Or Sphere",
						    boxLog,
						    worldPhys, //OR worldLog,
						    kXAxis,
						    7,
						    &myParam);
    // Copies 0, 1 & 2 will exist    

    return worldPhys;
}

//
// Test LocateGlobalPointAndSetup
//
G4bool testG4Navigator1(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(100,100,100),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// 
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-5,-5),0,false);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetCopyNo()== 0);
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),G4ThreeVector(0,-5,-5)));
    G4cout << " Local coords = " << myNav.CurrentLocalCoordinate() << G4endl;

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(100,0,5));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetCopyNo()== 1);
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),
                       G4ThreeVector(0. ,0., 5)));
    
    // Check that the rotation is correct
    //
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(105,0,0));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetCopyNo()== 1);
#if 0 
//  G4cout << " Local coords = " << myNav.GetCurrentLocalCoordinate() << G4endl;
    G4ThreeVector ExpectedPosition(5*std::cos(angle1),-5.*std::sin(angle1),0.);
    G4cout << " Expected     = " << ExpectedPosition << G4endl;
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),ExpectedPosition ));
    if(!ApproxEqual(myNav.CurrentLocalCoordinate(),ExpectedPosition ))
       {
          G4cout << " Error: The coordinates do not match " << G4endl;
       }
#endif 
    
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check parameterised volumes

// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(5,0,5));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 0);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,15,15));
    assert(located->GetName()=="World");

// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(105,0,5));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 1);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-17));
    assert(located->GetName()=="World");

// Replication 2
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(205,0,5));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 2);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-18));
    assert(located->GetName()=="World");
    
// Replication 3
    // Sphere 1,  radii: inner/outer=  5 to 10
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(307.5,0.0,0.0));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 3);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(315,15,-18));
    assert(located->GetName()=="World");

// Replication 4
    // Sphere 2,  radii: inner/outer= 10 to 20
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(410.,10.,10.));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 4);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(315,15,-18));
    assert(located->GetName()=="World");

// Replication 5
    // Sphere 3,  radii: inner/outer= 15 to 30
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(510.0,10.0,10.0));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 5);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(500,35,-10));
    assert(located->GetName()=="World");
    
// Replication 6
    // Sphere 4,  radii: inner/outer= 20 to 40
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(622.5,22.5,22.5));
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 6);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(600,45,-10));
    assert(located->GetName()=="World");

    return true;
}


//
// Test Stepping
//
G4bool testG4Navigator2(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4VPhysicalVolume *located;
    G4double Step,physStep,safety;
    G4ThreeVector xHat(1,0,0),yHat(0,1,0),zHat(0,0,1);
    G4ThreeVector mxHat(-1,0,0),myHat(0,-1,0),mzHat(0,0,-1);
    
    myNav.SetWorldVolume(pTopNode);
  
//
// Test location & Step computation
//  
    G4ThreeVector  StartPoint(-50*mm,0,-5*mm);
    located=myNav.LocateGlobalPointAndSetup( StartPoint ); 
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( StartPoint, mxHat,physStep,safety);  // -x dir
    assert(ApproxEqual(Step,950*mm));
    // assert(ApproxEqual(safety,40*mm));
    // assert(safety>=0);

    StartPoint= G4ThreeVector(-15*mm,0,-5*mm);
    located=myNav.LocateGlobalPointAndSetup( StartPoint ); 
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( StartPoint,xHat,physStep,safety); // +x dir
    assert(ApproxEqual(Step,5));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);
    myNav.SetGeometricallyLimitedStep();
    G4ThreeVector EndPoint = StartPoint + Step * xHat;
    located=myNav.LocateGlobalPointAndSetup(EndPoint,0,true);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 0);


    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-40));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-40),zHat,physStep,safety); // +z
    assert(ApproxEqual(Step,30));
    assert(safety>=0);
    // Now locate the endpoint
    myNav.SetGeometricallyLimitedStep();
    EndPoint = StartPoint + Step * xHat;
    located=myNav.LocateGlobalPointAndSetup(EndPoint,0,true);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0, 40));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,40),mzHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);


//
// Test moving through series of volumes
//
    StartPoint= G4ThreeVector(-20,0,0);
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-20,0,0));
    assert(located->GetName()=="World");
    
    // Replication 0 block
    //
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,10));
    EndPoint= StartPoint + Step * xHat;   //  Should be  -10, 0, 0
    assert( ApproxEqual( EndPoint, G4ThreeVector(-10,0,0) ) );
    assert(safety<=10);
    
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint) ;
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 0);
    Step=myNav.ComputeStep(EndPoint,xHat,physStep,safety); // +x
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * xHat;   //  Should be   +10, 0, 0
    assert(ApproxEqual( EndPoint, G4ThreeVector(10,0,0) ));
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Replication 1 block
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,90.-10./std::cos(angle1)));
    EndPoint= StartPoint + Step * xHat;   //  Should be near  90, 0, 0
    assert(safety==0.);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint) ;
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 1);
    
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,20./std::cos(angle1)));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * xHat;   //  Should be near 110, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(100.+10./std::cos(angle1),0,0)));
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Replication 2 block
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,100.-10.*(1./std::cos(angle1)+1./std::cos(2.*angle1))));
    EndPoint= StartPoint + Step * xHat;   //  Should near  0, 190, 0
    assert(safety<=Step);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Box");
    assert(located->GetCopyNo()== 2);
    
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,20./std::cos(2.*angle1)));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * xHat;   //  Should be near 210, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(200.+10./std::cos(2.*angle1),0,0)));
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");


    // Replication 3 :  sphere #1
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,(100.-10./std::cos(2.*angle1)-10.)*mm));
    EndPoint= StartPoint + Step * xHat;   //  Should be    290, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(290.*mm,0,0)));
    assert(safety==0.); // Started from a surface
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 3);
    
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,5.));
    assert(ApproxEqual(safety,0.));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * xHat;   //  Should be near 295, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(295*mm,0,0)));
    // Now Hit inner surface of spherical shell
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");
    
    // Cross "empty" inner sphere
    //
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,(10.*mm)));
    EndPoint= StartPoint + Step * xHat;   //  Should be    290, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(305.*mm,0,0)));
    assert(ApproxEqual(safety,0.)); // Started from a surface
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(EndPoint);
    assert(located->GetName()=="Rotating Block Or Sphere");
    assert(located->GetLogicalVolume()->GetSolid()->GetName()=="Test Sphere");
    assert(located->GetCopyNo()== 3);
    
    // Now Hit outer surface of spherical shell
    StartPoint= EndPoint;
    physStep=kInfinity;
    Step=myNav.ComputeStep(StartPoint,xHat,physStep,safety);
    assert(ApproxEqual(Step,5.));
    assert(ApproxEqual(safety,0.));
    myNav.SetGeometricallyLimitedStep();
    EndPoint += Step * xHat;   //  Should be near 310, 0, 0
    assert(ApproxEqual(EndPoint,G4ThreeVector(310*mm,0,0)));
    located=myNav.LocateGlobalPointAndSetup( EndPoint );
    assert(located->GetName()=="World");

    // Continue the test later ...

    return true;
}

int main()
{
    G4VPhysicalVolume *myTopNode;
    myTopNode=BuildGeometry();	// Build the geometry
    G4GeometryManager::GetInstance()->CloseGeometry(false);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4Navigator1(myTopNode);
    testG4Navigator2(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}




