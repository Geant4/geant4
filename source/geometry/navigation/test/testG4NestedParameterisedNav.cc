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
// $Id: testG4NestedParameterisedNav.cc,v 1.1 2005-02-18 18:23:16 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//   Locate & Step within simple boxlike geometry, both
//   with and without voxels. Parameterised volumes are included.

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4VNestedParameterisation.hh"

#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

// Sample First level Parameterisation -- host to nested 2nd
class XTopParam: public G4VPVParameterisation
{
 public: 
  XTopParam( G4int numRowsX, G4double xFullWidth,  G4double yFullWidth, G4double zFullWidth) 
    : fNumRows(numRowsX), 
      fXfullWidth(xFullWidth), 
      fYfullWidth(yFullWidth),
      fZfullWidth(zFullWidth) {} ; 

  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector( (n-((fNumRows+1.0)/2.))*fXfullWidth, 0., 0.) );
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int n,
				 const G4VPhysicalVolume*,
				 const G4VTouchable * ) const
  {
    pBox.SetXHalfLength(fXfullWidth*0.5);
    pBox.SetYHalfLength(fYfullWidth*0.5);
    pBox.SetZHalfLength(fZfullWidth*0.5);
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
  virtual void ComputeDimensions(G4Sphere &,
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
    G4int    fNumRows;
    G4double fXfullWidth, fYfullWidth, fZfullWidth; 

} ;

// Sample Nested Parameterisation
class YSecondNestedParam: public G4VNestedParameterisation
{
  // 
  //   This parameterisation is nested inside another
  //   It creates boxes in a checker-board manner
  //    with different sizes on the odd-even diagonals.
  // 
public:
  YSecondNestedParam( G4int numCols, G4double fullBoxSize, G4double smallBoxSize ) 
   : fNumCols( numCols ) , 
     fFullBoxWidth(fullBoxSize),
     fFullBoxHalfWidth(fullBoxSize*0.5), 
     fSmallBoxHalfWidth(smallBoxSize*0.5)
   {}

  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep,
				     const G4VTouchable * ) const
  {
    pRep->SetTranslation(G4ThreeVector(0., (n-((fNumCols+1.0)/2.))*fFullBoxWidth, 0.) );
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int no_lev,
				 const G4VPhysicalVolume*,
				 const G4VTouchable *parentTouch) const
  {
    // Get the information about the parent volume
    G4int no_parent= parentTouch->GetCopyNumber(); 

    // Rule: Odd ones are small, even ones are full size
    G4int num, odd; 
    num= no_lev + no_parent;
    odd= ( num % 2 ); 

    G4double boxHalfWidth= 0.0; 
    if( odd == 1 ) { 
      boxHalfWidth= fSmallBoxHalfWidth; 
    } else {
      boxHalfWidth= fFullBoxHalfWidth; 
    } 
    pBox.SetXHalfLength(boxHalfWidth);
    pBox.SetYHalfLength(boxHalfWidth);
    pBox.SetZHalfLength(boxHalfWidth);
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
  virtual void ComputeDimensions(G4Sphere &,
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
  G4int    fNumCols; 
  G4double fFullBoxHalfWidth, fFullBoxWidth;
  G4double fSmallBoxHalfWidth; 

} ; // level2NestedParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myWorldBox= new G4Box ("WorldBox",1000.*cm,1000.*cm,1000.*cm);
    G4Box *myTopBox=new G4Box("cube",100.*cm,100.*cm,100.*cm);
    G4Box *mySlab= new G4Box("slab",10.0*cm,100.*cm,100.*cm);
    G4Box *myVariableBox=new G4Box("Variable Box",10.*cm,10.*cm,10.*cm);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myWorldBox,0,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume *topLog=new G4LogicalVolume(myTopBox,0,
						 "Level0 Top-LV"); // ,0,0,0);
    G4LogicalVolume *slabLog=new G4LogicalVolume(mySlab,0,
						 "Level1 Slab-LV"); // ,0,0,0);
    G4LogicalVolume *variLog=new G4LogicalVolume(myVariableBox,0,
						"Level2 Variable Box-LV"); // ,0,0,0);

 
    // Place two 'Top' Boxes
    new G4PVPlacement(0, G4ThreeVector(-250.*cm, 0., 0.),
		      "Target 1", topLog, worldPhys, false, 0);

    new G4PVPlacement(0, G4ThreeVector( 250.*cm, 0., 0.), 
		      "Target 2", topLog, worldPhys, false, 0);


    // Place slabs inside Top Box
    // --------------------------

    G4int    numSlabs= 10; 
    G4double xTopFullWidth=100.*cm, yTopFullWidth=100.*cm, zTopFullWidth=100.*cm; 

    G4double xSlabFullWidth= xTopFullWidth / numSlabs; 
    XTopParam* pFirstLevelParam = 
      new XTopParam( numSlabs, xSlabFullWidth, yTopFullWidth, zTopFullWidth ); 
 
   // myFirstLevParam = new XTopParam( 10, 10.*cm, 100.*cm, 100.*cm ); 

    // G4PVParameterised *paramLevelOnePhys=
                                 new G4PVParameterised("Slab Blocks in X",
						       slabLog,
						       topLog,
						       kXAxis,
						       numSlabs,
						       pFirstLevelParam);

    // Place inner-boxes inside Slabs Box
    // ----------------------------------

    G4int    numBoxesY= 10; 

    G4double ySlabFullWidth= yTopFullWidth;
    G4double zSlabFullWidth= zTopFullWidth;
    // G4double xBoxFullWidth==100.*cm, yBoxFullWidth=100.*cm, zBoxFullWidth=100.*cm; 

    G4double xBoxFullWidth= xSlabFullWidth; 
    G4double yBoxFullWidth= ySlabFullWidth / numBoxesY; 
    G4double zBoxFullWidth= zSlabFullWidth;

    G4VNestedParameterisation* pSecondLevParam =  
      new YSecondNestedParam( numBoxesY, zBoxFullWidth, 0.9*zBoxFullWidth ); 

    // G4PVParameterised *paramLevelTwoPhys=
                                 new G4PVParameterised("Slab Blocks in X",
						       topLog,
						       variLog,
						       kYAxis,
						       numBoxesY,
						       pSecondLevelParam);

    G4cout << " Slab dimensions (full-width) are: " << G4endl
	   << " x= " << xSlabFullWidth/cm << " cm  "
	   << " y= " << ySlabFullWidth/cm << " cm  "
	   << " z= " << ySlabFullWidth/cm << " cm  " << G4endl << G4endl; 

    G4cout << " Box dimensions (full-width) are: " << G4endl
	   << " x= " << xBoxFullWidth/cm << " cm  "
	   << " y= " << yBoxFullWidth/cm << " cm  "
	   << " z= " << yBoxFullWidth/cm << " cm  " << G4endl << G4endl; 

  
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
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check relative search that causes backup one level and then search down:
// Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10),0,false);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,-15,20));
    assert(located->GetName()=="Target 5");
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),G4ThreeVector(0,0,10)));
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check parameterised volumes

// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-15,-16));
    assert(located->GetName()=="Target 3");

// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,0,-17));
    assert(located->GetName()=="Target 3");

// Replication 2
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-10));
    assert(located->GetName()=="Vari' Blocks");
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,15,-18));
    assert(located->GetName()=="Target 3");

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
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mxHat,physStep,safety);
    assert(ApproxEqual(Step,25));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),xHat,physStep,safety);
    assert(ApproxEqual(Step,5));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(5,0,-10),0,true);
    assert(located->GetName()=="Vari' Blocks");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),zHat,physStep,safety);
    assert(ApproxEqual(Step,30));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,-10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(0,0,-10),mzHat,physStep,safety);
    assert(ApproxEqual(Step,10));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);


//
// Test stepping through common boundaries
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,-20));
    assert(located->GetName()=="Target 1");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,-20),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,0));
    assert(located->GetName()=="Target 4");
    Step=myNav.ComputeStep(G4ThreeVector(-7,7,0),zHat,physStep,safety);
    assert(ApproxEqual(Step,20));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-7,7,20));
    assert(!located);

//
// Test mother limited Step
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-25,0,10));
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(-25,0,10),xHat,physStep,safety);
    assert(ApproxEqual(Step,50));
    assert(ApproxEqual(safety,0));

//
// Test stepping through parameterised volumes
//
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-25,-10),0,false);
    assert(located->GetName()=="Target 3");
    physStep=kInfinity;
    Step=myNav.ComputeStep(G4ThreeVector(15,-25,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,5));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-20,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-20,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,10));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-10,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,-10,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,4));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,-6,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,-6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,12));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,6,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,6,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,2));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,8,-10));
    assert(located->GetName()=="Vari' Blocks");
    Step=myNav.ComputeStep(G4ThreeVector(15,8,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,14));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,22,-10));
    assert(located->GetName()=="Target 3");
    Step=myNav.ComputeStep(G4ThreeVector(15,22,-10),yHat,physStep,safety);
    assert(ApproxEqual(Step,3));
    assert(ApproxEqual(safety,0));
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(15,25,-10));
    assert(!located);

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
