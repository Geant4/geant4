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
// $Id: testG4NestedParameterisedNav.cc,v 1.5 2005-11-24 18:01:11 japost Exp $
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

#include "G4Material.hh"
#include "G4Element.hh"

G4Material *darkMaterial, *brightMaterial, *defaultMaterial;  // Chessboard

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
    // G4cout << " Transf for n= " << n << " Offset x= " 
    //   << (n-((fNumRows-1.0)/2.))*fXfullWidth << G4endl;
    pRep->SetTranslation(G4ThreeVector( (n-((fNumRows-1.0)/2.))*fXfullWidth, 0., 0.) );
  }
  
  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int,
				 const G4VPhysicalVolume*) const
  {
    pBox.SetXHalfLength(fXfullWidth*0.5);
    pBox.SetYHalfLength(fYfullWidth*0.5);
    pBox.SetZHalfLength(fZfullWidth*0.5);
  }

  virtual void ComputeDimensions(G4Tubs &, const G4int ,
                                 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trd &,  const G4int,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Cons &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Trap &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Hype &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Orb &,  const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Sphere &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Torus &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Para &,  const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polycone &, const G4int ,
				 const G4VPhysicalVolume*) const {}
  virtual void ComputeDimensions(G4Polyhedra &, const G4int ,
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
  YSecondNestedParam( G4int numCols, G4double ySliceHalfSize, G4double xBoxHalfWidth, G4double zBoxHalfWidth ) 
   : fNumCols( numCols ) , 
     fYBoxHalfWidth(ySliceHalfSize), 
     fYFullBoxWidth(ySliceHalfSize*2.0),
     fXBoxHalfWidth(xBoxHalfWidth), 
     fZBoxHalfWidth(zBoxHalfWidth)
   {}

  virtual void ComputeTransformation(const G4int n,
				     G4VPhysicalVolume* pRep) const
  {
    pRep->SetTranslation(G4ThreeVector(0., 
				       (n-((fNumCols-1)/2.))*fYFullBoxWidth, 
				       0.) );
    pRep->SetRotation(0);  
 }
  
  virtual G4Material* ComputeMaterial(G4VPhysicalVolume *currentVol,
				      const G4int no_lev, 
				      const G4VTouchable *parentTouch) 
  {
    G4Material *material;

    if( parentTouch == 0) {
      G4Exception( "YSecondNestedParam::ComputeMaterial()", 
		   "Null parent TouchHist",  FatalException,
		   " Null pointer as parent touchable pointer. " );
    }

    // Get the information about the parent volume
    G4int no_parent= parentTouch->GetReplicaNumber(); 

    // Rule: Odd ones are one material, even ones are another
    G4int num, odd; 
    num= no_lev + no_parent;
    odd= ( num % 2 ); 

    if( odd == 1 ) { 
      material= darkMaterial;
    } else {
      material= brightMaterial;
    } 
    G4LogicalVolume* currentLogVol= currentVol->GetLogicalVolume(); 
    currentLogVol->SetMaterial( material ); 

    return material;
  }

  G4int       GetNumberOfMaterials() const  { return 2; }

  G4Material* GetMaterial(G4int idx) const 
  {
    G4Material *mat;
    if (idx % 2 == 0){ 
      mat= darkMaterial;
    }else{
      mat= brightMaterial;
    }
    return mat; 
 }

  virtual void ComputeDimensions(G4Box &pBox,
				 const G4int,
				 const G4VPhysicalVolume*) const
  {
    pBox.SetXHalfLength(fXBoxHalfWidth);
    pBox.SetYHalfLength(fYBoxHalfWidth);
    pBox.SetZHalfLength(fZBoxHalfWidth);
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
  G4double fYBoxHalfWidth, fYFullBoxWidth;
  G4double fXBoxHalfWidth, fZBoxHalfWidth; 

} ; // level2NestedParam;

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{
    // Materials 
    // --------------------------------
    //   for use in world and parameterisation 
    G4double a, fractionmass, density; 
    G4int  z, ncomponents; 

    G4Element* N = new G4Element("Nitrogen", "N", z=7, a= 14.01*g/mole);
    G4Element* O = new G4Element("Oxygen"  , "O", z=8, a= 16.00*g/mole);

    G4Material* Air = 
      new G4Material("Air"  , density= 1.290*mg/cm3, ncomponents=2);
    Air->AddElement(N, fractionmass=0.7);
    Air->AddElement(O, fractionmass=0.3);

    //Lead
    G4Material* Pb = 
      new G4Material("Lead", z=82, a= 207.19*g/mole, density= 11.35*g/cm3);
    G4Material* Al = 
     new G4Material("Aluminium", z=13, a=26.98*g/mole, density=2.700*g/cm3);

    // Define standard materials
    darkMaterial= Pb;
    brightMaterial= Al;
    defaultMaterial= Air; 

    // Solids
    // --------------------------------
    G4Box *myWorldBox= new G4Box ("WorldBox",1000.*cm,1000.*cm,1000.*cm);
    G4Box *myTopBox=new G4Box("cube",100.*cm,100.*cm,100.*cm);

    G4LogicalVolume *worldLog=new G4LogicalVolume(myWorldBox,defaultMaterial,
						  "World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
    
    G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					       "World",worldLog,
					       0,false,0);
				// Note: no mother pointer set

    G4LogicalVolume *topLog=new G4LogicalVolume(myTopBox,defaultMaterial,
						 "Level0 Top-LV"); // ,0,0,0);

 
    // Place two 'Top' Boxes in world
    // ------------------------------
    new G4PVPlacement(0, G4ThreeVector(-250.*cm, 0., 0.),
		      "Top 1-pv", topLog, worldPhys, false, 0);

    new G4PVPlacement(0, G4ThreeVector( 250.*cm, 0., 0.), 
		      "Top 2-pv", topLog, worldPhys, false, 1);


    // Place slabs inside Top Box
    // --------------------------
    G4int    numSlabs= 10; 
    G4double xTopHalfWidth=100.*cm, yTopHalfWidth=100.*cm, zTopHalfWidth=100.*cm; 
    G4double xSlabHalfWidth= xTopHalfWidth / numSlabs; 
    G4double ySlabHalfWidth= yTopHalfWidth;
    G4double zSlabHalfWidth= zTopHalfWidth;

    G4Box *mySlab= new G4Box("slab", xSlabHalfWidth, yTopHalfWidth, zTopHalfWidth);
                           // Original: 10.0*cm, 100.*cm, 100.*cm);
    G4LogicalVolume *slabLog=new G4LogicalVolume(mySlab,defaultMaterial,
						 "Level1 Slab-LV"); // ,0,0,0);

    XTopParam* pFirstLevelParam = 
      new XTopParam( numSlabs, xSlabHalfWidth*2., yTopHalfWidth*2., zTopHalfWidth*2. ); 
 
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

    // G4double xBoxHalfWidth==100.*cm, yBoxHalfWidth=100.*cm, zBoxHalfWidth=100.*cm; 

    G4double xBoxHalfWidth= xSlabHalfWidth; 
    G4double yBoxHalfWidth= ySlabHalfWidth / numBoxesY; 
    G4double zBoxHalfWidth= zSlabHalfWidth;
    G4Box *mySmallestBox=new G4Box("Smallest Box", 
				   // 10.*cm, 10.*cm, 100.*cm);
				   xBoxHalfWidth, yBoxHalfWidth, zBoxHalfWidth); 
    G4LogicalVolume *variLog=new G4LogicalVolume(mySmallestBox,defaultMaterial,
						"Level2 Smallest Box-LV"); 
    G4VNestedParameterisation* pSecondLevelParam =  
      new YSecondNestedParam( numBoxesY, yBoxHalfWidth, xBoxHalfWidth, zBoxHalfWidth); 

    // G4PVParameterised *paramLevelTwoPhys=
    new G4PVParameterised("Level 2 blocks in y",
			  variLog,
			  slabLog, 
			  kYAxis,
			  numBoxesY,
			  pSecondLevelParam);

    G4cout << " Slab dimensions (half-width) are: " << G4endl
	   << " x= " << xSlabHalfWidth/cm << " cm  "
	   << " y= " << ySlabHalfWidth/cm << " cm  "
	   << " z= " << zSlabHalfWidth/cm << " cm  " << G4endl << G4endl; 

    G4cout << " Box dimensions (half-width) are: " << G4endl
	   << " x= " << xBoxHalfWidth/cm << " cm  "
	   << " y= " << yBoxHalfWidth/cm << " cm  "
	   << " z= " << zBoxHalfWidth/cm << " cm  " << G4endl << G4endl; 


    // Other volumes
    G4Box *myMediumBox=new G4Box("Med Box", 25.*cm,25.*cm,25.*cm);
    G4LogicalVolume *medLog=new G4LogicalVolume(myMediumBox,Al,
						 "medBox-LV"); // ,0,0,0);
    new G4PVPlacement(0, G4ThreeVector(-500.*cm,  500.*cm, 0.),
		      "Target-X+Y", medLog, worldPhys, false, 1);
    new G4PVPlacement(0, G4ThreeVector( 500.*cm, -500.*cm, 0.),
		      "Target+X-Y", medLog, worldPhys, false, 1);

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
    G4int copyNo= -1; 

#ifdef ALL_TESTS      
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="World");

    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));

// Check relative search that causes backup one level and then search down:
// Nonrel' finds Target 3, then rel' with point in Target 5 finds Target 5
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(151.*cm,0,-10),0,false);
    assert(located->GetName()=="Level 2 blocks in y");

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(500.*cm,-510.*cm,0));
    assert(located->GetName()=="Target+X-Y");
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),G4ThreeVector(0.,-10.*cm,0.)));
// Check that outside point causes stack to unwind
    assert(!myNav.LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0)));
// Check parameterised volumes
// ---------------------------------------------------------
// Replication 0
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-340.*cm,-95.*cm,-2.5*cm));
    assert(located->GetName()=="Level 2 blocks in y");
    copyNo= located->GetCopyNo();
    // G4cout << " Located ( -340.*cm, -95.*cm, -2.5*cm ) in " 
    //   << located->GetName() << " copy no " << copyNo << G4endl; 
    assert(located->GetCopyNo() == 0 ); 
    // Mother copy/replica number should be 0 

    //  Center of volume should be at ( -340 cm, -90 cm, 0 )
    G4ThreeVector localCoords ( 0.*cm,-5.*cm, -2.5*cm ); 
    // G4cout << " Local coordinates: " << G4endl
    //	     << "  Expected " << localCoords << G4endl
    //       << "  Obtained " << myNav.CurrentLocalCoordinate() << G4endl;
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(), localCoords )); 
            //       G4ThreeVector(0.*cm,-5.*cm,-2.5*cm)));
    assert(located->GetLogicalVolume()->GetMaterial()==brightMaterial); 

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(500.*cm,-510.*cm,0));
    assert(located->GetName()=="Target+X-Y");


// Replication 1
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(-340.*cm,-75.*cm,-2.5*cm));
    assert(located->GetName()=="Level 2 blocks in y");
    // copyNo= located->GetCopyNo();
    assert(located->GetCopyNo() == 1 ); 
    // Mother copy/replica number should be 0
    //
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),
		       G4ThreeVector(0.0*cm,-5.*cm,-2.5*cm)));
    assert(located->GetLogicalVolume()->GetMaterial()==darkMaterial); 

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(500.*cm,-510.*cm,0));
    assert(located->GetName()=="Target+X-Y");

// Replication 2
    /// ....
#endif 

    // Forward part

    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector( 280.*cm, 75.*cm, 25.*cm));
    // Position inside Top Box is     +30 cm,  75 cm, 25 cm
    // Box extent: -100 cm to +100 cm in X,  cut in 10 slabs
    // Slab number(X): #6, local position +0 cm,  75 cm, 25 cm
    // Slab extent: -10 cm to +10 cm in X, -100 cm to 100 cm in Y, Z
    // Slice number(Y, from 10): #8, local position:  0 cm, 5 cm, 25 cm.

    copyNo= located->GetCopyNo();
    G4cout << " Located ( 280.*cm, 75.*cm, 25.*cm ) in  '" 
           << located->GetName() << "' copy no " << copyNo << G4endl; 
    G4cout << " Local coordinates " << myNav.CurrentLocalCoordinate() << G4endl;
    G4cout << G4endl;
    assert(located->GetName()=="Level 2 blocks in y");
    assert(located->GetCopyNo() == 8 ); 
    // Mother copy/replica number should be 0
    //
    assert(ApproxEqual(myNav.CurrentLocalCoordinate(),
		       G4ThreeVector(0.0*cm, 5.*cm, 25.*cm)));
    assert(located->GetLogicalVolume()->GetMaterial()==brightMaterial); 

    /*---------------------------------------------------------------------------
    G4PhysicalTouchable *locPT= dynamic_cast<G4PhysicalTouchable*>(located); 
    if( locPT != 0 ){ 
      G4VPhysicalVolume *parent= locPT->GetParentTouchable()->GetVolume(); 
      G4cout << " **   Parent volume " << locPT << G4endl
	     << "  Expected '" << "Slab Blocks in X" << "'" << G4endl 
	     << "  Obtained '" << parent->GetName() << "' copy no " << parent->GetCopyNo()
	     << G4endl; 
      assert(parent->GetCopyNo() == 6 ); 

      G4VPhysicalVolume *parent2= locPT->GetParentTouchable()->GetVolume(1); 
      if( parent2 != 0){
	G4cout << " **** Parent 2 volume " << locPT << G4endl
	       << "  Expected " << "Top 2-pv" << G4endl
	       << "  Obtained " << parent2->GetName() << " copy no " << parent2->GetCopyNo()
	       << G4endl; 
	assert(parent2->GetName()=="Top 2-pv");
      }
    }
    *****************************************************************************/

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
    G4ThreeVector startXm4( -400.*cm, 0., -10.*cm ); 
    located=myNav.LocateGlobalPointAndSetup( startXm4 );
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( startXm4, xHat, physStep, safety);
    assert(ApproxEqual(Step, 50.*cm));
//    assert(ApproxEqual(safety,5));
    assert(safety>=0);

    G4ThreeVector startXm1( -100.*cm, 0., -10.*cm ); 
    located=myNav.LocateGlobalPointAndSetup(startXm1);
    assert(located->GetName()=="World");
    physStep=kInfinity;
    Step=myNav.ComputeStep( startXm1, mxHat, physStep, safety);
    assert(ApproxEqual(Step,50.*cm));
    assert(ApproxEqual(safety,50.*cm));
    // assert(safety>=0);

    G4ThreeVector newPoint = startXm1 + Step * mxHat;
    myNav.SetGeometricallyLimitedStep();
    located=myNav.LocateGlobalPointAndSetup(newPoint,0,true);
    assert(located->GetName()=="Level 2 blocks in y"); 

    return true;

    // The following tests depend on physical touchables -- obsolete
    // 
    /**********************************************************************
    G4PhysicalTouchable *locPT= dynamic_cast<G4PhysicalTouchable*>(located); 
    G4VPhysicalVolume *parent= locPT->GetParentTouchable()->GetVolume();  

    G4VPhysicalVolume *parent2= locPT->GetParentTouchable()->GetVolume(1); 
    G4cout << " Parent 2 volume " << locPT << G4endl
	   << "  Expected " << "Top 1-pv" << G4endl
    	   << "  Obtained " << parent2->GetName() << " copy no " << parent->GetCopyNo()
	   << G4endl; 

    assert(parent2->GetName()=="Top 1-pv");
    ***********************************************************************/

    G4cerr << " Testing in TestNavigator2() is ending line " << __LINE__ 
	   << " for the time being. " << G4endl; 
    return true;
// ------------------------------------------------------------------
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
    // Must end with geometry open

    return 0;
}
