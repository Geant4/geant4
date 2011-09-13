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
// $Id: testG4PathFinder.cc,v 1.8 2007-02-13 16:15:34 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $ 
//
// 
// Moving with PathFinder in simple boxlike geometry, 

//  TO DO //    both with and without voxels. 

#include <assert.h>
#include "ApproxEqual.hh"

// Global defs
#include "globals.hh"

// To create the geometry 
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

//   For the unit under test 
//                                 #include "G4Navigator.hh"
#include "G4PathFinder.hh"

#include <iomanip>

G4bool  printing= true;    //  false for minimum output, true to debug/check

// Build simple geometry:
// 4 small cubes + 1 slab (all G4Boxes) are positioned inside a larger cuboid
G4VPhysicalVolume* BuildGeometry()
{

    G4Box *myBigBox= new G4Box ("BixBox",250*cm,250*cm,200*cm);

    G4LogicalVolume *worldLog=
      new G4LogicalVolume(myBigBox,0, "World",0,0,0);
    // No material, field, sensitive detector or user limits
    G4PVPlacement *worldPhys=
      new G4PVPlacement(0,G4ThreeVector(0,0,0),
			"World",worldLog, 0,false,0);
    // No displacement, no mother pointer - world!

    G4Box *myBox=new G4Box("cube-10",10*cm,10*cm,10*cm);
    G4LogicalVolume *boxLog=
                new G4LogicalVolume(myBox,0, "Crystal Box",0,0,0);

    new G4PVPlacement(0,G4ThreeVector(15*cm,0,0),  "Target 01",boxLog,
		      worldPhys,false,1);
    new G4PVPlacement(0,G4ThreeVector(25*cm,0,0),  "Target 02",boxLog,
		      worldPhys,false,2);
    new G4PVPlacement(0,G4ThreeVector(-15*cm,0,0), "Target 03",boxLog,
		      worldPhys,false,3);
    new G4PVPlacement(0,G4ThreeVector(-25*cm,0,0), "Target 04",boxLog,
		      worldPhys,false,4);
    new G4PVPlacement(0,G4ThreeVector(0,15*cm,0),  "Target 11",boxLog,
		      worldPhys,false,11);
    new G4PVPlacement(0,G4ThreeVector(0,25*cm,0),  "Target 12",boxLog,
		      worldPhys,false,12);
    new G4PVPlacement(0,G4ThreeVector(0,-15*cm,0), "Target 13",boxLog,
		      worldPhys,false,13);
    new G4PVPlacement(0,G4ThreeVector(0,-25*cm,0), "Target 14",boxLog,
		      worldPhys,false,14);
    new G4PVPlacement(0,G4ThreeVector(0,0,15*cm),  "Target 21",boxLog,
		      worldPhys,false,21);
    new G4PVPlacement(0,G4ThreeVector(0,0,25*cm),  "Target 22",boxLog,
		      worldPhys,false,22);
    new G4PVPlacement(0,G4ThreeVector(0,0,-15*cm), "Target 23",boxLog,
		      worldPhys,false,23);
    new G4PVPlacement(0,G4ThreeVector(0,0,-25*cm), "Target 24",boxLog,
		      worldPhys,false,24);


    new G4PVPlacement(0,G4ThreeVector(-15*cm,15*cm,-10*cm), "Target 101",boxLog,
		      worldPhys,false,1);
    new G4PVPlacement(0,G4ThreeVector(-15*cm,-15*cm,-10*cm), "Target 102",boxLog,
		      worldPhys,false,2);

    G4Box *mySlab= new G4Box("slab",10*cm,25*cm,10*cm);
    G4LogicalVolume *slabLog=new G4LogicalVolume(mySlab,0,
						 "Crystal Slab",0,0,0);
    // G4PVPlacement *offYPhys=
    new G4PVPlacement(0,G4ThreeVector(15*cm,0,-10*cm), "Target 3",slabLog,
		      worldPhys,false,3);

    new G4PVPlacement(0,G4ThreeVector(0,15*cm,10*cm), "Target 4",boxLog,
		      worldPhys,false,4);
    
    new G4PVPlacement(0,G4ThreeVector(0,-15*cm,10*cm), "Target 5",boxLog,
		      worldPhys,false,5);


    return worldPhys;
}

#include "G4TransportationManager.hh"

G4int navIdExpected= 0;  // New convention in G4TransportationManager

//
// Test LocateGlobalPointAndSetup
//
G4PathFinder* setupPathFinder(G4VPhysicalVolume *pTopNode)
{
    MyNavigator myNav;
    G4Navigator* pNav;

    G4int navId;
    G4bool useMyNav= false;
    // G4bool overwriteNav= false;  // Default
    G4bool overwriteNav= true;  // Second option 'new'

    pNav= new G4Navigator();   // (&myNav); 
    // pNav= (&myNav); 

    if( useMyNav ){ 
       pNav= (&myNav); 
       myNav.SetWorldVolume(pTopNode);
    }else{
       pNav->SetWorldVolume(pTopNode);
    }
 
    G4PathFinder* pathFinder= G4PathFinder::GetInstance();
               //===========  --------------------------
    static G4TransportationManager* transportMgr=
      G4TransportationManager::GetTransportationManager(); 
  
    G4Navigator*
    origNav= transportMgr->GetNavigatorForTracking();
    origNav->SetWorldVolume(pTopNode);

    navId= transportMgr->ActivateNavigator( origNav ); 
    G4cout << " navId for original Navigator for tracking is " << navId << G4endl;
    
    if( overwriteNav ){ 
      transportMgr->DeActivateNavigator( origNav );  
      G4cout << " DeActivated Navigator " << origNav << G4endl;  

      // G4bool registered=  transportMgr->RegisterNavigator( pNav ); 
      // assert( registered ); 
      transportMgr->SetNavigatorForTracking(pNav);
      G4cout << " Setting new Navigator for Tracking " << pNav << G4endl;  
      // transportMgr->SetWorldVolume( pTopNode );

      navId= transportMgr->ActivateNavigator( pNav ); 
      G4cout << " navId for new  Navigator for tracking is " << navId << G4endl;
      assert ( navId == navIdExpected ); 
    }
 
    return pathFinder; 
}

G4bool testG4PathFinder1(G4VPhysicalVolume *) // pTopNode)
{
    G4VPhysicalVolume *located;
    static G4TransportationManager* transportMgr=
      G4TransportationManager::GetTransportationManager(); 
    static G4PathFinder* pathFinder= G4PathFinder::GetInstance();

    G4Navigator *pNav= transportMgr->GetNavigatorForTracking();

    G4int navId= navIdExpected;  // As checked/ asserted above !

    G4ThreeVector position( 0., 0., 0.), dirUx(1.,0.,0.); 
    pathFinder->PrepareNewTrack( position, dirUx ); 
  //==========  ---------------
    // Also locates !!

    // Need the volume information 
    // pathFinder->Locate( endPoint, endDirection ); 
     //==========  ------
    located= pathFinder->GetLocatedVolume( navId ); 

    G4double t0=0.0, Ekin=100.00*MeV, restMass= 0.511*MeV, charge= -1, magDipole=0.0, s=0.0; 
    G4ThreeVector Spin(0.,0.,0.); 
    G4FieldTrack  startFT( position, t0, dirUx, Ekin, restMass, charge, Spin, magDipole, s), 
                  endFT('a');     //  Default values for return

    G4int stepNo=0; 
    G4double steplen= 123.00 * mm, safetyRet=0.0; 
    ELimited  limited;
    G4cout << " test:  input FT = " << startFT << G4endl
	   << "                   " 
	   << " len = " << steplen
	   << " navId = " << navId << " stepNo = " << stepNo << G4endl; 
    G4double stepdone= 
    pathFinder->ComputeStep( startFT, steplen, navId, stepNo, safetyRet, limited, endFT, located );
  //==========  -----------


    G4ThreeVector endPoint = endFT.GetPosition(); 
    G4ThreeVector endDirection= endFT.GetMomentumDirection(); 
    pathFinder->Locate( endPoint, endDirection ); 
  //==========  ------
    located= pathFinder->GetLocatedVolume( navId ); 
    if (located && printing ) {
      G4cout << " Located in volume " << located->GetName()
	     << "  id= " << located->GetCopyNo() << G4endl; 
    }

    pathFinder->ComputeStep( startFT, steplen, navId, stepNo, safetyRet, limited, endFT, located );   
    // Should not move, since the 'stepNo' is the same !!
 
    G4double endSafety= pathFinder->ComputeSafety( endFT.GetPosition() );  

    G4cout.precision(5); 
    do{ 
      startFT= endFT; 

      stepdone= 
        pathFinder->ComputeStep( startFT, steplen, navId, ++stepNo, safetyRet, limited, endFT, located );   
      pathFinder->Locate( endFT.GetPosition(), endFT.GetMomentumDirection() ); 
      located= pathFinder->GetLocatedVolume( navId ); 

      if( std::fabs(safetyRet-endSafety) > 1.0e-9 * safetyRet ) { 
	G4cerr << " Problem in safety at new point " 
	       << " Last endpoint returned " << endSafety
	       << " while new computeStep gives " 
	       << G4endl;
      }

      endSafety= pathFinder->ComputeSafety( endFT.GetPosition() ); 

      G4double truestep= (stepdone < steplen ) ? stepdone : steplen; 
      if ( printing ) {
	G4cout << " Step " << std::setw(3) << stepNo << "   " 
	       << " start-safety= " << std::setw(7) << safetyRet << "  " 
	       << " step-length= "  << std::setw(7) << truestep / mm << " mm " 
	       << " to " << endFT.GetPosition() << "  " ; 
	if( located ) {
	    G4cout
	       << " new volume= " << std::setw(10) << located->GetName()
	       << "  copyNo= " << located->GetCopyNo(); 
	}
	G4cout << " end-safety= " << endSafety;
	if( located )  G4cout << G4endl;
      }
    } while ( located );

    G4cout << "  Last Step it exits the World. Located = " << located << G4endl; 
    G4cout << G4endl;
    // G4cout << " Step " << stepNo << " is the last one, it exits the World." << G4endl; 
    return true; 

    ////  OLD Checks ......
    assert(!pNav->LocateGlobalPointAndSetup(G4ThreeVector(kInfinity,0,0),0,false));
    located=pNav->LocateGlobalPointAndSetup(G4ThreeVector(0,0,0),0,false);
    assert(located->GetName()=="World");


    return true;
 
    // transportMgr->DeActivateNavigator( pNav ); 
    // ----- // transportMgr->DeRegisterNavigator( pNav );  // Cannot de-register Mass Navigator
    //  exit(1); 
    
    //--------------------   END                   --------------------------------
    // G4VPhysicalVolume *located;
    MyNavigator& myNav= dynamic_cast<MyNavigator&>(*pNav);

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
    G4PathFinder* pPathFinder;

    pPathFinder= setupPathFinder(myTopNode);

    testG4PathFinder1(myTopNode);
    // testG4Navigator2(myTopNode);
// Repeat tests but with full voxels
    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
    //  exit(1); 

    // --- Future 2nd test 

    G4GeometryManager::GetInstance()->CloseGeometry(true);
    testG4PathFinder1(myTopNode);
    // testG4Navigator2(myTopNode);

    G4GeometryManager::GetInstance()->OpenGeometry();
    return 0;
}
