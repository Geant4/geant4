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
// $Id: testG4PVDivision.cc,v 1.1 2003-06-16 15:11:43 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// test for G4PVDivision classes
//
// 26.05.03 - P.Arce Initial version
// ********************************************************************

#include <assert.h>
#include <fstream>
#include <vector>

#include "G4ios.hh"
#include "globals.hh"

#include "G4Navigator.hh"

#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVParameterised.hh"
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Para.hh"
#include "G4Polycone.hh"

#include "G4GeometryManager.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"

#include "Randomize.hh"
#include "G4PVReplica.hh"

#include "G4PVDivision.hh"

G4String theSolidTypeStr;
G4String thePVTypeStr;
G4double theWorldDim = 1*m;
G4int numberOfPoints = 1000;
G4int theNReplicas;
G4double theWidth;
G4double theOffset;
G4VSolid* theParentSolid;
std::vector<G4LogicalVolume*> theParentLogs;
std::vector<G4VPhysicalVolume*> theParentPhyss;
std::vector<G4VSolid*> theChildSolids;
std::vector<G4LogicalVolume*> theChildLogs;
std::vector<EAxis> theAxis;
std::vector<G4double> theWidths;
std::vector<G4double> theExtraPars;
G4int theDivType;

enum SolidType{g4box, g4trd, g4tube, g4tubs, g4cone, g4cons, g4polycone, g4polyhedra };
enum PVType{pvDivision, pvReplica, pvPlacement };

//--------------------------------------------------------------------------
void initialize();
void calculateParentSolid( SolidType solType );
void calculateChildSolids( SolidType solType, G4VSolid* pSolid );
void calculateAxis( SolidType solType );
void buildOutputName( SolidType& soltype, PVType pvtype );
void starttest( const std::vector<G4String>& vsarg );
G4int checkNumberOfArguments( const G4String& st, G4int narg );
PVType getPVType( const G4String& pvt );
SolidType getSolidType( const G4String& st );
void generateRandomPoints();
void generateScanPointsForBox();
void generateScanPointsForTube();
void generateScanPointsForTrd();
void generateScanPointsForPolycone();
G4VPhysicalVolume* BuildGeometry( SolidType solType, PVType pvType );
G4bool testG4Navigator1(G4VPhysicalVolume *pTopNode);
G4bool testG4Navigator2(G4VPhysicalVolume *pTopNode);
// World geometry is a box 1 X 1 X 3. 
// User select the type of solid
// Inside it three solids of the chosen type are placed along Z
// Each of these solids is divided along a different axis
//
// Then a set of points is generated and it is tested in which division copy they are and which is the DistanceToOut in the three directions (X,Y,Z)

//--------------------------------------------------------------------------
int main( G4int argc, char** argv ) 
{
  // first argument is type of divisioning (repli/divis), second is type of solid
  std::vector<G4String> vsarg;
  for( G4int jj = 0; jj < argc; jj ++ ) {
    vsarg.push_back( G4String(argv[jj] ) );
  } 

  starttest( vsarg );
}

//--------------------------------------------------------------------------
void starttest( const std::vector<G4String>& vsarg ) 
{
  G4int narg = vsarg.size();
  if( narg == 1 ) {
    thePVTypeStr = "divis";
    theSolidTypeStr = "box";
  } else if( narg == 2 ){
    // wrong number
    checkNumberOfArguments( " ", narg );
  } else { 
    G4int divTypeSet = checkNumberOfArguments( vsarg[2], narg );
    if( divTypeSet ) {
      theDivType = atoi( vsarg[divTypeSet] );
    } else {
      theDivType = 0;
    }
    thePVTypeStr = G4String(vsarg[1]);
    theSolidTypeStr = G4String(vsarg[2]);
  }
  if( narg > 3 ){
    for( G4int ii = 3; ii < narg; ii++ ){
      theExtraPars.push_back( atof( vsarg[ii] ) );
    }
  }
  PVType pvtype = getPVType( thePVTypeStr );
  SolidType soltype = getSolidType( theSolidTypeStr );

  //  PVType pvtype = pvDivision;
  initialize();
  if( soltype == g4tube || soltype == g4tubs ) {
    // generateRandomPoints();
    generateScanPointsForTube();
  }else if( soltype == g4cone || soltype == g4cons ) {
    // generateRandomPoints();
    generateScanPointsForTube();
  } else if( soltype == g4box ) {
    //     generateRandomPoints();
    generateScanPointsForBox();
  } else if( soltype == g4trd ) {
    //     generateRandomPoints();
    generateScanPointsForTrd();
  } else if( soltype == g4polycone ) {
    //     generateRandomPoints();
    generateScanPointsForPolycone();
  } else {
    generateRandomPoints();
  }

  G4VPhysicalVolume *myTopNode;
  myTopNode=BuildGeometry( soltype, pvtype );	// Build the geometry
  G4GeometryManager::GetInstance()->CloseGeometry(false);
  testG4Navigator1(myTopNode);
  testG4Navigator2(myTopNode);
  // Repeat tests but with full voxels
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4GeometryManager::GetInstance()->CloseGeometry(true);
  testG4Navigator1(myTopNode);
  testG4Navigator2(myTopNode);
  
  G4cout << " theParentSolid " << G4endl;
  G4cout << *theParentSolid << G4endl;
  for( size_t ii = 0; ii < theChildSolids.size(); ii++) {
    G4cout << " theChildSolid after tracking " << "  "<<  ii << G4endl;
    G4cout <<  *theChildSolids[ii] << G4endl;
  }
 
  G4GeometryManager::GetInstance()->OpenGeometry();
}

//--------------------------------------------------------------------------
void generateRandomPoints()
{
  std::ofstream fout("points.lis");
  G4double posX, posY, posZ;

  for( G4int ii = 0; ii < numberOfPoints; ii++ ) {
    posX = RandFlat::shoot(-theWorldDim, theWorldDim);
    posY = RandFlat::shoot(-theWorldDim, theWorldDim);
    posZ = RandFlat::shoot(-3*theWorldDim, 3*theWorldDim);
    fout << posX << " " << posY << " " << posZ << G4endl;
  }
}

//--------------------------------------------------------------------------
void generateScanPointsForBox()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerReplica = 2;
  G4int numberOfPoints = theNReplicas*nPointsPerReplica;
  // For division along X
  G4ThreeVector centre(0.,0.,-2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any Z, any Y
    G4ThreeVector pR( 0., theWorldDim/100., theWorldDim/100. );
    G4double X = -theWorldDim + (ii+0.001) * 2*theWorldDim/numberOfPoints;
    pR.setX( X );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Y
  centre = G4ThreeVector(0.,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any X, any Z
    G4ThreeVector pR( theWorldDim/100., 0., theWorldDim/100. );
    G4double Y = -theWorldDim + (ii+0.001) * 2*theWorldDim/numberOfPoints;
    pR.setY( Y );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Z
  centre = G4ThreeVector(0.,0.,2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any X, any Y
    G4ThreeVector pR( theWorldDim/100., 0., theWorldDim/100. );
    G4double Z = -theWorldDim + (ii+0.001) * 2*theWorldDim/numberOfPoints;
    pR.setZ( Z );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}

//--------------------------------------------------------------------------
void generateScanPointsForTrd()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerReplica = 2;
  G4int numberOfPoints = theNReplicas*nPointsPerReplica;
  // For division along X
  G4ThreeVector centre(0.,0.,-2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any Z, any Y
    G4ThreeVector pR( 0., theWorldDim/100., theWorldDim/100. );
    G4double X = -theWorldDim*0.75 + (ii+0.001) * 1.5*theWorldDim/numberOfPoints;
    pR.setX( X );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Y
  centre = G4ThreeVector(0.,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any X, any Z
    G4ThreeVector pR( theWorldDim/100., 0., theWorldDim/100. );
    G4double Y = -theWorldDim*0.75 + (ii+0.001) * 1.5*theWorldDim/numberOfPoints;
    pR.setY( Y );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Z
  centre = G4ThreeVector(0.,0.,2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any X, any Y
    G4ThreeVector pR( theWorldDim/100., 0., theWorldDim/100. );
    G4double Z = -theWorldDim + (ii+0.001) * 2*theWorldDim/numberOfPoints;
    pR.setZ( Z );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}

//--------------------------------------------------------------------------
void generateScanPointsForTube()
{
  std::ofstream fout("points.lis");
  G4int ii;

  G4int nPointsPerReplica = 2;
  G4int numberOfPoints = theNReplicas*nPointsPerReplica;
  G4cout << " numberOfPoints " << numberOfPoints << G4endl; 
  // For division along R
  G4ThreeVector centre(0.,0.,-2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    // any Z, phi close initial phi
    G4double phi = 1.*deg;
    //t?    if( theExtraPars.size() > 0 ) phi = (theExtraPars[0]+1.)*deg;
    G4ThreeVector pR( cos(phi), sin(phi), theWorldDim/100. );
    G4double rho = (ii+0.001) * theWorldDim/numberOfPoints;
    pR.setRho( rho );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along phi
  centre = G4ThreeVector(0.,0.,0.);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    G4double phi = (ii+0.001) * 360*deg/numberOfPoints;
    // any rho (=1), any Z
    G4ThreeVector pR( cos(phi), sin(phi), theWorldDim/100. );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }

  // For division along Z
  centre = G4ThreeVector(0.,0.,2*theWorldDim);
  for( ii = 0; ii < numberOfPoints; ii++ ) {
    //any rho (=1), phi close initial phi
    G4double phi = 1.*deg;
    //t?  if( theExtraPars.size() > 0 ) phi = (theExtraPars[0]+1.)*deg;
    G4ThreeVector pR( cos(phi), sin(phi),0.);
    G4double Z = -theWorldDim + (ii+0.001) * 2*theWorldDim / numberOfPoints;
    pR.setZ( Z );
    pR += centre;
    fout << pR.x() << " " << pR.y() << " " << pR.z() << G4endl;
  }
}

//--------------------------------------------------------------------------
void generateScanPointsForPolycone()
{
  generateScanPointsForTube();
}

//--------------------------------------------------------------------------
// Build simple geometry:
//  world is 
G4VPhysicalVolume* BuildGeometry( SolidType solType, PVType pvType )
{
  G4int ii;
  // The world volume
  //
  G4Box *worldSolid= new G4Box ("Big Cube", theWorldDim, theWorldDim, 3*theWorldDim);
  
  G4LogicalVolume *worldLog=new G4LogicalVolume(worldSolid,0,
						"World",0,0,0);
				// Logical with no material,field,
                                // sensitive detector or user limits
  
  G4PVPlacement *worldPhys=new G4PVPlacement(0,G4ThreeVector(0,0,0),
					     "World",worldLog,
					     0,false,0);
				// Note: no parent pointer set
  
    // build three boxes 
    // A set of boxes
  calculateParentSolid( solType );

  //parent logical volumes do not depend on SolidType, PVType
  theParentLogs.push_back( new G4LogicalVolume(theParentSolid,0, "Parent1",0,0,0) );
  theParentLogs.push_back( new G4LogicalVolume(theParentSolid,0, "Parent2",0,0,0) );
  theParentLogs.push_back( new G4LogicalVolume(theParentSolid,0, "Parent3",0,0,0) );
  
  //parent physical volumes positions do not depend on SolidType, PVType
  for( ii = 0; ii < 3; ii++ ) {
    theParentPhyss.push_back( new G4PVPlacement( 0, G4ThreeVector(0.,0.,(ii-1)*2*theWorldDim), theParentLogs[ii] , "parent", worldLog, 0, 0 ) );
  }

  // children
  calculateChildSolids( solType, theParentSolid );
  calculateAxis( solType );
    for( ii = 0; ii < 3; ii++ ) {
      theChildLogs.push_back( new G4LogicalVolume(theChildSolids[ii], 0, "child",0,0,0) );
    }

  if( pvType == pvDivision ) {
    for( ii = 0; ii < 3; ii++ ) {
      if( theDivType == 0 ) {
	new G4PVDivision("child", theChildLogs[ii], theParentLogs[ii], 
			 theAxis[ii],
			 theNReplicas, 
			 theWidths[ii], 
			 theOffset );
	G4cout << "division NDIVandWIDTH " << theNReplicas << " " << theWidths[ii]<< " " << theOffset << G4endl;
      } else if( theDivType == 1 ) {
	new G4PVDivision("child", theChildLogs[ii], theParentLogs[ii], 
			 theAxis[ii],
			 theWidths[ii], 
			 theOffset );
	G4cout << "division WIDTH " << theWidths[ii]<< " " << theOffset << G4endl;
      }else if( theDivType == 2 ) {
	new G4PVDivision("child", theChildLogs[ii], theParentLogs[ii], 
			 theAxis[ii],
			 theNReplicas, 
			 theOffset );
	G4cout << "division NDIV " << theNReplicas << " " << theOffset << G4endl;
      }
    }
  } else if( pvType == pvReplica ) {
 
    for( ii = 0; ii < 3; ii++ ) {
      new G4PVReplica("child", theChildLogs[ii], theParentLogs[ii], 
		    theAxis[ii],
		    theNReplicas, 
		    theWidths[ii], 
		    theOffset );
      G4cout << "replica " <<  theNReplicas << " " << theWidths[ii]<< " " << theOffset << G4endl;
    }
  }
  return worldPhys;
}

//--------------------------------------------------------------------------
//
// Test LocateGlobalPointAndSetup
//
G4bool testG4Navigator1(G4VPhysicalVolume *pTopNode)
{
    G4Navigator myNav;
    G4VPhysicalVolume *located;
    myNav.SetWorldVolume(pTopNode);
    
    G4String foutname = "points." + thePVTypeStr + "." + theSolidTypeStr + ".out";
    std::ofstream fout(foutname);
    std::ifstream fin("points.lis");
    G4double posX, posY, posZ;
    for( G4int ii = 0; ii < numberOfPoints; ii++ ) {
      fin >> posX >> posY >> posZ;
      if( fin.eof() ) break;
      located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ),0,false);
      G4ThreeVector pos(posX,posY,posZ);
      fout << ii << pos << " " << located->GetName() << " " << located->GetCopyNo();
      if( theSolidTypeStr == "tubs" || theSolidTypeStr == "tube" ){
	//	fout << " " << pos.phi()/deg << G4endl;
	fout << G4endl;
      } else {
	fout << G4endl;
      }
    }

    return true;
}

//--------------------------------------------------------------------------
//
// Test Stepping
//
G4bool testG4Navigator2(G4VPhysicalVolume *pTopNode)
{
  G4Navigator myNav;
  G4VPhysicalVolume *located;
  G4double Step,physStep,safety;
  G4ThreeVector* Hat = new G4ThreeVector[3];
  Hat[0] = G4ThreeVector(1,0,0);
  Hat[1] = G4ThreeVector(0,1,0);
  Hat[2] = G4ThreeVector(0,0,1);
  
  myNav.SetWorldVolume(pTopNode);
  
  G4String foutname = "step." + thePVTypeStr + "." + theSolidTypeStr + ".out";
  std::ofstream fout(foutname);
  std::ifstream fin("points.lis");
  G4double posX, posY, posZ;
  for( G4int ii = 0; ii < numberOfPoints; ii++ ) {
    if( fin.eof() ) break;
    fin >> posX >> posY >> posZ;
    located=myNav.LocateGlobalPointAndSetup(G4ThreeVector(posX,posY,posZ));
    physStep=kInfinity;
    for( G4int jj = 0; jj < 3; jj++ ) {
      Step=myNav.ComputeStep(G4ThreeVector(posX,posY,posZ),Hat[jj],physStep,safety);
      fout << Step << " ";
    }
    fout << G4endl; 
  }

  return true;
}

//--------------------------------------------------------------------------
void initialize()
{
 theNReplicas = 10;
 theOffset = 0.;
}

//--------------------------------------------------------------------------
void calculateParentSolid( SolidType solType )
{
  if( solType == g4box ) {
    theParentSolid = new G4Box("parent", theWorldDim, theWorldDim, theWorldDim);
  }else if( solType == g4trd ) {
    theParentSolid = new G4Trd("parent", theWorldDim/2., theWorldDim, theWorldDim/2., theWorldDim, theWorldDim);
  }else if( solType == g4tube ) {
    theParentSolid = new G4Tubs("parent", 0., theWorldDim, theWorldDim, 0., 360.*deg);
  }else if( solType == g4tubs ) {
    theParentSolid = new G4Tubs("parent", 0., theWorldDim, theWorldDim, theExtraPars[0]*deg, theExtraPars[1]*deg);
  }else if( solType == g4cone ) {
    theParentSolid = new G4Cons("parent", 1.E-6*mm, theWorldDim, 1.E-6*mm, theWorldDim/2., theWorldDim, 0., 360.*deg);
  }else if( solType == g4cons ) {
    theParentSolid = new G4Cons("parent", 1.E-6*mm, theWorldDim/2., 1.E-6*mm, theWorldDim, theWorldDim, theExtraPars[0]*deg, theExtraPars[1]*deg);
  }else if( solType == g4polycone ) {

    G4double* zPlane = new G4double(4);
    zPlane[0] = -theWorldDim; zPlane[1] = -theWorldDim/2.; zPlane[2] = theWorldDim/2; zPlane[3] = theWorldDim;
    G4double* rInner = new G4double(4);
    rInner[0] = theWorldDim/10.; rInner[1] = 0.; rInner[2] = 0.; rInner[3] = theWorldDim/10.;
    G4double* rOuter = new G4double(4);
    rOuter[0] = theWorldDim; rOuter[1] = theWorldDim*0.9; rOuter[2] = theWorldDim*0.9; rOuter[3] = theWorldDim;

    theParentSolid = new G4Polycone("parent", 0., 360.*deg, 4, zPlane, rInner, rOuter);
  }
}

//--------------------------------------------------------------------------
void calculateChildSolids( SolidType solType, G4VSolid* pSolid )
{
  if( solType == g4box ) {
    theWidths.push_back( 2*theWorldDim / theNReplicas );
    theWidths.push_back( 2*theWorldDim / theNReplicas );
    theWidths.push_back( 2*theWorldDim / theNReplicas );

    theChildSolids.push_back( new G4Box("child", theWorldDim/theNReplicas, theWorldDim, theWorldDim) );
    theChildSolids.push_back( new G4Box("child", theWorldDim, theWorldDim/theNReplicas, theWorldDim) );
    theChildSolids.push_back( new G4Box("child", theWorldDim, theWorldDim, theWorldDim/theNReplicas) );

  }else if( solType == g4trd ) {
    theWidths.push_back( 1.5*theWorldDim / theNReplicas );
    theWidths.push_back( 1.5*theWorldDim / theNReplicas );
    theWidths.push_back( 2*theWorldDim / theNReplicas );

    theChildSolids.push_back( new G4Trap("child", theWorldDim/2./theNReplicas, theWorldDim/theNReplicas, theWorldDim, theWorldDim, theWorldDim) ); // solid dimensions will be recalculated
    theChildSolids.push_back( new G4Trap("child", theWorldDim, theWorldDim, theWorldDim/2./theNReplicas, theWorldDim/theNReplicas, theWorldDim) ); // solid dimensions will be recalculated
    theChildSolids.push_back( new G4Trd("child", theWorldDim, theWorldDim, theWorldDim, theWorldDim, theWorldDim/theNReplicas) );

  }else if( solType == g4tube || solType == g4tubs ) {
    G4Tubs* pTubs = reinterpret_cast<G4Tubs*>( pSolid );

    theWidths.push_back( (pTubs->GetInnerRadius()+pTubs->GetOuterRadius()) / theNReplicas );
    theWidths.push_back( pTubs->GetDeltaPhiAngle() / theNReplicas );
    theWidths.push_back( 2*pTubs->GetZHalfLength() / theNReplicas );

    theChildSolids.push_back( new G4Tubs("child", pTubs->GetInnerRadius(), pTubs->GetInnerRadius()+theWidths[0], pTubs->GetZHalfLength(), pTubs->GetStartPhiAngle(), pTubs->GetDeltaPhiAngle() ));
    theChildSolids.push_back( new G4Tubs("child", pTubs->GetInnerRadius(), pTubs->GetOuterRadius(), pTubs->GetZHalfLength(), pTubs->GetStartPhiAngle(), theWidths[1] ));
    theChildSolids.push_back( new G4Tubs("child", pTubs->GetInnerRadius(), pTubs->GetOuterRadius(), theWidths[2]/2., pTubs->GetStartPhiAngle(), pTubs->GetDeltaPhiAngle() ));
  }else if( solType == g4cone || solType == g4cons ) {
    G4Cons* pCons = reinterpret_cast<G4Cons*>( pSolid );

    theWidths.push_back( (pCons->GetInnerRadiusMinusZ()+pCons->GetOuterRadiusMinusZ()) / theNReplicas );
    theWidths.push_back( pCons->GetDeltaPhiAngle() / theNReplicas );
    theWidths.push_back( 2*pCons->GetZHalfLength() / theNReplicas );

    G4double fwidthPlus = 0.;
    if( pCons->GetInnerRadiusMinusZ() != 0. ) fwidthPlus = theWidths[0] * pCons->GetInnerRadiusPlusZ() / pCons->GetInnerRadiusMinusZ();
    theChildSolids.push_back( new G4Cons("child", pCons->GetInnerRadiusMinusZ(), pCons->GetInnerRadiusMinusZ()+theWidths[0], pCons->GetInnerRadiusPlusZ(), pCons->GetInnerRadiusPlusZ()+fwidthPlus, pCons->GetZHalfLength(), pCons->GetStartPhiAngle(), pCons->GetDeltaPhiAngle() ));
    theChildSolids.push_back( new G4Cons("child", pCons->GetInnerRadiusMinusZ(), pCons->GetOuterRadiusMinusZ(), pCons->GetInnerRadiusPlusZ(), pCons->GetOuterRadiusPlusZ(), pCons->GetZHalfLength(), pCons->GetStartPhiAngle(), theWidths[1] ));
    theChildSolids.push_back( new G4Cons("child", pCons->GetInnerRadiusMinusZ(), pCons->GetOuterRadiusMinusZ(), pCons->GetInnerRadiusPlusZ(), pCons->GetOuterRadiusPlusZ(), theWidths[2]/2., pCons->GetStartPhiAngle(), pCons->GetDeltaPhiAngle() ));
 
  }else if( solType == g4polycone ) {
    G4Polycone* pPCone = reinterpret_cast<G4Polycone*>( pSolid );
    G4PolyconeHistorical* origparam = pPCone->GetOriginalParameters() ;
    theWidths.push_back( (origparam->Rmax[0] - origparam->Rmin[0] ) / theNReplicas );
    theWidths.push_back( (origparam->Rmax[0] - origparam->Rmin[0] ) / theNReplicas );
    theWidths.push_back( (origparam->Rmax[0] - origparam->Rmin[0] ) / theNReplicas );
    //    theWidths.push_back( pTubs->GetDeltaPhiAngle() / theNReplicas );
    //  theWidths.push_back( 2*pTubs->GetZHalfLength() / theNReplicas );

    G4double* zPlane = new G4double(4);
    zPlane[0] = -theWorldDim; zPlane[1] = -theWorldDim/2.; zPlane[2] = theWorldDim/2; zPlane[3] = theWorldDim;
    G4double* rInner = new G4double(4);
    rInner[0] = theWorldDim/10.; rInner[1] = 0.; rInner[2] = 0.; rInner[3] = theWorldDim/10.;
    G4double* rOuter = new G4double(4);
    rOuter[0] = theWorldDim; rOuter[1] = theWorldDim*0.9; rOuter[2] = theWorldDim*0.9; rOuter[3] = theWorldDim;

    theChildSolids.push_back( new G4Polycone("child", 0., 360.*deg, 4, zPlane, rInner, rOuter) );
    theChildSolids.push_back( new G4Polycone("child", 0., 360.*deg, 4, zPlane, rInner, rOuter) );
    theChildSolids.push_back( new G4Polycone("child", 0., 360.*deg, 4, zPlane, rInner, rOuter) );

  }
  for( size_t ii = 0; ii < theChildSolids.size(); ii++) {
    G4cout << theChildSolids[0] << "child solid built  0 " << *theChildSolids[0] << G4endl;
    G4cout << " theChildSolid built " <<  ii << G4endl;
    G4cout <<  *theChildSolids[ii] << G4endl;
  }
}

//--------------------------------------------------------------------------
void calculateAxis( SolidType solType )
{
  if( solType == g4box ) {
    theAxis.push_back( kXAxis );
    theAxis.push_back( kYAxis );
    theAxis.push_back( kZAxis );
  }else if( solType == g4trd ) {
    theAxis.push_back( kXAxis );
    theAxis.push_back( kYAxis );
    theAxis.push_back( kZAxis );
  }else if( solType == g4tube || solType == g4tubs ) {
    theAxis.push_back( kRho );
    theAxis.push_back( kPhi );
    theAxis.push_back( kZAxis );
  }else if( solType == g4cone || solType == g4cons || solType == g4polycone ) {
    theAxis.push_back( kRho );
    theAxis.push_back( kPhi );
    theAxis.push_back( kZAxis );
  }
}

//--------------------------------------------------------------------------
SolidType getSolidType( const G4String& st )
{
  SolidType stype = g4box;

  if( st == "box" ) {
    stype = g4box;
  }else if( st == "trd" ){
    stype = g4trd;
  }else if( st == "tube" ){
    stype = g4tube;
  }else if( st == "tubs" ){
    stype = g4tubs;
  }else if( st == "cone" ){
    stype = g4cone;
  }else if( st == "cons" ){
    stype = g4cons;
  }else if( st == "polycone" ){
    stype = g4polycone;
  }else if( st == "polyhedra" ){
    stype = g4polyhedra;
  } else {
    G4Exception(" Input solid type can only be 'box', 'trd', 'tube', 'tubs', 'cons', 'polycone', 'polyhedra' " );
  }
  return stype;
}

//--------------------------------------------------------------------------
G4int checkNumberOfArguments( const G4String& st, G4int narg )
{
  G4bool nok = false;
  if( narg == 1 ) return 0;
  G4int nArgExcess;
  if( st == "box" || st == "trd" || st == "tube" || st == "cone" ){
    nArgExcess = narg - 3;
  }else if( st == "tubs" || st == "cons" ){
    nArgExcess = narg - 5;
  }else if( st == "polycone" ){
    nArgExcess = narg - 3;
  }else if( st == "polyhedra" ){
    nArgExcess = narg - 3;
  } else {
    nArgExcess = 10;
  }

  G4int divTypeSet = 0;
  if( nArgExcess == 0 ) {
    nok = true;
  } else if( nArgExcess == 1 ) {
    divTypeSet = narg - 1; //last argument sets the division type 
    nok = true;
  }

  if( !nok ){
    G4cerr << " Number Of arguments " << narg << " for solid type " << st << G4endl;
    G4Exception(" Number of arguments incorrect for input solid type, please check method checkNumberOfArguments " );
  }

  return divTypeSet;
}

//--------------------------------------------------------------------------
PVType getPVType( const G4String& pvt )
{
  G4cout << pvt << G4endl;
  PVType vtype = pvPlacement;

  if( pvt == "divis"){
    vtype = pvDivision; 
  }else if( pvt == "repli") {
    vtype = pvReplica; 
  }else if( pvt == "place") {
    vtype = pvPlacement; 
  } else {
    G4Exception(" Input PV type can only be 'divis', 'repli' or 'place' " );
  }
  return vtype;
}
