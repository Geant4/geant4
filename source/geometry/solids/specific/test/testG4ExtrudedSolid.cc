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
// $Id$
//
// testG4ExtrudedSolid
//
// Test file for class G4ExtrudedSolid.
// In the functions createSolidN(...), there are defined several 
// test cases of extruded solid. These functions also fill in 
// the vectors with explicitely defined points inside, on the surface 
// and outside the solid.
// All the test results for the definesd solids and
// points can be printed via PrintResults() function.
// The tests are then defined in testXYZ() functions
// using assert() on the comparison with the expected
// result value.
//
// The functions DistanceToIn, DistanceToOut on surface
// point do not give always expected values, that's why
// they are not yet included in the tests with assert.
// To be added tests for SurfaceNormal(p) function.
// Ensure asserts are compiled in.
//
// Author:
//   Ivana Hrivnacova, IPN Orsay
// 

#undef NDEBUG
#include <assert.h>
#include <cmath>
#include <vector>
#include <iomanip>

#include "globals.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4ExtrudedSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4Box.hh"
#include "G4UnionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4Timer.hh"

G4double kCarTolerance = G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

G4ThreeVector dirx(1,0,0);
G4ThreeVector diry(0,1,0);
G4ThreeVector dirz(0,0,1);

//_____________________________________________________________________________
G4VSolid* createSolid0(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// Create extruded solid with triangular polygon
// and fill vectors with points

  inside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,   5.0*cm));
  inside_points.push_back(G4ThreeVector(-10.0*cm, -5.0*cm, -15.0*cm));  
  
  surface_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  30.0*cm));
  surface_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,  30.0*cm));
  surface_points.push_back(G4ThreeVector( -5.0*cm, -5.0*cm, -30.0*cm));
  surface_points.push_back(G4ThreeVector(-15.0*cm,  0.0*cm,  10.0*cm));
  surface_points.push_back(G4ThreeVector(+15.0*cm,  0.0*cm, -10.0*cm));

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm, -5.0*cm, -40.0*cm));
  outside_points.push_back(G4ThreeVector(-20.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector(+20.0*cm,  0.0*cm, -10.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm,-35.0*cm, -20.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm, 40.0*cm,  10.0*cm));

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(  0.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  
  return new G4ExtrudedSolid(
               "extrudedSolid1", polygon, 30.*cm, 
               G4TwoVector(), 1.0, G4TwoVector(), 1.0);
}                             

//_____________________________________________________________________________
G4VSolid* createSolid1(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// Create box defined as extruded solid
// and fill vectors with points
// The same solid can be defined as G4TessellatedSolid or G4Box,
// when uncommenting the appropriate lines below.


  inside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,   5.0*cm));
  inside_points.push_back(G4ThreeVector(-10.0*cm, -5.0*cm, -15.0*cm));  
  
  surface_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  30.0*cm));
  surface_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,  30.0*cm));
  surface_points.push_back(G4ThreeVector( 30.0*cm,  0.0*cm,   0.0*cm));
  surface_points.push_back(G4ThreeVector( 30.0*cm, -5.0*cm,  -5.0*cm));
  surface_points.push_back(G4ThreeVector(  0.0*cm,-30.0*cm,   0.0*cm));
  surface_points.push_back(G4ThreeVector(  5.0*cm,-30.0*cm,  -5.0*cm));

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm, -5.0*cm, -40.0*cm));
  outside_points.push_back(G4ThreeVector(-35.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector(+35.0*cm,  0.0*cm, -10.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm,-40.0*cm, -20.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm, 40.0*cm,  10.0*cm));

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  
  G4ExtrudedSolid* extruded
    = new G4ExtrudedSolid(
               "extrudedSolid1", polygon, 30.*cm, 
               G4TwoVector(), 1.0, G4TwoVector(), 1.0);

  //G4TessellatedSolid* tessellated
  //  = new G4TessellatedSolid(*extruded);
  //return tessellated;

  //G4Box* box 
  //  = new G4Box("box", 30.*cm, 30*cm, 30*cm);
  // return box;
  
  return extruded;
}                             

//_____________________________________________________________________________
G4VSolid* createSolid2(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// Create extruded solid with concave polygon with 2 z- sections
// and fill vectors with points

  inside_points.push_back(G4ThreeVector( 10.0*cm, 25.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -20.0*cm));
  inside_points.push_back(G4ThreeVector( 15.0*cm,-15.0*cm,  20.0*cm));  
  
  surface_points.push_back(G4ThreeVector( 20.0*cm, 30.0*cm, -25.0*cm));
  surface_points.push_back(G4ThreeVector(  0.0*cm, 10.0*cm,  25.0*cm));
  surface_points.push_back(G4ThreeVector(-40.0*cm,  0.0*cm,   0.0*cm));
  surface_points.push_back(G4ThreeVector( 20.0*cm, -5.0*cm,   0.0*cm));
  surface_points.push_back(G4ThreeVector( 20.0*cm,  5.0*cm,   0.0*cm));

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  30.0*cm));
  outside_points.push_back(G4ThreeVector( 10.0*cm,  5.0*cm, -40.0*cm));
  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,   0.0*cm));
  outside_points.push_back(G4ThreeVector(-40.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( 40.0*cm,  0.0*cm, -10.0*cm));

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));
  
  return new G4ExtrudedSolid(
               "extrudedSolid3", polygon, 25.*cm, 
                G4TwoVector(-20.*cm, 10.*cm), 1.5, G4TwoVector(), 0.5);
}                             


//_____________________________________________________________________________
G4VSolid* createSolid3(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
  // Extruded solid with the same polygon as solid3 bit with 4 z-sections

  inside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -35.0*cm));
  inside_points.push_back(G4ThreeVector( 10.0*cm, 10.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector(-15.0*cm,  0.0*cm,  12.0*cm));
  inside_points.push_back(G4ThreeVector( 35.0*cm, -5.0*cm,  30.0*cm));  
   
  surface_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -40.0*cm));
  surface_points.push_back(G4ThreeVector( 15.0*cm,  0.0*cm,  10.0*cm));
  surface_points.push_back(G4ThreeVector( -5.0*cm, 10.5*cm,  15.0*cm));
  surface_points.push_back(G4ThreeVector( 45.0*cm, 33.5*cm,  40.0*cm));

  outside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -50.0*cm));
  outside_points.push_back(G4ThreeVector( 25.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm,  5.0*cm,  15.0*cm));
  outside_points.push_back(G4ThreeVector( 45.0*cm, 40.0*cm,  45.0*cm));

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));
  
  std::vector<G4ExtrudedSolid::ZSection> zsections;
  zsections.push_back(G4ExtrudedSolid::ZSection(-40.*cm, G4TwoVector(-20.*cm, 10.*cm), 1.5));
  zsections.push_back(G4ExtrudedSolid::ZSection( 10.*cm, G4TwoVector(  0.*cm,  0.*cm), 0.5));
  zsections.push_back(G4ExtrudedSolid::ZSection( 15.*cm, G4TwoVector(  0.*cm,  0.*cm), 0.7));
  zsections.push_back(G4ExtrudedSolid::ZSection( 40.*cm, G4TwoVector( 20.*cm, 20.*cm), 0.9));

  G4ExtrudedSolid* extruded
   = new G4ExtrudedSolid("extrudedSolid4", polygon, zsections);

  return extruded;

//  G4TessellatedSolid* tessellated
//    = new G4TessellatedSolid(*extruded);
//   return tessellated;
}  
  

//_____________________________________________________________________________
G4VSolid* createSolid4(std::vector<G4ThreeVector>& /*inside_points*/,
                       std::vector<G4ThreeVector>& /*surface_points*/,
                       std::vector<G4ThreeVector>& /*outside_points*/)
{
  // Extruded solid with 4 z-sections, with 2 sections with the same z position
  // defined via union solid

  std::vector<G4TwoVector> polygon;
  polygon.push_back(G4TwoVector(-30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector(-30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm,  30.*cm));
  polygon.push_back(G4TwoVector( 30.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm, -30.*cm));
  polygon.push_back(G4TwoVector( 15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm,  15.*cm));
  polygon.push_back(G4TwoVector(-15.*cm, -30.*cm));
  
  G4ExtrudedSolid* xtruS1 
    = new G4ExtrudedSolid("XtruS1", polygon, 25.*cm, 
                 G4TwoVector(-20.*cm, 10.*cm), 1.5, G4TwoVector(), 0.5);
    
  G4ExtrudedSolid* xtruS2 
    = new G4ExtrudedSolid("XtruS2", polygon, 15.*cm, 
                 G4TwoVector(), 0.7, G4TwoVector(20.*cm, 20.*cm), 0.9);
  
  
  return new G4UnionSolid(
              "unionExtrudedSolid", xtruS1, xtruS2, 
              0, G4ThreeVector(0., 0., 40.*cm));
}              


//_____________________________________________________________________________
G4VSolid* createSolidAtlas1(std::vector<G4ThreeVector>& inside_pts,
                       std::vector<G4ThreeVector>& surface_pts,
                       std::vector<G4ThreeVector>& outside_pts)
{
  // Create Atlas problem shape

  G4TwoVector afllhpA( 2749.0        ,    0.0         ); 
  G4TwoVector afllhpB(   70.0        , 2412.182434674 );
  G4TwoVector afllhpC(   70.0        , 2431.021292888 ); 
  G4TwoVector afllhpD( 2769.922671698,    0.0         ); 

  // double tan42deg = 0.9004040442978399451204772038853717020764662112994852824270
  // double sin42deg = 0.6691306063588582138262733306867804735995832189597956768174
  // double cos42deg = 0.7431448254773942350146970489742569771891138734980263860401
  
  std::vector<G4TwoVector>  PointVector;
  PointVector.push_back( afllhpA ); 
  PointVector.push_back( afllhpB ); 
  PointVector.push_back( afllhpC ); 
  PointVector.push_back( afllhpD ); 
  
  G4TwoVector offSet1( 0.0,  0.0 ); 
  G4TwoVector offSet2( 0.0,  0.0 ); 
  G4double    dZ= 250.0;
  
  G4ExtrudedSolid* afslhpSolid
    = new G4ExtrudedSolid( "JFSH_AFrame_Leg_LowerHorizontalPlate",
			   PointVector,
			   dZ,		  // halfZ
			   offSet1, 1.0,	  	    
			   offSet2, 1.0); 	    

  inside_pts.push_back(G4ThreeVector( 2750.0*mm,     1.0*mm,     0.0*mm));  // inside [0]
  // inside_pts.push_back(G4ThreeVector( 2750.0*mm,     1.0*mm,  -233.0*mm));
  inside_pts.push_back(G4ThreeVector( 2750.0*mm,     1.0*mm,   249.0*mm));  // inside [1] 
  inside_pts.push_back(G4ThreeVector( 2767.9*mm,     1.0*mm,     0.0*mm));  // inside [2] 
  inside_pts.push_back(G4ThreeVector( 2767.9*mm,     1.0*mm,  -233.0*mm));  // inside [3] 
  // inside_pts.push_back(G4ThreeVector( 2767.9*mm,     1.0*mm,   249.0*mm)); 
  inside_pts.push_back(G4ThreeVector(   70.1*mm,  2412.182434674*mm, 0.0*mm));  // inside [4] 
  inside_pts.push_back(G4ThreeVector(   70.1*mm,  2430.000000000*mm, 0.0*mm));  // inside [5] 

  // inside_pts.push_back(G4ThreeVector(   70.1*mm,  2412.182434674*mm, 0.0*mm));  // inside [6] 
  
  // G4double  slopeE= afllhpC.y(); 
  const G4TwoVector   ABdir2v= afllhpB - afllhpA;
  const G4ThreeVector ABdir3v( ABdir2v.x(),  ABdir2v.y(), 0.0 ); 
  const G4TwoVector   ADdir2v= afllhpD - afllhpA;
  const G4ThreeVector ADdir3v( ADdir2v.x(),  ADdir2v.y(), 0.0 ); 

  const G4ThreeVector cornerAtop=    G4ThreeVector( afllhpA.x(), afllhpA.y(), dZ ); 
  const G4ThreeVector Zdirection=    G4ThreeVector( 0.0, 0.0, 1.0 ); 

  // G4double      fractionLine = 0.001;   // must be in    [ 0.0, 1.0 ] 
  G4ThreeVector pointOnEdge;
  pointOnEdge = cornerAtop + ABdir3v * 0.013579; 
  surface_pts.push_back( pointOnEdge );    //  No 0

  pointOnEdge = cornerAtop + ABdir3v * 0.1234567890; 
  surface_pts.push_back( pointOnEdge );    //  No 1

  pointOnEdge = cornerAtop + ADdir3v * 0.2345678901; 
  surface_pts.push_back( pointOnEdge );    //  No 2

  const G4TwoVector   ACdiag2v= afllhpC - afllhpA;
  const G4ThreeVector ACdiag3v( ACdiag2v.x(),  ACdiag2v.y(), 0.0 );   // Diagonal

  G4ThreeVector pointOnDiag= cornerAtop + ACdiag3v * 0.3456789;
  surface_pts.push_back( pointOnDiag );    //  No 3

  double      fraction = 0.3456789;    // must be in  [ 0.0, 1.0 ]
  double  factorHeight = 0.2030405;    // must be from ( -1.0 to +1.0 )
  const G4ThreeVector ApositionXY( afllhpA.x(), afllhpA.y(), 0.0 );   
  G4ThreeVector  pointOnSurf(0.0, 0.0, 0.0);
  pointOnSurf= ApositionXY + fraction * ABdir3v + factorHeight * dZ * Zdirection ; 
  surface_pts.push_back( pointOnSurf );    //  No 4
  //  if( 1 ) G4cout << " Surface point: " << pointOnSurf << std::endl;

  pointOnSurf= ApositionXY + fraction * ADdir3v + (1.0-factorHeight) * dZ * Zdirection ; 
  surface_pts.push_back( pointOnSurf );    //  No 5

  const G4ThreeVector CpositionXY( afllhpC.x(), afllhpC.y(), 0.0 );   
  // G4double  slopeE= afllhpC.y(); 
  const G4TwoVector   CDvec2v= afllhpD - afllhpC;
  const G4ThreeVector CDvec3v( CDvec2v.x(),  CDvec2v.y(), 0.0 ); 
  factorHeight = 0.405607809; 
  pointOnSurf= CpositionXY + fraction * CDvec3v - factorHeight       * dZ * Zdirection ; 
  surface_pts.push_back( pointOnSurf );    //  No 6

  G4ThreeVector BCvec2v= afllhpC - afllhpB;
  G4ThreeVector BCvec3v( BCvec2v.x(),  BCvec2v.y(), 0.0 ); 
  pointOnSurf= CpositionXY - fraction * BCvec3v + (1.0-factorHeight) * dZ * Zdirection ; 
  surface_pts.push_back( pointOnSurf );    //  No 7

  surface_pts.push_back( G4ThreeVector( afllhpD.x(), afllhpD.y(), 0.0*dZ ) );  // No 8 alt
  G4ThreeVector DcornerLowerMid( afllhpD.x(), afllhpD.y(), -0.5*dZ ); 
  surface_pts.push_back( DcornerLowerMid );    //  No 8

  // G4ThreeVector DcornerBottom( afllhpD.x(), afllhpD.y(), -dZ ); 
  // surface_pts.push_back( DcornerBottom );    //  No 9

  // G4ThreeVector pointOnSurface= cornerAtop  + fraction * ABdir2v; 
  // G4ThreeVector cornerAbottom= G4ThreeVector( afllhpA.x(), afllhpA.y(), -dZ ); 

  outside_pts.push_back(G4ThreeVector( 2748.99*mm,   0.0*mm,    0.0*mm));  // [0]
  outside_pts.push_back(G4ThreeVector( 2748.90*mm,   0.0*mm, -233.0*mm));  // [1]
  outside_pts.push_back(G4ThreeVector( 2748.0*mm,    0.0*mm,  249.0*mm));  // [2]
  outside_pts.push_back(G4ThreeVector( 2740.0*mm,    0.0*mm,  250.1*mm));  // [3]
  outside_pts.push_back(G4ThreeVector( 2748.9999*mm, 0.0*mm,  250.0*mm));  // [4] in top    plane
  outside_pts.push_back(G4ThreeVector( 2748.9999*mm, 0.0*mm, -250.0*mm));  // [5] in bottom plane
  outside_pts.push_back(G4ThreeVector( 2755.0*mm,   -1.0*mm,  125.0*mm));  // [6] below in Y

  G4ThreeVector outsidePt;
  factorHeight = -0.706050409; 
  outsidePt = ApositionXY + 1.002 * ABdir3v + factorHeight * dZ * Zdirection ; // [7]
  outside_pts.push_back(outsidePt);
  factorHeight =  0.500; 
  outsidePt = ApositionXY + 1.003 * ABdir3v + factorHeight * dZ * Zdirection 
                          + 0.015 * BCvec3v                                  ; // [8]
  outside_pts.push_back(outsidePt);

  outsidePt = CpositionXY + 0.010 * BCvec3v + 0.5          * dZ * Zdirection ; // [9]
  outside_pts.push_back(outsidePt);

  return afslhpSolid;
}                             

#if 0 
// EXTEND
// #include "AtlasExtrudedSolids.cc"

//_____________________________________________________________________________
G4VSolid* createSolidAtlas2(std::vector<G4ThreeVector>& inside_pts,
                       std::vector<G4ThreeVector>& surface_pts,
                       std::vector<G4ThreeVector>& outside_pts)
{
  // G4ExtrudedSolid *extrudedSolid=0; 
  G4VSolid  *ptrSolid= 0;
  G4double         dZ= -1.0;

  std::vector <G4TwoVector>  polygon;
  ptrSolid= CreateExtrudedSolidsFromAtlas( 2, polygon, dZ );  // Get 2nd solid


  return ptrSolid;  // extrudedSolid;
}
#endif


//_____________________________________________________________________________
G4VSolid* createSolid(G4int testCase,
                      std::vector<G4ThreeVector>& inside_points,
                      std::vector<G4ThreeVector>& surface_points,
                      std::vector<G4ThreeVector>& outside_points)
{
// Create selected test solid and fill vectors with points 

  switch ( testCase ) {
    case 0: return createSolid0(inside_points, surface_points, outside_points);                  
    case 1: return createSolid1(inside_points, surface_points, outside_points);                  
    case 2: return createSolid2(inside_points, surface_points, outside_points);                  
    case 3: return createSolid3(inside_points, surface_points, outside_points);                  
    case 4: return createSolid4(inside_points, surface_points, outside_points);
    case 5: return createSolidAtlas1( inside_points, surface_points, outside_points);
    // case 6: return createSolidAtlas2( inside_points, surface_points, outside_points);
    default: return 0;
  }
}                      

//_____________________________________________________________________________
void printResults(G4int testCase)
{
  
  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);

  // Set precision
  G4cout << std::setprecision(20) << G4endl;

  //
  // Test Inside
  //
  G4cout << " Inside = " << kInside << ",  kSurface= " << kSurface << ",  kOutside= " << kOutside << std::endl;

  for (G4int i=0; i<G4int(inside_points.size()); ++i) { 
    G4cout << i << "th inside_point Inside(p): " 
           << solid->Inside(inside_points[i]) << G4endl;
  } 
  G4cout << G4endl;

  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
    G4cout << i << "th surface_point Inside(p): " 
           << solid->Inside(surface_points[i]) << G4endl;
    EInside insideSurf = solid->Inside(surface_points[i]);
    if( insideSurf != kSurface ) 
	G4cerr << " Expected surface point " << surface_points[i] << " returned " << insideSurf << G4endl;
  } 
  G4cout << G4endl;

  for (G4int i=0; i<G4int(outside_points.size()); ++i) { 
    G4cout << i << "th outside_point Inside(p): " 
           << solid->Inside(outside_points[i]) << G4endl;
  } 
  G4cout << G4endl;

  //
  // Test DistanceToIn(p, v)
  //

  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
    G4cout << i << "th surface_point DistanceToIn(p, vx) " 
           << solid->DistanceToIn(surface_points[i], dirx) << G4endl;
    G4cout << i << "th surface_point DistanceToIn(p, -vx) " 
           << solid->DistanceToIn(surface_points[i], -dirx) << G4endl;
    G4cout << i << "th surface_point DistanceToIn(p, vy) " 
           << solid->DistanceToIn(surface_points[i], diry) << G4endl;
    G4cout << i << "th surface_point DistanceToIn(p, -vy) " 
           << solid->DistanceToIn(surface_points[i], -diry) << G4endl;
    G4cout << i << "th surface_point DistanceToIn(p, vz) " 
           << solid->DistanceToIn(surface_points[i], dirz) << G4endl;
    G4cout << i << "th surface_point DistanceToIn(p, -vz) " 
           << solid->DistanceToIn(surface_points[i], -dirz) << G4endl;
  }         
  G4cout << G4endl;

  for (G4int i=0; i<G4int(outside_points.size()); ++i) { 
    G4cout << i << "th outside_point DistanceToIn(p, vx) " 
           << solid->DistanceToIn(outside_points[i], dirx) << G4endl;
    G4cout << i << "th outside_point DistanceToIn(p, -vx) " 
           << solid->DistanceToIn(outside_points[i], -dirx) << G4endl;
    G4cout << i << "th outside_point DistanceToIn(p, vy) " 
           << solid->DistanceToIn(outside_points[i], diry) << G4endl;
    G4cout << i << "th outside_point DistanceToIn(p, -vy) " 
           << solid->DistanceToIn(outside_points[i], -diry) << G4endl;
    G4cout << i << "th outside_point DistanceToIn(p, vz) " 
           << solid->DistanceToIn(outside_points[i], dirz) << G4endl;
    G4cout << i << "th outside_point DistanceToIn(p, -vz) " 
           << solid->DistanceToIn(outside_points[i], -dirz) << G4endl;
  }         
  G4cout << G4endl;

  //
  // Test DistanceToOut(p, v) function.
  //
           
  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
    G4cout << i << "th surface_point DistanceToOut(p, vx) " 
           << solid->DistanceToOut(surface_points[i], dirx, false, 0, 0) << G4endl;
    G4cout << i << "th surface_point DistanceToOut(p, -vx) " 
           << solid->DistanceToOut(surface_points[i], -dirx, false, 0, 0) << G4endl;
    G4cout << i << "th surface_point DistanceToOut(p, vy) " 
           << solid->DistanceToOut(surface_points[i], diry, false, 0, 0) << G4endl;
    G4cout << i << "th surface_point DistanceToOut(p, -vy) " 
           << solid->DistanceToOut(surface_points[i], -diry, false, 0, 0) << G4endl;
    G4cout << i << "th surface_point DistanceToOut(p, vz) " 
           << solid->DistanceToOut(surface_points[i], dirz, false, 0, 0) << G4endl;
    G4cout << i << "th surface_point DistanceToOut(p, -vz) " 
           << solid->DistanceToOut(surface_points[i], -dirz, false, 0, 0) << G4endl;
  }         
  G4cout << G4endl;

  for (G4int i=0; i<G4int(inside_points.size()); ++i) { 
    G4cout << i << "th inside_point DistanceToOut(p, vx) " 
           << solid->DistanceToOut(inside_points[i], dirx, false, 0, 0) << G4endl;
    G4cout << i << "th inside_point DistanceToOut(p, -vx) " 
           << solid->DistanceToOut(inside_points[i], -dirx, false, 0, 0) << G4endl;
    G4cout << i << "th inside_point DistanceToOut(p, vy) " 
           << solid->DistanceToOut(inside_points[i], diry, false, 0, 0) << G4endl;
    G4cout << i << "th inside_point DistanceToOut(p, -vy) " 
           << solid->DistanceToOut(inside_points[i], -diry, false, 0, 0) << G4endl;
    G4cout << i << "th inside_point DistanceToOut(p, vz) " 
           << solid->DistanceToOut(inside_points[i], dirz, false, 0, 0) << G4endl;
    G4cout << i << "th inside_point DistanceToOut(p, -vz) " 
           << solid->DistanceToOut(inside_points[i], -dirz, false, 0, 0) << G4endl;
  }         
  G4cout << G4endl;
  
  //
  // Test surface area
  //
  G4cout << "Surface: " << solid->GetSurfaceArea() << G4endl;
  G4cout << G4endl;

  //
  // Test volume
  //
  G4int ntrial = 10;
  G4double sum = 0;
  G4double min = DBL_MAX;
  G4double max = -DBL_MAX;
  
  G4cout << "Evaluating volume ..." << G4endl;
  G4Timer time;
  time.Start();
  for (G4int i=0; i<ntrial; ++i ) {
    G4VSolid* solid0 = createSolid(testCase, inside_points, surface_points, outside_points);
    G4double value = solid0->GetCubicVolume();
    sum += value;
    if ( value < min ) min = value;
    if ( value > max ) max = value;
    delete solid0;
  }
  time.Stop();
  G4cout << "Average volume after " << ntrial << " trials: " << sum/ntrial << G4endl;
  G4cout << "               in the interval: " << max-min << G4endl;
  G4cout << "Time taken was: " << time.GetRealElapsed() << " seconds. " << G4endl;
  G4cout << G4endl;
 
  delete solid;
}
  

//_____________________________________________________________________________
void testInside(G4int testCase)
{
// Test Inside

  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);
  
  //
  // Test Inside
  //

  for (G4int i=0; i<G4int(inside_points.size()); ++i) { 
    assert( solid->Inside(inside_points[i]) == kInside );
  } 

  // if( testCase == 5 ) G4cout << " Surface points " << G4endl;
  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
    // if( testCase == 5 ) G4cout << " point " << i << G4endl;
    assert( solid->Inside(surface_points[i]) == kSurface );
  } 

  for (G4int i=0; i<G4int(outside_points.size()); ++i) { 
    assert( solid->Inside(outside_points[i]) == kOutside );
  } 
  
  delete solid;
}

//_____________________________________________________________________________
void testDistanceToInPV(G4int testCase)
{
// Test DistanceToIn  

  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);
  
  if ( testCase == 0 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirx) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[4], -dirx) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[4], -diry) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[5],  diry) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[6], -diry) - 200.0 )< kCarTolerance );
  
    // For points on surface we get testDistanceToIn = 9e+99, is it ok ?
  }  
  else if ( testCase == 1 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirx) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[4], -dirx) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[5],  diry) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[6], -diry) - 100.0 )< kCarTolerance );
  
    // Add points on surface
  }  
  else if ( testCase == 2 ) {
    //assert( std::fabs(solid->DistanceToIn(outside_points[0], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1],  dirz) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirx) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  diry) - 200.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3], -dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[4], -dirx) - 180.0 )< kCarTolerance );
  
    // Add points on surface
  }  
  else if ( testCase == 3 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1],  dirz) - 88.4615384615384670045 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirz) - 500.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirx) - 155.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -dirx) -  55.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  diry) -  55.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirz) -  80.882352941176478112 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3], -dirz) -  50.0 )< kCarTolerance );
 
    // Add points on surface
  }  
  else if ( testCase == 5 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirx) - 0.01 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[0], -dirx) == kInfinity );


    
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  diry) - 0.01*std::tan(42.*deg) )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[0], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[0],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[0], -dirz) == kInfinity );

    assert( std::fabs(solid->DistanceToIn(outside_points[1],  dirx) - 0.1 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[1], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(outside_points[1],  diry) - 0.1*std::tan(42.*deg) )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[1], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[1],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[1], -dirz) == kInfinity );

    assert( std::fabs(solid->DistanceToIn(outside_points[2],  dirx) - 1.0 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[2], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(outside_points[2],  diry) - 1.0*std::tan(42.*deg) )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[2], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[2],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[2], -dirz) == kInfinity );

    assert( solid->DistanceToIn(outside_points[3], -dirx) == kInfinity );
    assert( solid->DistanceToIn(outside_points[3],  dirx) == kInfinity );
    assert( solid->DistanceToIn(outside_points[3], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[3],  diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[3], -dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[3],  dirz) == kInfinity );
    G4ThreeVector vectorToEdge = surface_points[0] - outside_points[3];
    G4ThreeVector dirToEdge= vectorToEdge.unit();

#ifdef EXTRA_PRINTS
    int idPt=  -1;
    G4double distIn=    -9.0;
    G4double expected= -10.0;
#endif
#if 0
    int idPt= 3;
    G4cout << "Checking 'outside' point = " << outside_points[idPt] << G4endl;
    G4cout << "         direction       = " << dirToEdge << G4endl;
    G4cout << "       aimed at surf point " << surface_points[0] << G4endl;
    G4double distIn= solid->DistanceToIn(outside_points[idPt], dirToEdge);
    G4cout << " Evaluating DistanceToIn ( pt, +y ) = " << G4endl
	   << " Obtained     = " << distIn << G4endl;
    G4double expected = vectorToEdge.mag(); 
    G4cout << " Expected     = " << expected << G4endl;
    G4cout << " Difference   = " << distIn-expected << " absolute (obtain-expect) " << G4endl;
    G4cout << " Relative diff= " << (distIn-expected)/expected << " (obtain-expect)/expect " << G4endl;
#endif

    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirToEdge) - vectorToEdge.mag() ) < kCarTolerance );
    G4ThreeVector vectorToDiag = surface_points[3]- outside_points[3];
    G4ThreeVector dirToDiag= vectorToDiag.unit();
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirToDiag) - vectorToDiag.mag() ) < kCarTolerance );

    assert( std::fabs(solid->DistanceToIn(outside_points[4],  dirx) - 0.0001 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[4], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(outside_points[4],  diry) - 0.0001*std::tan(42.*deg) )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[4], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[4],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[4], -dirz) == kInfinity );

    assert( std::fabs(solid->DistanceToIn(outside_points[5],  dirx) - 0.0001 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[5], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(outside_points[5],  diry) - 0.0001*std::tan(42.*deg) )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[5], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[5],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[5], -dirz) == kInfinity );

    assert( solid->DistanceToIn(outside_points[6],  dirx) == kInfinity );
    assert( solid->DistanceToIn(outside_points[6], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(outside_points[6],  diry) - 1.0 )< kCarTolerance );

    assert( solid->DistanceToIn(outside_points[6], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[6],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[6], -dirz) == kInfinity );

    double  sin45= std::sin(45.0*deg);   //    = 1.0 / std::sqrt( 2.0 ); 
    G4ThreeVector dirXYplus( sin45, sin45, 0.0 ); 
    assert( std::fabs(solid->DistanceToIn(outside_points[6],  dirXYplus) - 1.0/sin45 )< kCarTolerance );
    // G4ThreeVector DcornerBottom= surface_points[8];
    G4ThreeVector DcornerLowerMid= surface_points[8]; 
    G4ThreeVector DcornerBottom=   DcornerLowerMid + dirz * DcornerLowerMid.z();
    G4ThreeVector vecToDbot=  outside_points[6] - DcornerBottom;
    assert( solid->DistanceToIn(outside_points[6],  vecToDbot.unit() ) > vecToDbot.mag() - kCarTolerance );

    assert( std::fabs(solid->DistanceToIn(outside_points[7],  dirx) - 0.002*2679.0 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[7], -dirx) == kInfinity );
    assert( solid->DistanceToIn(outside_points[7],  diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[7], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[7],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[7], -dirz) == kInfinity );
    G4ThreeVector dirX42mY( std::cos(42.0*deg), -std::sin(42.0*deg) );      // along B->A
    // Check the line that should graze solid is seen to enter 
    assert( std::fabs(solid->DistanceToIn(outside_points[7],  dirX42mY) 
		      - 0.002*2679.0/std::cos(42.0*deg)                 )< kCarTolerance );

    assert( std::fabs(solid->DistanceToIn(outside_points[8],  dirx) - 0.003*2679.0 )< kCarTolerance );
    assert( solid->DistanceToIn(outside_points[8], -dirx) == kInfinity );
    assert( solid->DistanceToIn(outside_points[8],  diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[8], -diry) == kInfinity );
    assert( solid->DistanceToIn(outside_points[8],  dirz) == kInfinity );
    assert( solid->DistanceToIn(outside_points[8], -dirz) == kInfinity );
    
    // assert( std::fabs(solid->DistanceToIn(outside_points[8],  dirx) - 0.003*2679.0 )< kCarTolerance );
 

    assert(     solid->DistanceToIn(outside_points[9],  dirx) == kInfinity );
    assert(     solid->DistanceToIn(outside_points[9], -dirx) == kInfinity );
    assert(     solid->DistanceToIn(outside_points[9],  diry) == kInfinity );
    assert( 
      std::fabs(solid->DistanceToIn(outside_points[9], -diry) - 0.010*14./std::cos(42.*deg) )< kCarTolerance );
    assert(     solid->DistanceToIn(outside_points[9],  dirz) == kInfinity );
    assert(     solid->DistanceToIn(outside_points[9], -dirz) == kInfinity );

    // Add points on surface
    // ---------------------
    assert( std::fabs(solid->DistanceToIn(surface_points[0],  dirx) - 0.0 )< kCarTolerance );
    assert( solid->DistanceToIn(surface_points[0], -dirx) == kInfinity );
    assert( std::fabs(solid->DistanceToIn(surface_points[0],  diry) - 0.0 )< kCarTolerance );
    assert( solid->DistanceToIn(surface_points[0], -diry) == kInfinity );
    assert( solid->DistanceToIn(surface_points[0],  dirz) == kInfinity );

#if 0
    // G4double distIn, expected; 
    idPt= 0;
    G4cout << "Checking 'surface' point [" << idPt << "]= " << surface_points[idPt] << G4endl;
    G4cout << "         direction          = " << -dirz << G4endl;   // dirToEdge << G4endl;
    // G4cout << "       aimed at surf point " << surface_points[0] << G4endl;
    distIn= solid->DistanceToIn(surface_points[idPt], -dirz); 
    G4cout << " Evaluating DistanceToIn ( pt, -z ) = " << G4endl
	   << " Obtained     = " << distIn << G4endl;
    expected = 0.0;
    G4cout << " Expected     = " << expected << G4endl;
    G4cout << " Difference   = " << distIn-expected << " absolute (obtain-expect) " << G4endl;
    // G4cout << " Relative diff= " << (distIn-expected)/expected << " (obtain-expect)/expect " << G4endl;

    assert( std::fabs(solid->DistanceToIn(surface_points[0], -dirz) - 0.0 )< kCarTolerance );
    //   Fails !!
#endif

    assert( std::fabs(solid->DistanceToIn(surface_points[0], 
					  G4ThreeVector(0.6, 0.8, 0.0) ) 
		      - 0.0 ) < kCarTolerance );
    assert( solid->DistanceToIn(surface_points[0], G4ThreeVector(0.6, -0.8, 0.0) ) 
	       == kInfinity );
  }  

  delete solid;
}


//_____________________________________________________________________________
void testDistanceToOutPV(G4int testCase)
{
// Test DistanceToOutPV  

  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);
  
  if ( testCase == 0 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) -  75.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 175.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 275.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) -  75.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 150.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 1 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) - 200.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 150.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 2 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 500.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 500.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 125.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 110.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 410.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 430.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) -  10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) -  80.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 340.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) -  20.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 450.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 3 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 110.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 410.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 430.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) -  70.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) -  35.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 130.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 290.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 141.66666666666668561 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) -  63.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) -  24.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 174.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 174.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 137.1428571428571388 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) -  20.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3],  dirx) -  16.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3], -dirx) - 107.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3],  diry) - 416.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3], -diry) -  76.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3],  dirz) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3], -dirz) -  15.384615384615381473 )< kCarTolerance );

    // For Add points on surface
  }  
  else if ( testCase == 5 ) {
    double tan42 = 0.9004040442978399451204772038853717020764662112994852824270;
    double sin42 = 0.6691306063588582138262733306867804735995832189597956768174;
    double cos42 = 0.7431448254773942350146970489742569771891138734980263860401;
    const double x0plus=  ( 14.0/cos42 - 1. ) / tan42 - 1.0;
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - x0plus )< kCarTolerance );
    const double x0minus=   14.0/sin42 - x0plus;
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - x0minus )< kCarTolerance );
    const double y0plus=  ( 14.0/sin42 - 1.0) * tan42 - 1.0; 
    assert( std::fabs(solid->DistanceToOut(inside_points[0],   diry) - y0plus  )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  -diry) -   1.0 )< kCarTolerance );
    const double dZ= 250.0;
    assert( std::fabs(solid->DistanceToOut(inside_points[0],   dirz) - dZ )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  -dirz) - dZ )< kCarTolerance );

    // assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 0.1 )< kCarTolerance );
    // assert( solid->DistanceToOut(inside_points[1], -dirx) == KInfinity );
    // assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 0.1*tan(42.*deg) )< kCarTolerance );
    // assert( solid->DistanceToOut(inside_points[1], -diry) == KInfinity );
    // G4cout << " DistToOut( insidePt1, +z ) = " << solid->DistanceToOut(inside_points[1],  dirz) << std::endl;

    assert( solid->DistanceToOut(inside_points[1],  dirz) == 1.0 );
    assert( solid->DistanceToOut(inside_points[1], -dirz) == (2.0*dZ-1.0) );

    assert( std::fabs(solid->DistanceToOut(inside_points[2],  -diry) - 1.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],   dirz) - dZ )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  -dirz) - dZ )< kCarTolerance );
    G4ThreeVector  cornerAmid( 2749.0,   0.0,  0.0  ); 
    G4ThreeVector  vecIn2toA = cornerAmid - inside_points[2] ; 
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  vecIn2toA.unit()) - vecIn2toA.mag() )< kCarTolerance );

    G4ThreeVector vectorToEdge = surface_points[0] - inside_points[3];
    G4ThreeVector dirToEdge= vectorToEdge.unit();

#if 0
    int idPt= 3;
    G4cout << "Checking 'inside' point = " << inside_points[idPt] << G4endl;
    G4cout << "         direction       = " << dirToEdge << G4endl;
    G4cout << "       aimed at surf point " << surface_points[0] << G4endl;
    G4double distIn= solid->DistanceToOut(inside_points[idPt], dirToEdge);
    G4cout << " Evaluating DistanceToOut ( pt, +y ) = " << G4endl
	   << " Obtained     = " << distIn << G4endl;
    G4double expected = vectorToEdge.mag(); 
    G4cout << " Expected     = " << expected << G4endl;
    G4cout << " Difference   = " << distIn-expected << " absolute (obtain-expect) " << G4endl;
    G4cout << " Relative diff= " << (distIn-expected)/expected << " (obtain-expect)/expect " << G4endl;
#endif 

    assert( std::fabs(solid->DistanceToOut(inside_points[3],  dirToEdge) - vectorToEdge.mag() ) < kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3],  dirz) - (233.0+dZ) ) < kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[3], -dirz) - (dZ-233.0) ) < kCarTolerance );

    G4ThreeVector vectorToDiag = surface_points[3]- inside_points[3];
    G4ThreeVector dirToDiag= vectorToDiag.unit();
    assert( std::fabs(solid->DistanceToOut(inside_points[3],  dirToDiag) - vectorToDiag.mag() ) < kCarTolerance );

    assert( std::fabs(solid->DistanceToOut(inside_points[4],  dirx) - ( (14.0/sin42) - 0.1) )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[4], -dirx) - 0.1 )< kCarTolerance );
    double  distOut_ptIn4Yplus  = (14.0/sin42 - 0.1) * tan42;  //  - 0.1;
    double  distOut_ptIn4Yminus = 0.1 * tan42; // (14.0/cos42) - distOut_ptIn4Yplus;

    // G4cout << " DistToOut ( ptIn4, +y ) = " << solid->DistanceToOut(inside_points[4],  diry) << std::endl;
    // G4cout << " Expected                = " << distOut_ptIn4Yplus << std::endl;

    assert( std::fabs(solid->DistanceToOut(inside_points[4],  diry) - distOut_ptIn4Yplus)  < 1.0e-04); // kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[4], -diry) - distOut_ptIn4Yminus) < kCarTolerance );
	    // ( (14.0/sin42 - 0.1) * tan42 -0.1 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[4],   dirz) - dZ )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[4],  -dirz) - dZ )< kCarTolerance );

    // assert( std::fabs(solid->DistanceToIn(inside_points[8],  dirx) - 0.003*2679.0 )< kCarTolerance );
 
    // Add points on surface
    // ---------------------
    assert( std::fabs(solid->DistanceToOut(surface_points[0],  dirx) - 14.0/sin42 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[0], -dirx) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[0],  diry) - 14.0/cos42 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[0], -diry) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[0],  dirz) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[0], -dirz) - 2.*dZ )< kCarTolerance );

    assert( std::fabs(solid->DistanceToOut(surface_points[1],  dirx) - 14.0/sin42 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[1], -dirx) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[1],  diry) - 14.0/cos42 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[1], -diry) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[1],  dirz) - 0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[1], -dirz) - 2.*dZ )< kCarTolerance );

    G4ThreeVector surfPt2= surface_points[2]; 
    double s2x = ( surfPt2.x() - 2749.0 );
    // double d2x = s2x ; 
    assert( std::fabs(solid->DistanceToOut(surface_points[2],  dirx) - (14.0/sin42-s2x) )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[2], -dirx) - s2x )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[2],  diry) - (14.0/sin42-s2x)*tan42 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[2], -diry) -  0.0  )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[2],  dirz) -  0.0  )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(surface_points[2], -dirz) - 2.*dZ )< kCarTolerance );

    G4ThreeVector surfPt3= surface_points[3];
    // pointOnEdge = cornerAtop + ADdir3v * 0.2345678901;  // Surface points
    // G4ThreeVector vecAtoSu3 = ( surfPt3 - G4ThreeVector( 2749.0, 0.0, 0.0 ) ) ;
    // double s3x = ( surfPt3.x() - 2749.0 );

    //****   Continue here ...   JA 2011.09.29
#if 0
    // G4double distIn, expected; 
    idPt= 0;
    G4cout << "Checking 'surface' point [" << idPt << "]= " << surface_points[idPt] << G4endl;
    G4cout << "         direction          = " << -dirz << G4endl;   // dirToEdge << G4endl;
    // G4cout << "       aimed at surf point " << surface_points[idPt] << G4endl;
    distIn= solid->DistanceToOut(surface_points[idPt], -dirz); 
    G4cout << " Evaluating DistanceToOut ( pt, -z ) = " << G4endl
	   << " Obtained     = " << distIn << G4endl;
    expected = 0.0;
    G4cout << " Expected     = " << expected << G4endl;
    G4cout << " Difference   = " << distIn-expected << " absolute (obtain-expect) " << G4endl;
    // G4cout << " Relative diff= " << (distIn-expected)/expected << " (obtain-expect)/expect " << G4endl;
#endif 

//    assert( std::fabs(solid->DistanceToOut(surface_points[0], 
//					  G4ThreeVector(0.6, 0.8, 0.0) ) 
//		      - 0.0 ) < kCarTolerance );
//    assert( solid->DistanceToOut(surface_points[0], G4ThreeVector(0.6, -0.8, 0.0) ) 
//	       == kInfinity );
  }  



  delete solid;
}

//_____________________________________________________________________________
void testSurface(G4int testCase) 
{
// Test surface

  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);
  
  if ( testCase == 0 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 1524984.4718999243341 ) < 1e-6 );
  }
  if ( testCase == 1 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 2160000 ) < 1e-6 );
  }
  if ( testCase == 2 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 2506922.4391292142682 ) < 1e-6 );
  }
  if ( testCase == 3 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 3638542.9775616745465 ) < 1e-6 );
  }

  delete solid;
}  
    
//_____________________________________________________________________________
void testVolume(G4int testCase) 
{
// Test volume
// The volume is evaluated via G4VSolid, that's why the precision is very low

  std::vector<G4ThreeVector> inside_points;
  std::vector<G4ThreeVector> surface_points;
  std::vector<G4ThreeVector> outside_points;
  G4VSolid* solid = createSolid(testCase, inside_points, surface_points, outside_points);
  
  if ( testCase == 0 ) {
    assert( std::fabs(solid->GetCubicVolume() - 108.0e+6 ) < 1.0e+6 );
  }
  if ( testCase == 1 ) {
    assert( std::fabs(solid->GetCubicVolume() - 216.0e+6 ) < 1.0e+6 );
  }
  if ( testCase == 2 ) {
    assert( std::fabs(solid->GetCubicVolume() - 121.7e+6 ) < 1.0e+6 );
  }
  
  if ( testCase == 3 ) {
    assert( std::fabs(solid->GetCubicVolume() - 162.1e+6 ) < 1.0e+6 );
  }

  delete solid;
}  
    

//_____________________________________________________________________________
int main()
{
  // Uncomment this line to print the results for a tested solid case
  // printResults(5);

  for  ( G4int testCase = 0; testCase < 6; ++testCase ) { 
    testInside(testCase);
    testDistanceToInPV(testCase);
    testDistanceToOutPV(testCase);
    testSurface(testCase);
    testVolume(testCase);
  }
  
  return 0;
}
