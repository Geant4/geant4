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
// $Id: testG4ExtrudedSolid.cc,v 1.4 2008-02-27 12:33:20 ivana Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include <assert.h>
#include <cmath>
#include <vector>
#include <iomanip>

#include "globals.hh"

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

  G4TessellatedSolid* tessellated
    = new G4TessellatedSolid(*extruded);
   return tessellated;
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

  for (G4int i=0; i<G4int(inside_points.size()); ++i) { 
    G4cout << i << "th inside_point Inside(p): " 
           << solid->Inside(inside_points[i]) << G4endl;
  } 
  G4cout << G4endl;

  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
    G4cout << i << "th surface_point Inside(p): " 
           << solid->Inside(surface_points[i]) << G4endl;
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

  for (G4int i=0; i<G4int(surface_points.size()); ++i) { 
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
  // printResults(1);

  for  ( G4int testCase = 0; testCase < 4; ++testCase ) { 
    testInside(testCase);
    testDistanceToInPV(testCase);
    testDistanceToOutPV(testCase);
    testSurface(testCase);
    testVolume(testCase);
  }
  
  return 0;
}
