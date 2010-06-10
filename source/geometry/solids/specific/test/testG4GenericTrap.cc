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
//
//
// testG4GenericTrap
//
// Test file for class G4GenericTrap, derived from testG4ExtrudedSolid 
// In the functions createSolidN(...), there are defined several 
// test cases of arbitrary trapezoid solid. These functions also fill in 
// the vectors with explicitely defined points inside, on the surface 
// and outside the solid.
// All the test results for the definesd solids and
// points can be printed via PrintResults() function.
// The tests are then defined in testXYZ() functions
// using assert() on the comparison with the expected*
// result value.
//
// The functions DistanceToIn, DistanceToOut on surface
// point do not give always expected values, that's why
// they are not yet included in the tests with assert.
// To be added tests for SurfaceNormal(p) function.
// Ensure asserts are compileed in.
//
// Author:
//   Ivana Hrivnacova, IPN Orsay*
//
// Adapted for G4GenericTrap by  Tatiana Nikitina 

#include <assert.h>
#include <cmath>
#include <vector>
#include <iomanip>

#include "globals.hh"

#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4GenericTrap.hh"
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
  // 1 down vertex: 0=1=2=3  
  inside_points.push_back(G4ThreeVector(  0.0*cm, -45.0*cm, 70.0*cm));
  
  surface_points.push_back(G4ThreeVector( 0.0*cm, -45.0*cm,  75.0*cm));

  outside_points.push_back(G4ThreeVector(  0.0*cm, -45.0*cm,  80.0*cm));
  outside_points.push_back(G4ThreeVector( 50.0*cm, -20.0*cm,  30.0*cm));
  outside_points.push_back(G4ThreeVector(-40.0*cm, -20.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector(  0.0*cm,   0.0*cm,  30.0*cm));
  outside_points.push_back(G4ThreeVector(-80.0*cm,   0.0*cm, -30.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( 15.*cm, -15.*cm));
  
  return new G4GenericTrap("genTrap0", 75.*cm, vertices);
}                             

//_____________________________________________________________________________
G4VSolid* createSolid1(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 2 down vertices: 0=1 2=3

  inside_points.push_back(G4ThreeVector( -15.0*cm, -74.0 *cm, 74.0 *cm));
  inside_points.push_back(G4ThreeVector(  -5.0*cm, -74.0*cm,  70.0*cm));  
  
  surface_points.push_back(G4ThreeVector(   0.0*cm,  -75*cm,     0.0*cm));
  surface_points.push_back(G4ThreeVector(  35.0*cm,  -15.0*cm,  75.0*cm));
  surface_points.push_back(G4ThreeVector( -15.0*cm,  -75.0*cm, -75.0*cm));
 

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  0.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm,  5.0*cm,  40.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm, -5.0*cm, -40.0*cm));
  outside_points.push_back(G4ThreeVector(-35.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector(+35.0*cm,  0.0*cm, -10.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm,-40.0*cm, -20.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm, 40.0*cm,  10.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector(   0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(   0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( -30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( -30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(  45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(   0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( -30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(  15.*cm, -15.*cm));
  
  return new G4GenericTrap("genTrap1", 75.*cm, vertices);
}                             

//_____________________________________________________________________________
G4VSolid* createSolid2(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 1 up vertex: 4=5=6=7


   
    // Add points on surface
  inside_points.push_back(G4ThreeVector( 30.0*cm, -20.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector( 10.0*cm, -50.0*cm, -20.0*cm));
  inside_points.push_back(G4ThreeVector( 44.0*cm, -16.0*cm,  20.0*cm));  
  
  surface_points.push_back(G4ThreeVector( 45.0*cm, -15.0*cm, -25.0*cm));
  surface_points.push_back(G4ThreeVector( 45.0*cm, -15.0*cm,  75.0*cm));
  surface_points.push_back(G4ThreeVector(-10.0*cm, -75.0*cm, -75.0*cm));
 

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,  30.0*cm));
  outside_points.push_back(G4ThreeVector( 10.0*cm,  5.0*cm, -40.0*cm));
  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,   0.0*cm));
  outside_points.push_back(G4ThreeVector(-40.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( 40.0*cm,  0.0*cm, -10.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( 15.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  
  return new G4GenericTrap("genTrap2", 75.*cm, vertices);
}                             


//_____________________________________________________________________________
G4VSolid* createSolid3(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 2 up vertices: 4=5 6=7

  inside_points.push_back(G4ThreeVector( -10.0*cm,  -50.0*cm, -35.0*cm));
  inside_points.push_back(G4ThreeVector( -10.0*cm,  -74.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector( -15.0*cm,  -74.0*cm, -74.0*cm));
    
  surface_points.push_back(G4ThreeVector( 0.0*cm,  -75.0*cm,   0.0*cm));
  surface_points.push_back(G4ThreeVector( 15.0*cm, -75.0*cm,  75.0*cm));
  surface_points.push_back(G4ThreeVector( -15.0*cm, -75.0*cm,   -75.0*cm));

  outside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -50.0*cm));
  outside_points.push_back(G4ThreeVector( 25.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( -5.0*cm,  5.0*cm,  15.0*cm));
  outside_points.push_back(G4ThreeVector( 45.0*cm, 40.0*cm,  45.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( 15.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  
  return new G4GenericTrap("genTrap3", 75.*cm, vertices);
}  
  

//_____________________________________________________________________________
G4VSolid* createSolid4(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 4 down vertex: 
// 4 up vertices: Box like Solid

  inside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm,    0.0*cm));
  inside_points.push_back(G4ThreeVector(-10.0*cm,  0.0*cm,    5.0*cm));
  inside_points.push_back(G4ThreeVector( 15.0*cm,-15.0*cm,   74.0*cm));
    
  surface_points.push_back(G4ThreeVector( 45.0*cm, 45.0*cm, -75.0*cm));
  surface_points.push_back(G4ThreeVector( 15.0*cm, 15.0*cm,  75.0*cm));
  surface_points.push_back(G4ThreeVector(-15.0*cm, -15.0*cm, 75.0*cm));

  outside_points.push_back(G4ThreeVector(  0.0*cm,  0.0*cm, -80.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm,  0.0*cm,  80.0*cm));
  outside_points.push_back(G4ThreeVector( 55.0*cm,  0.0*cm,  15.0*cm));
  outside_points.push_back(G4ThreeVector(  0.0*cm, 55.0*cm, -15.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( -45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( -45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( -45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( -45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm, -45.*cm));
  
  return new G4GenericTrap("genTrap4", 75.*cm, vertices);
}              


//_____________________________________________________________________________
G4VSolid* createSolid5(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 4 down vertices:  
// 4 up vertex: Trap like Solid 
  inside_points.push_back(G4ThreeVector(   0.0*cm,  0.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector( -10.0*cm,  0.0*cm,   5.0*cm));
  inside_points.push_back(G4ThreeVector(  15.0*cm,-15.0*cm,  74.0*cm));
    
  surface_points.push_back(G4ThreeVector( 35.0*cm, 35.0*cm,   75.0*cm));
  surface_points.push_back(G4ThreeVector( 35.0*cm,-35.0*cm,  -75.0*cm));
  surface_points.push_back(G4ThreeVector( 30.0*cm, 30.0*cm,   75.0*cm));

  outside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -50.0*cm));
  outside_points.push_back(G4ThreeVector( 55.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( 0.0*cm,  5.0*cm,  80.0*cm));
  outside_points.push_back(G4ThreeVector( 5.0*cm, 0.0*cm,  -80.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( -45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( -45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector(  45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( -35.*cm, -35.*cm));
  vertices.push_back(G4TwoVector( -35.*cm,  35.*cm));
  vertices.push_back(G4TwoVector(  35.*cm,  35.*cm));
  vertices.push_back(G4TwoVector(  35.*cm, -35.*cm));
  
  
  return new G4GenericTrap("arbTrap5", 75.*cm, vertices);
}              

//_____________________________________________________________________________
G4VSolid* createSolid6(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
  // all up & down vertices different, twisted
  inside_points.push_back(G4ThreeVector(   0.0*cm,  0.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector( -10.0*cm,  0.0*cm,   5.0*cm));
  inside_points.push_back(G4ThreeVector(  15.0*cm, 15.0*cm,  74.0*cm));
    
  surface_points.push_back(G4ThreeVector( 35.0*cm, 35.0*cm,   75.0*cm));
  surface_points.push_back(G4ThreeVector( 35.0*cm,-35.0*cm,  -75.0*cm));
  surface_points.push_back(G4ThreeVector( 30.0*cm, 30.0*cm,   75.0*cm));

  outside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -50.0*cm));
  outside_points.push_back(G4ThreeVector( 55.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector(  0.0*cm,  5.0*cm,  80.0*cm));
  outside_points.push_back(G4ThreeVector(  5.0*cm, 0.0*cm,  -80.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector(-45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector(-45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector( 45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector(-35.*cm, -35.*cm));
  vertices.push_back(G4TwoVector(-35.*cm,  35.*cm));
  vertices.push_back(G4TwoVector( 35.*cm,  35.*cm));
  vertices.push_back(G4TwoVector( 35.*cm,  15.*cm));
  
  return new G4GenericTrap("genTrap6", 75.*cm, vertices);
}              

//_____________________________________________________________________________
G4VSolid* createSolid7(std::vector<G4ThreeVector>& inside_points,
                       std::vector<G4ThreeVector>& surface_points,
                       std::vector<G4ThreeVector>& outside_points)
{
// 3 up vertices: 7=8 (twisted )
 
  inside_points.push_back(G4ThreeVector(   0.0*cm,  0.0*cm,   0.0*cm));
  inside_points.push_back(G4ThreeVector( -10.0*cm,  0.0*cm,   5.0*cm));
  inside_points.push_back(G4ThreeVector(  15.0*cm,-15.0*cm,  74.0*cm));
    
  surface_points.push_back(G4ThreeVector( 35.0*cm, 35.0*cm,   75.0*cm));
  surface_points.push_back(G4ThreeVector( 35.0*cm,-35.0*cm,  -75.0*cm));
  surface_points.push_back(G4ThreeVector( 30.0*cm, 30.0*cm,   75.0*cm));

  outside_points.push_back(G4ThreeVector(-50.0*cm, 10.0*cm, -50.0*cm));
  outside_points.push_back(G4ThreeVector( 55.0*cm,  0.0*cm,  10.0*cm));
  outside_points.push_back(G4ThreeVector( 0.0*cm,  5.0*cm,  80.0*cm));
  outside_points.push_back(G4ThreeVector( 5.0*cm, 0.0*cm,  -80.0*cm));

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector(-45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector(-45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector( 45.*cm,  45.*cm));
  vertices.push_back(G4TwoVector( 45.*cm, -45.*cm));
  vertices.push_back(G4TwoVector( 35.*cm, -35.*cm));
  vertices.push_back(G4TwoVector(-35.*cm, -35.*cm));
  vertices.push_back(G4TwoVector( 35.*cm,  35.*cm));
  vertices.push_back(G4TwoVector( 35.*cm,  35.*cm));
  
  return new G4GenericTrap("genTrap7", 75.*cm, vertices);
}              

//_____________________________________________________________________________
G4VSolid* createSolid8(std::vector<G4ThreeVector>& /*inside_points*/,
                       std::vector<G4ThreeVector>& /*surface_points*/,
                       std::vector<G4ThreeVector>& /*outside_points*/)
{
// 3 up vertices: 4=5 (twisted )

  std::vector<G4TwoVector> vertices;
  vertices.push_back(G4TwoVector( 45.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( 15.*cm, -15.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(  0.*cm, -75.*cm));
  vertices.push_back(G4TwoVector(-30.*cm, -75.*cm));
  vertices.push_back(G4TwoVector( 15.*cm, -15.*cm));
  
  return new G4GenericTrap("genTrap8", 75.*cm, vertices);
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
    case 5: return createSolid5(inside_points, surface_points, outside_points);                  
    case 6: return createSolid6(inside_points, surface_points, outside_points);                  
    case 7: return createSolid7(inside_points, surface_points, outside_points);                  
    case 8: return createSolid8(inside_points, surface_points, outside_points);                  
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
  //assert( false );


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
  //assert(false);
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
    assert( std::fabs(solid->DistanceToIn(outside_points[0], -dirz) - 50.0 )< kCarTolerance );
  
    // For points on surface we get testDistanceToIn = 9e+99, is it ok ?
  }  
  else if ( testCase == 1 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0], -diry) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -diry) - 340.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -diry) - 560.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[5], -diry) - 130.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[6], -diry) - 810.0)< kCarTolerance );
  
    // Add points on surface
  }  
  else if ( testCase == 3 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -diry) - 560.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1],  diry) -   0.0 )< kCarTolerance );
    
   
    // Add points on surface
  }  
  else if ( testCase == 2 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[1],  -diry) - 360.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[4],  -diry) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0],  -dirx) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0],   diry) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0],  -diry) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1],  -dirz) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2],   diry) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2],   dirz) -   0.0 )< kCarTolerance );
   
    // Add points on surface
  }  

  else if ( testCase == 4 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -dirx) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3], -diry) - 100.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0], -dirx) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1], -dirz) -   0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2], -dirz) -   0.0 )< kCarTolerance );
   
   
    // Add points on surface
  }  

  else if ( testCase == 5 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirx) - 66.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirx) - 156.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -dirz) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirz) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0], -dirx) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1],  dirz) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2], -dirz) -  0.0 )< kCarTolerance );
   
   
    // Add points on surface
  } 
 
  else if ( testCase == 6 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirx) - 66.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirx) - 156.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -dirz) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirz) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0], -dirx) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1],  dirz) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2], -dirz) -  0.0 )< kCarTolerance );
   
   
    // Add points on surface
  }  

   else if ( testCase == 7 ) {
    assert( std::fabs(solid->DistanceToIn(outside_points[0],  dirx) - 100.37037037037037 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[1], -dirx) - 159.94397759103646 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[2], -dirz) - 216.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(outside_points[3],  dirz) -  50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[0], -dirx) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[1],  dirz) -  0.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToIn(surface_points[2], -dirz) -  0.0 )< kCarTolerance );
   
   
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
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 50.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 1 ) {
   
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 157.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 142.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 190.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 1465.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 57.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 242.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 1425.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 3 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 287.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 12.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 16.666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 250.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 475.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 107.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 192.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 256.66666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 725.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) -  750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) -  157.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) -  142.5)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 190.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) -  10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) -  1465.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 10.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  else if ( testCase == 2 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 112.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 37.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 50.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 150.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 187.49999999999997158 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) -  750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) -  87.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 102.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 30.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 75.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 550.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) -  2.4999999999999 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) -  107.5 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 3.3333333333333333 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 525.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) -  950.0 )< kCarTolerance );
   

    // For Add points on surface
  }   else if ( testCase == 4 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 550.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 350.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 450.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 700.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 800.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) - 600.0)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 600.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 300.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 1490.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
    else if ( testCase == 5 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 496.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 296.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 396.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 396.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 700.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 800.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 200.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) - 500.6666666666667)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 500.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 200.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 1490.0 )< kCarTolerance );
  
    // For Add points on surface
  }  

   else if ( testCase == 6 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 275.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 496.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 296.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 396.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 296.9467787114846 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 700.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 800.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 200.6666666666667)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) - 500.6666666666667)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 200.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 146.1070975918848 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) - 10.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 1490.0 )< kCarTolerance );
  
    // For Add points on surface
  }  
  
     else if ( testCase == 7 ) {
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirx) - 400.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirx) - 361.1111111111111 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  diry) - 225.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -diry) - 361.1111111111111 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0],  dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[0], -dirz) - 750.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirx) - 496.6666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirx) - 275.92592592592592 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  diry) - 162.94117647058823 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -diry) - 310.41666666666667 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1],  dirz) - 437.97272126752415 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[1], -dirz) - 800.0 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirx) - 196.39821029082771)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirx) - 304.32023010546504)< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  diry) - 301.71673003802277 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -diry) - 198.98689677213173 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2],  dirz) -   9.99999999999999 )< kCarTolerance );
    assert( std::fabs(solid->DistanceToOut(inside_points[2], -dirz) - 1490.0 )< kCarTolerance );
  
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
    assert( std::fabs(solid->GetSurfaceArea() - 1779486.9168738 ) < 1e-6 );
  }
  if ( testCase == 1 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 2239664.8326421054) < 1e-6 );
  }
  if ( testCase == 2 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 1779486.9168738 ) < 1e-6 );
  }
  if ( testCase == 3 ) {
    assert( std::fabs(solid->GetSurfaceArea() - 2239664.8326421054 ) < 1e-6 );
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
    assert( std::fabs(solid->GetCubicVolume() - 33138720.00020987168 ) < 1.0e+6 );
  }
  if ( testCase == 1 ) {
    assert( std::fabs(solid->GetCubicVolume() - 8408880.00005325303 ) < 1.0e+6 );
  }
  if ( testCase == 2 ) {
    assert( std::fabs(solid->GetCubicVolume() - 33191640.00021020696) < 1.0e+6 );
  }
  
  if ( testCase == 3 ) {
    assert( std::fabs(solid->GetCubicVolume() - 8743680.00005537 ) < 1.0e+6 );
  }

  delete solid;
}  
    

//_____________________________________________________________________________
int main()
{
  // Uncomment this line to print the results for a tested solid case
  //   printResults(0);

  for  ( G4int testCase = 0; testCase < 8; ++testCase ) { 
    //G4cout<<"Test ="<<testCase<<G4endl;
    testInside(testCase);
    testDistanceToInPV(testCase);
    testDistanceToOutPV(testCase);
    testSurface(testCase);
    testVolume(testCase);
  }
  
  return 0;
}
