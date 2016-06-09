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
// $Id: G4GeomTestVolume.cc,v 1.6 2007-11-16 09:39:14 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestVolume
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestVolume.hh"

#include "G4GeomTestLogger.hh"
#include "G4GeomTestVolPoint.hh"
#include "G4GeomTestSegment.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4VSolid.hh"

#include <vector>
#include <set>
#include <algorithm>
#include <iomanip>

//
// Constructor
//
G4GeomTestVolume::G4GeomTestVolume( const G4VPhysicalVolume *theTarget,
                                          G4GeomTestLogger *theLogger,
                                          G4double theTolerance )
  : target(theTarget),
    logger(theLogger),
    tolerance(theTolerance),
    extent(theTarget->GetLogicalVolume()->GetSolid()->GetExtent())
{;}

//
// Destructor
//
G4GeomTestVolume::~G4GeomTestVolume() {;}

//
// Get error tolerance
//
G4double G4GeomTestVolume::GetTolerance() const
{
  return tolerance;
}

//
// Set error tolerance
//
void G4GeomTestVolume::SetTolerance(G4double val)
{
  tolerance = val;
}

//
// TestCartGridXYZ
//
void G4GeomTestVolume::TestCartGridXYZ( G4int nx, G4int ny, G4int nz )
{
  TestCartGridX( ny, nz );
  TestCartGridY( nz, nx );
  TestCartGridZ( nx, ny );
}

//
// TestCartGridX
//
void G4GeomTestVolume::TestCartGridX( G4int ny, G4int nz )
{
  TestCartGrid( G4ThreeVector(0,1,0), G4ThreeVector(0,0,1),
                G4ThreeVector(1,0,0), ny, nz );
}

//
// TestCartGridY
//
void G4GeomTestVolume::TestCartGridY( G4int nz, G4int nx )
{
  TestCartGrid( G4ThreeVector(0,0,1), G4ThreeVector(1,0,0),
                G4ThreeVector(0,1,0), nz, nx );
}

//
// TestCartGridZ
//
void G4GeomTestVolume::TestCartGridZ( G4int nx, G4int ny )
{
  TestCartGrid( G4ThreeVector(1,0,0), G4ThreeVector(0,1,0),
                G4ThreeVector(0,0,1), nx, ny );
}

//
// TestRecursiveCartGrid
//
void G4GeomTestVolume::TestRecursiveCartGrid( G4int nx, G4int ny, G4int nz,
                                              G4int slevel, G4int depth )
{
  // If reached requested level of depth (i.e. set to 0), exit.
  // If not depth specified (i.e. set to -1), visit the whole tree.
  // If requested initial level of depth is not zero, visit from beginning
  //
  if (depth == 0) return;
  if (depth != -1) depth--;
  if (slevel != 0) slevel--;

  //
  // As long as we aren't a replica and we reached the requested
  // initial level of depth, test ourselves
  //
  if ( (!target->IsReplicated()) && (slevel==0) )
  {
    TestCartGridXYZ( nx, ny, nz );
    ReportErrors();
  }

  //
  // Loop over unique daughters
  //
  std::set<const G4LogicalVolume *> tested;

  const G4LogicalVolume *logical = target->GetLogicalVolume();
  G4int nDaughter = logical->GetNoDaughters();
  G4int iDaughter;
  for( iDaughter=0; iDaughter<nDaughter; ++iDaughter )
  {
    const G4VPhysicalVolume *daughter =
          logical->GetDaughter(iDaughter);
    const G4LogicalVolume *daughterLogical =
          daughter->GetLogicalVolume();
    
    //
    // Skip empty daughters
    //
    if (daughterLogical->GetNoDaughters() == 0) continue;
    
    //
    // Tested already?
    //
    std::pair<std::set<const G4LogicalVolume *>::iterator,G4bool>
           there = tested.insert(daughterLogical);
    if (!there.second) continue;

    //
    // Recurse
    //
    G4GeomTestVolume vTest( daughter, logger, tolerance );
    vTest.TestRecursiveCartGrid( nx,ny,nz,slevel,depth );
  }
}

//
// TestRecursiveCylinder
//
void
G4GeomTestVolume::TestRecursiveCylinder( G4int nPhi, G4int nZ, G4int nRho,
                                         G4double fracZ, G4double fracRho,
                                         G4bool usePhi,
                                         G4int slevel, G4int depth )
{
  // If reached requested level of depth (i.e. set to 0), exit.
  // If not depth specified (i.e. set to -1), visit the whole tree.
  // If requested initial level of depth is not zero, visit from beginning
  //
  if (depth == 0) return;
  if (depth != -1) depth--;
  if (slevel != 0) slevel--;

  //
  // As long as we aren't a replica and we reached the requested
  // initial level of depth, test ourselves
  //
  if ( (!target->IsReplicated()) && (slevel==0) )
  {
    TestCylinder( nPhi, nZ, nRho, fracZ, fracRho, usePhi );
    ReportErrors();
  }

  //
  // Loop over unique daughters
  //
  std::set<const G4LogicalVolume *> tested;

  const G4LogicalVolume *logical = target->GetLogicalVolume();
  G4int nDaughter = logical->GetNoDaughters();
  G4int iDaughter;
  for( iDaughter=0; iDaughter<nDaughter; ++iDaughter )
  {
    const G4VPhysicalVolume *daughter =
          logical->GetDaughter(iDaughter);
    const G4LogicalVolume *daughterLogical =
          daughter->GetLogicalVolume();
    
    //
    // Skip empty daughters
    //
    if (daughterLogical->GetNoDaughters() == 0) continue;
    
    //
    // Tested already?
    //
    std::pair<std::set<const G4LogicalVolume *>::iterator,G4bool>
           there = tested.insert(daughterLogical);
    if (!there.second) continue;

    //
    // Recurse
    //
    G4GeomTestVolume vTest( daughter, logger, tolerance );
    vTest.TestRecursiveCylinder( nPhi, nZ, nRho,
                                 fracZ, fracRho, usePhi,
                                 slevel, depth );
  }
}

//
// TestCylinder
//
void G4GeomTestVolume::TestCylinder( G4int nPhi, G4int nZ, G4int nRho,
                                     G4double fracZ, G4double fracRho,
                                     G4bool usePhi    )
{
  //
  // Get size of our volume
  //
  G4double xMax = std::max(extent.GetXmax(),-extent.GetXmin());
  G4double yMax = std::max(extent.GetYmax(),-extent.GetYmin());
  G4double rhoMax = std::sqrt(xMax*xMax + yMax*yMax);
  
  G4double zMax = extent.GetZmax();
  G4double zMin = extent.GetZmin();
  G4double zHalfLength = 0.5*(zMax-zMin);
  G4double z0 = 0.5*(zMax+zMin);
  
  //
  // Loop over phi
  //
  G4double deltaPhi = 2*pi/G4double(nPhi);
  
  G4double phi = 0;
  G4int iPhi = nPhi;
  if ((iPhi&1) == 0) iPhi++;  // Also use odd number phi slices
  do
  {
    G4double cosPhi = std::cos(phi);
    G4double sinPhi = std::sin(phi);
    G4ThreeVector vPhi1(sinPhi,-cosPhi,0);

    //
    // Loop over rho
    //
    G4double rho = rhoMax;
    G4int iRho = nRho;
    do
    {
      G4ThreeVector p(rho*cosPhi,rho*sinPhi,0);
      static G4ThreeVector v(0,0,1);
      
      TestOneLine( p, v );
      
      if (usePhi)
      {
        //
        // Loop over z
        //
        G4double zScale = 1.0;
        G4int iZ=nZ;
        do
        {
          p.setZ( z0 + zScale*zHalfLength );
          TestOneLine(p,vPhi1);
          p.setZ( z0 - zScale*zHalfLength );
          TestOneLine(p,vPhi1);
        } while( zScale *= fracZ, --iZ );
      }
    } while( rho *= fracRho, --iRho );

    //
    // Loop over z
    //
    G4ThreeVector p(0,0,0);
    G4ThreeVector vPhi2(cosPhi,sinPhi,0);
    
    G4double zScale = 1.0;
    G4int iZ=nZ;
    do
    {
      p.setZ( z0 + zScale*zHalfLength );
      
      TestOneLine(p,vPhi2);
      
      p.setZ( z0 - zScale*zHalfLength );
      
      TestOneLine(p,vPhi2);
    } while( zScale *= fracZ, --iZ );
    
  } while( phi += deltaPhi, --iPhi );
}

//
// TestCartGrid
//
void G4GeomTestVolume::TestCartGrid( const G4ThreeVector &theG1,
                                     const G4ThreeVector &theG2,
                                     const G4ThreeVector &theV,
                                     G4int n1, G4int n2 )
{
  if (n1 <= 0 || n2 <= 0) 
    G4Exception( "G4GeomTestVolume::TestCartGrid()", "InvalidParameters",
                 FatalException, "Arguments n1 and n2 must be >= 1" );
    
  G4ThreeVector xMin( extent.GetXmin(), extent.GetYmin(),
                      extent.GetZmin() );
  G4ThreeVector xMax( extent.GetXmax(), extent.GetYmax(),
                      extent.GetZmax() );
  
  G4ThreeVector g1(theG1.unit());
  G4ThreeVector g2(theG2.unit());
  G4ThreeVector v(theV.unit());
  
  G4double gMin1 = g1.dot(xMin);
  G4double gMax1 = g1.dot(xMax);
  
  G4double gMin2 = g2.dot(xMin);
  G4double gMax2 = g2.dot(xMax);
  
  G4double delta1 = (gMax1-gMin1)/G4double(n1);
  G4double delta2 = (gMax2-gMin2)/G4double(n2);
    
  G4int i1, i2;
  
  for(i1=0;i1<=n1;++i1)
  {
    G4ThreeVector p1 = (gMin1 + G4double(i1)*delta1)*g1;
    
    for(i2=0;i2<=n2;++i2)
    {
      G4ThreeVector p2 = (gMin2 + G4double(i2)*delta2)*g2;
      
      TestOneLine( p1+p2, v );
    }
  }
}  

//
// TestRecursiveLine
//
void
G4GeomTestVolume::TestRecursiveLine( const G4ThreeVector& p,
                                     const G4ThreeVector& v,
                                     G4int slevel, G4int depth )
{
  // If reached requested level of depth (i.e. set to 0), exit.
  // If not depth specified (i.e. set to -1), visit the whole tree.
  // If requested initial level of depth is not zero, visit from beginning
  //
  if (depth == 0) return;
  if (depth != -1) depth--;
  if (slevel != 0) slevel--;

  //
  // As long as we aren't a replica and we reached the requested
  // initial level of depth, test ourselves
  //
  if ( (!target->IsReplicated()) && (slevel==0) )
  {
    TestOneLine( p, v );
    ReportErrors();
  }

  //
  // Loop over unique daughters
  //
  std::set<const G4LogicalVolume *> tested;

  const G4LogicalVolume *logical = target->GetLogicalVolume();
  G4int nDaughter = logical->GetNoDaughters();
  G4int iDaughter;
  for( iDaughter=0; iDaughter<nDaughter; ++iDaughter )
  {
    const G4VPhysicalVolume *daughter =
          logical->GetDaughter(iDaughter);
    const G4LogicalVolume *daughterLogical =
          daughter->GetLogicalVolume();
    
    //
    // Skip empty daughters
    //
    if (daughterLogical->GetNoDaughters() == 0) continue;
    
    //
    // Tested already?
    //
    std::pair<std::set<const G4LogicalVolume *>::iterator,G4bool>
           there = tested.insert(daughterLogical);
    if (!there.second) continue;

    //
    // Recurse
    //
    G4GeomTestVolume vTest( daughter, logger, tolerance );
    vTest.TestRecursiveLine( p, v, slevel, depth );
  }
}

//
// TestOneLine
//
// Test geometry consistency along a single line
//
void G4GeomTestVolume::TestOneLine( const G4ThreeVector &p,
                                    const G4ThreeVector &v )
{
  //
  // Keep track of intersection points
  //
  std::vector<G4GeomTestVolPoint> points;
  
  //
  // Calculate intersections with the mother volume
  //
  G4GeomTestSegment targetSegment( target->GetLogicalVolume()->GetSolid(),
                                   p, v, logger );
  
  //
  // Add them to our list
  // 
  G4int n = targetSegment.GetNumberPoints();
  G4int i;
  for(i=0;i<n;++i)
  {
    points.push_back( G4GeomTestVolPoint( targetSegment.GetPoint(i), -1 ) );
  } 

  //
  // Loop over daughter volumes
  //
  const G4LogicalVolume *logical = target->GetLogicalVolume();
  G4int nDaughter = logical->GetNoDaughters();
  G4int iDaughter;
  for( iDaughter=0; iDaughter<nDaughter; ++iDaughter)
  {
    const G4VPhysicalVolume *daughter = 
          logical->GetDaughter(iDaughter);
    
    //
    // Convert coordinates to daughter local coordinates
    //
    const G4RotationMatrix *rotation = daughter->GetFrameRotation();
    const G4ThreeVector &translation = daughter->GetFrameTranslation();
    
    G4ThreeVector pLocal = translation + p;
    G4ThreeVector vLocal = v;
    
    if (rotation)
    {
      pLocal = (*rotation)*pLocal;
      vLocal = (*rotation)*vLocal;
    }
    
    //
    // Find intersections
    //
    G4GeomTestSegment
       daughterSegment( daughter->GetLogicalVolume()->GetSolid(), 
                        pLocal, vLocal, logger );
    
    //
    // Add them to the list
    //
    n = daughterSegment.GetNumberPoints();
    for(i=0;i<n;++i)
    {
      points.push_back( G4GeomTestVolPoint( daughterSegment.GetPoint(i),
            iDaughter, translation, rotation ) );
    } 
  }
  
  //
  // Now sort the list of intersection points
  //
  std::sort( points.begin(), points.end() );
  
  //
  // Search for problems:
  //    1. Overlap daughters will be indicated by intersecting
  //       points that lie within another daughter
  //    2. Daughter extending outside the mother will be
  //       indicated by intersecting points outside the mother
  //
  // The search method is always looking forward when
  // resolving discrepencies due to roundoff error. Once
  // one instance of a daughter intersection is analyzed,
  // it is removed from further consideration
  //
  n = points.size();
  
  //
  // Set true if this point has been analyzed
  //
  std::vector<G4bool> checked( n, false );
  
  for(i=0;i<n;++i)
  {
    if (checked[i]) continue;
  
    G4int iDaug = points[i].GetDaughterIndex();
    if (iDaug < 0) continue;
  
    //
    // Intersection found. Find matching pair.
    //
    G4double iS = points[i].GetDistance();
    G4int j = i;
    while(++j<n)
    {
      if (iDaug == points[j].GetDaughterIndex()) break;
    }
    if (j>=n)
    {
      //
      // Unmatched? This shouldn't happen
      //
      logger->SolidProblem( logical->GetDaughter(iDaug)->
                            GetLogicalVolume()->GetSolid(),
                            "Unmatched intersection point",
          points[i].GetPosition() );
      continue;
    }
    
    checked[j] = true;
    
    G4double jS = points[j].GetDistance();

    //
    // Otherwise, we could have a problem
    //
    G4int k = i;
    while(++k<j)
    {
      if (checked[k]) continue;
      
      G4bool kEntering = points[k].Entering();
      G4double kS = points[k].GetDistance();
      
      
      //
      // Problem found: catagorize
      //
      G4int kDaug = points[k].GetDaughterIndex();
      if (kDaug < 0)
      {
        //
        // Ignore small overshoots if they are within tolerance
        //
        if (kS-iS < tolerance &&   kEntering ) continue;
        if (jS-kS < tolerance && (!kEntering)) continue;

        //
        // We appear to extend outside the mother volume
        //
        std::map<G4long,G4GeomTestOvershootList>::iterator overshoot =
          overshoots.find(iDaug);
        if (overshoot == overshoots.end())
        {
          std::pair<std::map<G4long,G4GeomTestOvershootList>::iterator,G4bool>
            result =
              overshoots.insert( std::pair<const G4long,G4GeomTestOvershootList>
                               (iDaug,G4GeomTestOvershootList(target,iDaug)) );
          assert(result.second);
          overshoot = result.first;
        }

        if (kEntering)
          (*overshoot).second.AddError( points[i].GetPosition(),
                                        points[k].GetPosition() );
        else
          (*overshoot).second.AddError( points[k].GetPosition(),
                                        points[j].GetPosition() );
      }
      else
      {
        //
        // Ignore small overlaps if they are within tolerance
        //
        if (kS-iS < tolerance && (!kEntering)) continue;
        if (jS-kS < tolerance &&   kEntering ) continue;
        
        //
        // We appear to overlap with another daughter
        //
        G4long key = iDaug < kDaug ?
               (iDaug*nDaughter + kDaug) : (kDaug*nDaughter + iDaug);
        
        std::map<G4long,G4GeomTestOverlapList>::iterator overlap =
          overlaps.find(key);
        if (overlap == overlaps.end())
        {
          std::pair<std::map<G4long,G4GeomTestOverlapList>::iterator,G4bool>
            result =
            overlaps.insert( std::pair<const G4long,G4GeomTestOverlapList>
                           (key,G4GeomTestOverlapList(target,iDaug,kDaug)) );
          assert(result.second);
          overlap = result.first;
        }

        if (points[k].Entering())
          (*overlap).second.AddError( points[k].GetPosition(),
                                      points[j].GetPosition() );
        else
          (*overlap).second.AddError( points[i].GetPosition(),
                                      points[k].GetPosition() );
      }
    }
  }
}

//
// ReportErrors
//
void G4GeomTestVolume::ReportErrors()
{
  //
  // Report overshoots
  //
  if (overshoots.empty())
    logger->NoProblem("GeomTest: no daughter volume extending outside mother detected.");
  else
  {
    std::map<G4long,G4GeomTestOvershootList>::iterator overshoot =
      overshoots.begin();
    while( overshoot != overshoots.end() )
    {
      logger->OvershootingDaughter( &(*overshoot).second );
      ++overshoot;
    }
  }

  //
  // Report overlaps
  //
  if (overlaps.empty())
    logger->NoProblem("GeomTest: no overlapping daughters detected.");
  else
  {
    std::map<G4long,G4GeomTestOverlapList>::iterator overlap =
      overlaps.begin();
    while( overlap != overlaps.end() )
    {
      logger->OverlappingDaughters( &(*overlap).second );
      ++overlap;
    }
  }
}

//
// ClearErrors
//
void G4GeomTestVolume::ClearErrors()
{
  overshoots.clear();
  overlaps.clear();
}
