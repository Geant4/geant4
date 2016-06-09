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
// $Id: G4GeomTestSegment.cc,v 1.4 2004/09/02 14:44:11 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// --------------------------------------------------------------------
// GEANT 4 class source file
//
// G4GeomTestSegment
//
// Author: D.C.Williams, UCSC (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4GeomTestSegment.hh"

#include "G4VSolid.hh"
#include "G4GeomTestLogger.hh"

//
// Constructor
//
G4GeomTestSegment::G4GeomTestSegment( const G4VSolid *theSolid,
                                      const G4ThreeVector &theP,
                                      const G4ThreeVector &theV,
                                            G4GeomTestLogger *logger )
  : solid(theSolid),
    p0(theP),
    v(theV)
{
  FindPoints(logger);
}


//
// Simple accessors
//
const G4VSolid * G4GeomTestSegment::GetSolid() const { return solid; }

const G4ThreeVector &G4GeomTestSegment::GetP() const { return p0; }

const G4ThreeVector &G4GeomTestSegment::GetV() const { return v; }


//
// Return number points
//
G4int G4GeomTestSegment::GetNumberPoints() const
{
  return points.size();
}


//
// Return ith point
//
const G4GeomTestPoint &G4GeomTestSegment::GetPoint( G4int i ) const
{
  return points[i];
}


//
// FindPoints
//
// Do the dirty work
//
void G4GeomTestSegment::FindPoints( G4GeomTestLogger *logger )
{
  FindSomePoints( logger, true );
  FindSomePoints( logger, false );
  
  PatchInconsistencies( logger );
}


//
// PatchInconsistencies
//
// Search for inconsistancies in the hit list and patch
// them up, if possible
//
void G4GeomTestSegment::PatchInconsistencies(  G4GeomTestLogger *logger )
{
  if (points.size() == 0) return;
  
  //
  // Sort
  //
  std::sort( points.begin(), points.end() );
  
  //
  // Loop over entering/leaving pairs
  //
  std::vector<G4GeomTestPoint>::iterator curr = points.begin();
  do {
    std::vector<G4GeomTestPoint>::iterator next = curr + 1;
  
    //
    // Is the next point close by?
    //
    while (next != points.end() && 
           next->GetDistance()-curr->GetDistance() < kCarTolerance) {
      //
      // Yup. If we have an entering/leaving pair, all is okay
      //
      if (next->Entering() != curr->Entering()) {
        curr = next + 1;
        next = curr + 1;
        break;
      }
      
      //
      // Otherwise, they are duplicate, and just delete one
      //
      next = points.erase(next);
      curr = next - 1;
    }
    
    if (curr == points.end()) break;
  
    //
    // The next point should be entering...
    //
    if (!curr->Entering()) {
      //
      // Point is out of sequence:
      // Test solid just before this point
      //
      G4double s = curr->GetDistance();
      G4ThreeVector p = p0 + s*v;
      G4ThreeVector p1 = p - 10*kCarTolerance*v;
      
      if (solid->Inside(p1) == kOutside) {
        //
        // We are missing an entrance point near the current
        // point. Add one.
        //
        curr = points.insert(curr, G4GeomTestPoint( p, s, true ) );
        ++curr;
        break;
      }
      
      //
      // Test solid just after last point
      //
      if (curr != points.begin()) {
        std::vector<G4GeomTestPoint>::iterator prev = curr - 1;

        s = prev->GetDistance();
        p = p0 + s*v;
        p1 = p + 10*kCarTolerance*v;

        if (solid->Inside(p1) == kOutside) {
          //
          // We are missing an entrance point near the previous
          // point. Add one.
          //
          curr = points.insert(curr, G4GeomTestPoint( p, s, true ) );
          ++curr;
          break;
        }
      }
      
      //
      // Oh oh. No solution. Delete the current point and complain.
      //
      logger->SolidProblem( solid, "Spurious exiting intersection point", p );
      curr = points.erase(curr);
      break;
    }  

    //
    // The following point should be leaving
    //
    if (next == points.end() || next->Entering() ) {
      //
      // But is missing. Check immediately after this point.
      //
      G4double s = curr->GetDistance();
      G4ThreeVector p = p0 + s*v;
      G4ThreeVector p1 = p + 10*kCarTolerance*v;
      
      if (solid->Inside(p1) == kOutside) {
        //
        // We are missing an exit point near the current
        // point. Add one.
        //
        curr = points.insert(next, G4GeomTestPoint( p, s, false ) );
        break;
      }
      
      if (next != points.end()) {
        //
        // Check just before next point
        //
        s = next->GetDistance();
        p = p0 + s*v;
        p1 = p - 10*kCarTolerance*v;
        if (solid->Inside(p1) == kOutside) {
          //
          // We are missing an exit point before the next
          // point. Add one.
          //
          curr = points.insert(next, G4GeomTestPoint( p, s, false ) );
          break;
        }
      }        
      
      //
      // Oh oh. Hanging entering point. Delete and complain.
      //
      logger->SolidProblem( solid, "Spurious entering intersection point", p );
      curr = points.erase(curr);
    }
    
    curr = next + 1;
  } while( curr != points.end() );
  
  //
  // Double check
  //
  if (points.size()&1) 
    logger->SolidProblem( solid,
                          "Solid has odd number of intersection points", p0 );
  
  return;
}  


//
// Find intersections either in front or behind the point
//
void G4GeomTestSegment::FindSomePoints( G4GeomTestLogger *logger,
                                        G4bool forward )
{
  G4double sign = forward ? +1 : -1;

  G4ThreeVector p(p0);
  G4ThreeVector vSearch(sign*v);
  G4double s(0);
  G4bool entering;
  
  //
  // Look for nearest intersection point in the specified
  // direction and return if there isn't one
  //
  G4double dist;
  switch(solid->Inside(p)) {
    case kInside:
      dist = solid->DistanceToOut(p,vSearch);
      if (dist >= kInfinity) {
        logger->SolidProblem( solid,
                "DistanceToOut(p,v) = kInfinity for point inside", p );
        return;
      }
      
      s += sign*dist;
      entering = false;
      break;
    case kOutside:
      dist = solid->DistanceToIn(p,vSearch);
      if (dist >= kInfinity) return;
      
      s += sign*dist;
      entering = true;
      break;
    case kSurface:
      entering = (v.dot(solid->SurfaceNormal(p)) < 0);
      break;
    default:
      logger->SolidProblem( solid,
              "Inside returns illegal enumerated value", p );
      return;
  }
  
  //
  // Loop
  //
  // nzero = the number of consecutive zeros
  //
  G4int nzero=0;
  
  for(;;) {
    //
    // Locate point
    //
    p = p0 + s*v;
    
    if (nzero > 2) {
      //
      // Oops. In a loop. Probably along a spherical or cylindrical surface.
      // Let's give the tool a little help with a push
      //
      G4double push = 1E-6;
      s += sign*push;
      for(;;) {
        p = p0 + s*v;
        EInside inside = solid->Inside(p);
        if (inside == kInside) {
          entering = true;
          break;
        }
        else if (inside == kOutside) {
          entering = false;
          break;
        }

        push = 2*push;
        if (push > 0.1) {
          logger->SolidProblem( solid,
                  "Push fails to fix geometry inconsistency", p );
          return;
        }
        s += sign*push;
      }
    }
    else {

      //
      // Record point
      //
      points.push_back( G4GeomTestPoint( p, s, entering==forward ) );
    }
    
    //
    // Find next intersection
    //
    G4double dist;
    if (entering) {
      dist = solid->DistanceToOut(p,vSearch);
      if (dist >= kInfinity) {
        logger->SolidProblem( solid,
                "DistanceToOut(p,v) = kInfinity for point inside", p );
        return;
      }
      
      if ( (dist > kCarTolerance)
        && (solid->Inside(p + (dist*0.999)*vSearch) == kOutside) ) {
        //
        // This shouldn't have happened
        //
        if (solid->Inside(p) == kOutside)
          logger->SolidProblem(solid,
                  "Entering point is outside (possible roundoff error)",p);
        else
          logger->SolidProblem(solid,
                  "DistanceToOut(p,v) brings trajectory well outside solid",p);
        return;
      }
      
      entering = false;
    }
    else {
      dist = solid->DistanceToIn(p,vSearch);
      if (dist >= kInfinity) return;
      
      if ( (dist > kCarTolerance)
        && (solid->Inside(p + (dist*0.999)*vSearch) == kInside) ) {
        //
        // This shouldn't have happened
        //
        if (solid->Inside(p) == kInside)
          logger->SolidProblem(solid,
                  "Exiting point is inside (possible roundoff error)", p);
        else
          logger->SolidProblem(solid,
                  "DistanceToIn(p,v) brings trajectory well inside solid", p);
        return;
      }
      
      entering = true;
    }
    
    //
    // Update s
    //
    if (dist <= 0) {
      nzero++; 
    }
    else {
      nzero=0;
      s += sign*dist;
    }
  }
}
