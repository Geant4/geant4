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
// $Id: G4UnionSolid.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// Implementation of methods for the class G4UnionSolid
//
// History:
//
// 12.09.98 V.Grichine: first implementation
// 28.11.98 V.Grichine: fix while loops in DistToIn/Out 
// 27.07.99 V.Grichine: modifications in DistToOut(p,v,...), while -> do-while
// 16.03.01 V.Grichine: modifications in CalculateExtent()
// 17.03.17 E.Tcherniaev: revision of SurfaceNormal()
//
// --------------------------------------------------------------------

#include <sstream>

#include "G4UnionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4VPVParameterisation.hh"
#include "G4GeometryTolerance.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"

///////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid

G4UnionSolid:: G4UnionSolid( const G4String& pName,
                                   G4VSolid* pSolidA ,
                                   G4VSolid* pSolidB )
  : G4BooleanSolid(pName,pSolidA,pSolidB)
{
}

/////////////////////////////////////////////////////////////////////
//
// Constructor
 
G4UnionSolid::G4UnionSolid( const G4String& pName,
                                  G4VSolid* pSolidA ,
                                  G4VSolid* pSolidB ,
                                  G4RotationMatrix* rotMatrix,
                            const G4ThreeVector& transVector )
  : G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)

{
}

///////////////////////////////////////////////////////////
//
// Constructor
 
G4UnionSolid::G4UnionSolid( const G4String& pName,
                                  G4VSolid* pSolidA ,
                                  G4VSolid* pSolidB ,
                            const G4Transform3D& transform )
  : G4BooleanSolid(pName,pSolidA,pSolidB,transform)
{
} 

//////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4UnionSolid::G4UnionSolid( __void__& a )
  : G4BooleanSolid(a)
{
}

///////////////////////////////////////////////////////////
//
// Destructor

G4UnionSolid::~G4UnionSolid()
{
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4UnionSolid::G4UnionSolid(const G4UnionSolid& rhs)
  : G4BooleanSolid (rhs)
{
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4UnionSolid& G4UnionSolid::operator = (const G4UnionSolid& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4BooleanSolid::operator=(rhs);

  return *this;
}  

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UnionSolid::BoundingLimits(G4ThreeVector& pMin,
                                  G4ThreeVector& pMax) const
{
  G4ThreeVector minA,maxA, minB,maxB;
  fPtrSolidA->BoundingLimits(minA,maxA);
  fPtrSolidB->BoundingLimits(minB,maxB);

  pMin.set(std::min(minA.x(),minB.x()),
           std::min(minA.y(),minB.y()),
           std::min(minA.z(),minB.z()));

  pMax.set(std::max(maxA.x(),maxB.x()),
           std::max(maxA.y(),maxB.y()),
           std::max(maxA.z(),maxB.z()));

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UnionSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
     
G4bool 
G4UnionSolid::CalculateExtent( const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pMin,
                                     G4double& pMax ) const 
{
  G4bool   touchesA, touchesB, out ;
  G4double minA =  kInfinity, minB =  kInfinity, 
           maxA = -kInfinity, maxB = -kInfinity; 

  touchesA = fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, 
                                          pTransform, minA, maxA);
  touchesB = fPtrSolidB->CalculateExtent( pAxis, pVoxelLimit, 
                                          pTransform, minB, maxB);
  if( touchesA || touchesB )
  {
    pMin = std::min( minA, minB ); 
    pMax = std::max( maxA, maxB );
    out  = true ; 
  }
  else
  {
    out = false ;
  }

  return out ;  // It exists in this slice if either one does.
}
 
/////////////////////////////////////////////////////
//
// Important comment: When solids A and B touch together along flat
// surface the surface points will be considered as kSurface, while points 
// located around will correspond to kInside

EInside G4UnionSolid::Inside( const G4ThreeVector& p ) const
{
  EInside positionA = fPtrSolidA->Inside(p);
  if (positionA == kInside)  { return kInside; }

  static const G4double rtol
    = 1000*G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  EInside positionB = fPtrSolidB->Inside(p);

  if( positionB == kInside  ||
    ( positionA == kSurface && positionB == kSurface &&
        ( fPtrSolidA->SurfaceNormal(p) + 
          fPtrSolidB->SurfaceNormal(p) ).mag2() < rtol ) )
  {
    return kInside;
  }
  else
  {
    if( ( positionB == kSurface ) || ( positionA == kSurface ) )
      { return kSurface; }
    else
      { return kOutside; } 
  }
}

//////////////////////////////////////////////////////////////
//
// Get surface normal

G4ThreeVector 
G4UnionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  EInside positionA = fPtrSolidA->Inside(p);
  EInside positionB = fPtrSolidB->Inside(p);

  if (positionA == kSurface &&
      positionB == kOutside) return fPtrSolidA->SurfaceNormal(p);

  if (positionA == kOutside &&
      positionB == kSurface) return fPtrSolidB->SurfaceNormal(p);

  if (positionA == kSurface &&
      positionB == kSurface)
  {
    if (Inside(p) == kSurface)
    {
      G4ThreeVector normalA = fPtrSolidA->SurfaceNormal(p);
      G4ThreeVector normalB = fPtrSolidB->SurfaceNormal(p);
      return (normalA + normalB).unit(); 
    }
  }
#ifdef G4BOOLDEBUG
  G4String surf[3] = { "OUTSIDE", "SURFACE", "INSIDE" };
  std::ostringstream message;
  G4int oldprc = message.precision(16);
  message << "Invalid call of SurfaceNormal(p) for union solid: "
          << GetName() << " !"
          << "\nPoint p" << p << " is " << surf[Inside(p)] << " !!!"; 
  message.precision(oldprc);
  G4Exception("G4UnionSolid::SurfaceNormal()", "GeomMgt0001",
              JustWarning, message);
#endif
  return fPtrSolidA->SurfaceNormal(p);
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4UnionSolid::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{
#ifdef G4BOOLDEBUG
  if( Inside(p) == kInside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cout << "          v = " << v << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
  }
#endif

  return std::min(fPtrSolidA->DistanceToIn(p,v),
                    fPtrSolidB->DistanceToIn(p,v) ) ;
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the union of
// two solids

G4double 
G4UnionSolid::DistanceToIn( const G4ThreeVector& p) const 
{
#ifdef G4BOOLDEBUG
  if( Inside(p) == kInside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
  }
#endif
  G4double distA = fPtrSolidA->DistanceToIn(p) ;
  G4double distB = fPtrSolidB->DistanceToIn(p) ;
  G4double safety = std::min(distA,distB) ;
  if(safety < 0.0) safety = 0.0 ;
  return safety ;
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4UnionSolid::DistanceToOut( const G4ThreeVector& p,
           const G4ThreeVector& v,
           const G4bool calcNorm,
                 G4bool *validNorm,
                 G4ThreeVector *n      ) const 
{
  G4double  dist = 0.0, disTmp = 0.0 ;
  G4ThreeVector normTmp;
  G4ThreeVector* nTmp= &normTmp;

  if( Inside(p) == kOutside )
  {
#ifdef G4BOOLDEBUG
      G4cout << "Position:"  << G4endl << G4endl;
      G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
      G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
      G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
      G4cout << "Direction:" << G4endl << G4endl;
      G4cout << "v.x() = "   << v.x() << G4endl;
      G4cout << "v.y() = "   << v.y() << G4endl;
      G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
      G4cout << "WARNING - Invalid call in "
             << "G4UnionSolid::DistanceToOut(p,v)" << G4endl
             << "  Point p is outside !" << G4endl;
      G4cout << "          p = " << p << G4endl;
      G4cout << "          v = " << v << G4endl;
      G4cerr << "WARNING - Invalid call in "
             << "G4UnionSolid::DistanceToOut(p,v)" << G4endl
             << "  Point p is outside !" << G4endl;
      G4cerr << "          p = " << p << G4endl;
      G4cerr << "          v = " << v << G4endl;
#endif
  }
  else
  {
    EInside positionA = fPtrSolidA->Inside(p) ;
    // EInside positionB = fPtrSolidB->Inside(p) ;

    if( positionA != kOutside )
    { 
      do  // Loop checking, 13.08.2015, G.Cosmo
      {
        disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp);
        dist += disTmp ;

        if(fPtrSolidB->Inside(p+dist*v) != kOutside)
        { 
          disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                             validNorm,nTmp);
          dist += disTmp ;
        }
      }
      while( (fPtrSolidA->Inside(p+dist*v) != kOutside)
          && (disTmp > 0.5*kCarTolerance) );
    }
    else // if( positionB != kOutside )
    {
      do  // Loop checking, 13.08.2015, G.Cosmo
      {
        disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp); 
        dist += disTmp ;

        if(fPtrSolidA->Inside(p+dist*v) != kOutside)
        { 
          disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                             validNorm,nTmp);
          dist += disTmp ;
        }
      }
      while( (fPtrSolidB->Inside(p+dist*v) != kOutside)
          && (disTmp > 0.5*kCarTolerance) );
    }
  }
  if( calcNorm )
  { 
     *validNorm = false ;
     *n         = *nTmp ;   
  }
  return dist ;
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4UnionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4double distout = 0.0;
  if( Inside(p) == kOutside )
  {
#ifdef G4BOOLDEBUG
    G4cout << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4UnionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
#endif
  }
  else
  {
    EInside positionA = fPtrSolidA->Inside(p) ;
    EInside positionB = fPtrSolidB->Inside(p) ;
  
    //  Is this equivalent ??
    //    if( ! (  (positionA == kOutside)) && 
    //             (positionB == kOutside))  ) 
    if((positionA == kInside  && positionB == kInside  ) ||
       (positionA == kInside  && positionB == kSurface ) ||
       (positionA == kSurface && positionB == kInside  )     )
    {     
      distout= std::max(fPtrSolidA->DistanceToOut(p),
                          fPtrSolidB->DistanceToOut(p) ) ;
    }
    else
    {
      if(positionA == kOutside)
      {
        distout= fPtrSolidB->DistanceToOut(p) ;
      }
      else
      {
        distout= fPtrSolidA->DistanceToOut(p) ;
      }
    }
  }
  return distout;
}

//////////////////////////////////////////////////////////////
//
//

G4GeometryType G4UnionSolid::GetEntityType() const 
{
  return G4String("G4UnionSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UnionSolid::Clone() const
{
  return new G4UnionSolid(*this);
}

//////////////////////////////////////////////////////////////
//
//

void 
G4UnionSolid::ComputeDimensions(       G4VPVParameterisation*,
                                 const G4int,
                                 const G4VPhysicalVolume* ) 
{
}

/////////////////////////////////////////////////
//
//                    

void 
G4UnionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4UnionSolid::CreatePolyhedron () const 
{
  HepPolyhedronProcessor processor;
  // Stack components and components of components recursively
  // See G4BooleanSolid::StackPolyhedron
  G4Polyhedron* top = StackPolyhedron(processor, this);
  G4Polyhedron* result = new G4Polyhedron(*top);
  if (processor.execute(*result)) { return result; }
  else { return 0; }
}
