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
// Implementation of methods for the class G4IntersectionSolid
//
// 12.09.98 V.Grichine: first implementation
// --------------------------------------------------------------------

#include <sstream>

#include "G4IntersectionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"

/////////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid
//

G4IntersectionSolid::G4IntersectionSolid( const G4String& pName,
                                                G4VSolid* pSolidA ,
                                                G4VSolid* pSolidB   )
  : G4BooleanSolid(pName,pSolidA,pSolidB)
{
} 

///////////////////////////////////////////////////////////////////
//

G4IntersectionSolid::G4IntersectionSolid( const G4String& pName,
                                                G4VSolid* pSolidA,
                                                G4VSolid* pSolidB,
                                                G4RotationMatrix* rotMatrix,
                                          const G4ThreeVector& transVector  )
  : G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)
{
}

//////////////////////////////////////////////////////////////////
//
// 
 
G4IntersectionSolid::G4IntersectionSolid( const G4String& pName,
                                                G4VSolid* pSolidA,
                                                G4VSolid* pSolidB,
                                          const G4Transform3D& transform )
  : G4BooleanSolid(pName,pSolidA,pSolidB,transform)
{
} 

//////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4IntersectionSolid::G4IntersectionSolid( __void__& a )
  : G4BooleanSolid(a)
{
}

///////////////////////////////////////////////////////////////
//
//

G4IntersectionSolid::~G4IntersectionSolid()
{
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4IntersectionSolid::G4IntersectionSolid(const G4IntersectionSolid& rhs)
  : G4BooleanSolid (rhs)
{
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4IntersectionSolid&
G4IntersectionSolid::operator = (const G4IntersectionSolid& rhs) 
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

void
G4IntersectionSolid::BoundingLimits(G4ThreeVector& pMin,
                                    G4ThreeVector& pMax) const
{
  G4ThreeVector minA,maxA, minB,maxB;
  fPtrSolidA->BoundingLimits(minA,maxA);
  fPtrSolidB->BoundingLimits(minB,maxB);

  pMin.set(std::max(minA.x(),minB.x()),
           std::max(minA.y(),minB.y()),
           std::max(minA.z(),minB.z()));

  pMax.set(std::min(maxA.x(),maxB.x()),
           std::min(maxA.y(),maxB.y()),
           std::min(maxA.z(),maxB.z()));

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4IntersectionSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
     
G4bool 
G4IntersectionSolid::CalculateExtent(const EAxis pAxis,
                                     const G4VoxelLimits& pVoxelLimit,
                                     const G4AffineTransform& pTransform,
                                           G4double& pMin,
                                           G4double& pMax) const 
{
  G4bool   retA, retB, out;
  G4double minA, minB, maxA, maxB; 

  retA = fPtrSolidA
          ->CalculateExtent( pAxis, pVoxelLimit, pTransform, minA, maxA);
  retB = fPtrSolidB
          ->CalculateExtent( pAxis, pVoxelLimit, pTransform, minB, maxB);

  if( retA && retB )
  {
    pMin = std::max( minA, minB ); 
    pMax = std::min( maxA, maxB );
    out  = (pMax > pMin); // true;
  }
  else
  {
    out = false;
  }

  return out; // It exists in this slice only if both exist in it.
}
 
/////////////////////////////////////////////////////
//
// Touching ? Empty intersection ?

EInside G4IntersectionSolid::Inside(const G4ThreeVector& p) const
{
  EInside positionA = fPtrSolidA->Inside(p);
  if(positionA == kOutside) return positionA; // outside A

  EInside positionB = fPtrSolidB->Inside(p);
  if(positionA == kInside)  return positionB;

  if(positionB == kOutside) return positionB; // outside B
  return kSurface;                            // surface A & B
}

//////////////////////////////////////////////////////////////
//

G4ThreeVector 
G4IntersectionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector normal;
  EInside insideA, insideB;
  
  insideA = fPtrSolidA->Inside(p);
  insideB = fPtrSolidB->Inside(p);

#ifdef G4BOOLDEBUG
  if( (insideA == kOutside) || (insideB == kOutside) )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
  }
#endif

  // On the surface of both is difficult ... treat it like on A now!
  //
  if( insideA == kSurface )
  {
    normal = fPtrSolidA->SurfaceNormal(p) ;
  }
  else if( insideB == kSurface )
  {
    normal = fPtrSolidB->SurfaceNormal(p) ;
  } 
  else  // We are on neither surface, so we should generate an exception
  {
    if(fPtrSolidA->DistanceToOut(p) <= fPtrSolidB->DistanceToOut(p) )
    { 
      normal= fPtrSolidA->SurfaceNormal(p) ;   
    }
    else
    {
      normal= fPtrSolidB->SurfaceNormal(p) ;   
    }
#ifdef G4BOOLDEBUG
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is out of surface !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is out of surface !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
#endif
    }

  return normal;
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4IntersectionSolid::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{
  G4double dist = 0.0;
  if( Inside(p) == kInside )
  {
#ifdef G4BOOLDEBUG
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cout << "          v = " << v << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
#endif
  }
  else // if( Inside(p) == kSurface ) 
  {
    EInside wA = fPtrSolidA->Inside(p);
    EInside wB = fPtrSolidB->Inside(p);

    G4ThreeVector pA = p,  pB = p;
    G4double      dA = 0., dA1=0., dA2=0.;
    G4double      dB = 0., dB1=0., dB2=0.;
    G4bool        doA = true, doB = true;

    static const size_t max_trials=10000;
    for (size_t trial=0; trial<max_trials; ++trial) 
    {
      if(doA) 
      {
        // find next valid range for A

        dA1 = 0.;

        if( wA != kInside ) 
        {
          dA1 = fPtrSolidA->DistanceToIn(pA, v);

          if( dA1 == kInfinity )   return kInfinity;
        
          pA += dA1*v;
        }
        dA2 = dA1 + fPtrSolidA->DistanceToOut(pA, v);
      }
      dA1 += dA;
      dA2 += dA;

      if(doB) 
      {
        // find next valid range for B

        dB1 = 0.;
        if(wB != kInside) 
        {
          dB1 = fPtrSolidB->DistanceToIn(pB, v);

          if(dB1 == kInfinity)   return kInfinity;
        
          pB += dB1*v;
        }
        dB2 = dB1 + fPtrSolidB->DistanceToOut(pB, v);
      }
      dB1 += dB;
      dB2 += dB;

       // check if they overlap

      if( dA1 < dB1 ) 
      {
        if( dB1 < dA2 )  return dB1;

        dA   = dA2;
        pA   = p + dA*v;  // continue from here
        wA   = kSurface;
        doA  = true;
        doB  = false;
      }
      else 
      {
        if( dA1 < dB2 )  return dA1;

        dB   = dB2;
        pB   = p + dB*v;  // continue from here
        wB   = kSurface;
        doB  = true;
        doA  = false;
      }
    }
  }
#ifdef G4BOOLDEBUG
  G4Exception("G4IntersectionSolid::DistanceToIn(p,v)",
              "GeomSolids0001", JustWarning,
              "Reached maximum number of iterations! Returning zero.");
#endif
  return dist ;  
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids

G4double 
G4IntersectionSolid::DistanceToIn( const G4ThreeVector& p) const 
{
#ifdef G4BOOLDEBUG
  if( Inside(p) == kInside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
  }
#endif
  EInside sideA = fPtrSolidA->Inside(p) ;
  EInside sideB = fPtrSolidB->Inside(p) ;
  G4double dist=0.0 ;

  if( sideA != kInside && sideB != kOutside )
  {
    dist = fPtrSolidA->DistanceToIn(p) ;
  }
  else
  {
    if( sideB != kInside && sideA != kOutside )
    {
      dist = fPtrSolidB->DistanceToIn(p) ;
    }
    else
    {
      dist =  std::min(fPtrSolidA->DistanceToIn(p),
                       fPtrSolidB->DistanceToIn(p) ) ; 
    }
  }
  return dist ;
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4IntersectionSolid::DistanceToOut( const G4ThreeVector& p,
                                    const G4ThreeVector& v,
                                    const G4bool calcNorm,
                                          G4bool *validNorm,
                                          G4ThreeVector *n ) const 
{
  G4bool         validNormA, validNormB;
  G4ThreeVector  nA, nB;

#ifdef G4BOOLDEBUG
  if( Inside(p) == kOutside )
  {
    G4cout << "Position:"  << G4endl << G4endl;
    G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
    G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
    G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
    G4cout << "Direction:" << G4endl << G4endl;
    G4cout << "v.x() = "   << v.x() << G4endl;
    G4cout << "v.y() = "   << v.y() << G4endl;
    G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToOut(p,v)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cout << "          v = " << v << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToOut(p,v)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
  }
#endif
  G4double distA = fPtrSolidA->DistanceToOut(p,v,calcNorm,&validNormA,&nA) ;
  G4double distB = fPtrSolidB->DistanceToOut(p,v,calcNorm,&validNormB,&nB) ;

  G4double dist = std::min(distA,distB) ; 

  if( calcNorm )
  {
    if ( distA < distB )
    {
       *validNorm = validNormA;
       *n =         nA;
    }
    else
    {   
       *validNorm = validNormB;
       *n =         nB;
    }
  }

  return dist ; 
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4IntersectionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
#ifdef G4BOOLDEBUG
  if( Inside(p) == kOutside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4IntersectionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
  }
#endif

  return std::min(fPtrSolidA->DistanceToOut(p),
                  fPtrSolidB->DistanceToOut(p) ) ; 

}

//////////////////////////////////////////////////////////////
//
// ComputeDimensions

void 
G4IntersectionSolid::ComputeDimensions( G4VPVParameterisation*,
                                  const G4int,
                                        const G4VPhysicalVolume* ) 
{
}

/////////////////////////////////////////////////
//
// GetEntityType

G4GeometryType G4IntersectionSolid::GetEntityType() const 
{
  return G4String("G4IntersectionSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4IntersectionSolid::Clone() const
{
  return new G4IntersectionSolid(*this);
}

/////////////////////////////////////////////////
//
// DescribeYourselfTo

void 
G4IntersectionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////
//
// CreatePolyhedron

G4Polyhedron* 
G4IntersectionSolid::CreatePolyhedron () const 
{
  HepPolyhedronProcessor processor;
  // Stack components and components of components recursively
  // See G4BooleanSolid::StackPolyhedron
  G4Polyhedron* top = StackPolyhedron(processor, this);
  G4Polyhedron* result = new G4Polyhedron(*top);
  if (processor.execute(*result)) { return result; }
  else { return nullptr; }
}
