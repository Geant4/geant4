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
// $Id: G4SubtractionSolid.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 14.10.98 V.Grichine: implementation of the first version 
// 19.10.98 V.Grichine: new algorithm of DistanceToIn(p,v)
// 02.08.99 V.Grichine: bugs fixed in DistanceToOut(p,v,...)
//                      while -> do-while & surfaceA limitations
// 13.09.00 V.Grichine: bug fixed in SurfaceNormal(p), p can be inside
// 22.07.11 T.Nikitina: add detection of Infinite Loop in DistanceToIn(p,v)
//
// --------------------------------------------------------------------

#include "G4SubtractionSolid.hh"

#include "G4SystemOfUnits.hh"
#include "G4VoxelLimits.hh"
#include "G4VPVParameterisation.hh"
#include "G4GeometryTolerance.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "HepPolyhedronProcessor.h"

#include <sstream>

///////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid

G4SubtractionSolid::G4SubtractionSolid( const G4String& pName,
                                              G4VSolid* pSolidA ,
                                              G4VSolid* pSolidB   )
  : G4BooleanSolid(pName,pSolidA,pSolidB)
{
}

///////////////////////////////////////////////////////////////
//
// Constructor
 
G4SubtractionSolid::G4SubtractionSolid( const G4String& pName,
                                              G4VSolid* pSolidA ,
                                              G4VSolid* pSolidB ,
                                              G4RotationMatrix* rotMatrix,
                                        const G4ThreeVector& transVector )
  : G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)
{
} 

///////////////////////////////////////////////////////////////
//
// Constructor

G4SubtractionSolid::G4SubtractionSolid( const G4String& pName,
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

G4SubtractionSolid::G4SubtractionSolid( __void__& a )
  : G4BooleanSolid(a)
{
}

///////////////////////////////////////////////////////////////
//
// Destructor

G4SubtractionSolid::~G4SubtractionSolid()
{
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4SubtractionSolid::G4SubtractionSolid(const G4SubtractionSolid& rhs)
  : G4BooleanSolid (rhs)
{
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4SubtractionSolid&
G4SubtractionSolid::operator = (const G4SubtractionSolid& rhs) 
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
G4SubtractionSolid::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  // Since it is unclear how the shape of the first solid will be changed
  // after subtraction, just return its original bounding box.
  //
  fPtrSolidA->BoundingLimits(pMin,pMax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4SubtractionSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
     
G4bool 
G4SubtractionSolid::CalculateExtent( const EAxis pAxis,
                                     const G4VoxelLimits& pVoxelLimit,
                                     const G4AffineTransform& pTransform,
                                           G4double& pMin,
                                           G4double& pMax ) const 
{
  // Since we cannot be sure how much the second solid subtracts 
  // from the first, we must use the first solid's extent!

  return fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, 
                                      pTransform, pMin, pMax );
}
 
/////////////////////////////////////////////////////
//
// Touching ? Empty subtraction ?

EInside G4SubtractionSolid::Inside( const G4ThreeVector& p ) const
{
  EInside positionA = fPtrSolidA->Inside(p);
  if (positionA == kOutside) return kOutside;

  EInside positionB = fPtrSolidB->Inside(p);
  
  if(positionA == kInside && positionB == kOutside)
  {
    return kInside ;
  }
  else
  {
    static const G4double rtol
      = 1000.0*G4GeometryTolerance::GetInstance()->GetRadialTolerance();
    if(( positionA == kInside && positionB == kSurface) ||
       ( positionB == kOutside && positionA == kSurface) ||
       ( positionA == kSurface && positionB == kSurface &&
         ( fPtrSolidA->SurfaceNormal(p) - 
           fPtrSolidB->SurfaceNormal(p) ).mag2() > rtol ) )
    {
      return kSurface;
    }
    else
    {
      return kOutside;
    }
  }
}

//////////////////////////////////////////////////////////////
//
// SurfaceNormal

G4ThreeVector 
G4SubtractionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector normal;

  EInside InsideA = fPtrSolidA->Inside(p);
  EInside InsideB = fPtrSolidB->Inside(p); 

  if( InsideA == kOutside )
  {
#ifdef G4BOOLDEBUG
    G4cout << "WARNING - Invalid call [1] in "
           << "G4SubtractionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call [1] in "
           << "G4SubtractionSolid::SurfaceNormal(p)" << G4endl
           << "  Point p is outside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
#endif
    normal = fPtrSolidA->SurfaceNormal(p) ;
  }
  else if( InsideA == kSurface && 
           InsideB != kInside      ) 
  {
    normal = fPtrSolidA->SurfaceNormal(p) ;
  }
  else if( InsideA == kInside && 
           InsideB != kOutside    )
  {
    normal = -fPtrSolidB->SurfaceNormal(p) ;
  }
  else 
  {
    if ( fPtrSolidA->DistanceToOut(p) <= fPtrSolidB->DistanceToIn(p) )
    {
      normal = fPtrSolidA->SurfaceNormal(p) ;
    }
    else
    {
      normal = -fPtrSolidB->SurfaceNormal(p) ;
    }
#ifdef G4BOOLDEBUG
    if(Inside(p) == kInside)
    {
      G4cout << "WARNING - Invalid call [2] in "
             << "G4SubtractionSolid::SurfaceNormal(p)" << G4endl
             << "  Point p is inside !" << G4endl;
      G4cout << "          p = " << p << G4endl;
      G4cerr << "WARNING - Invalid call [2] in "
             << "G4SubtractionSolid::SurfaceNormal(p)" << G4endl
             << "  Point p is inside !" << G4endl;
      G4cerr << "          p = " << p << G4endl;
    }
#endif
  }
  return normal;
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4SubtractionSolid::DistanceToIn(  const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{
  G4double dist = 0.0, dist2 = 0.0, disTmp = 0.0;
    
#ifdef G4BOOLDEBUG
  if( Inside(p) == kInside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cout << "          v = " << v << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToIn(p,v)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
  }
#endif

    // if( // ( fPtrSolidA->Inside(p) != kOutside) &&  // case1:p in both A&B 
    if ( fPtrSolidB->Inside(p) != kOutside )   // start: out of B
    {
      dist = fPtrSolidB->DistanceToOut(p,v) ; // ,calcNorm,validNorm,n) ;
      
      if( fPtrSolidA->Inside(p+dist*v) != kInside )
      {
        G4int count1=0;
        do   // Loop checking, 13.08.2015, G.Cosmo
        {
          disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;

          if(disTmp == kInfinity)
          {  
            return kInfinity ;
          }
          dist += disTmp ;

          if( Inside(p+dist*v) == kOutside )
          {
            disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v) ;
            dist2 = dist+disTmp;
            if (dist == dist2)  { return dist; }   // no progress
            dist = dist2 ;
            count1++;
            if( count1 > 1000 )  // Infinite loop detected
            {
              G4String nameB = fPtrSolidB->GetName();
              if(fPtrSolidB->GetEntityType()=="G4DisplacedSolid")
              {
                nameB = (dynamic_cast<G4DisplacedSolid*>(fPtrSolidB))
                        ->GetConstituentMovedSolid()->GetName();
              }
              std::ostringstream message;
              message << "Illegal condition caused by solids: "
                      << fPtrSolidA->GetName() << " and " << nameB << G4endl;
              message.precision(16);
              message << "Looping detected in point " << p+dist*v
                      << ", from original point " << p
                      << " and direction " << v << G4endl
                      << "Computed candidate distance: " << dist << "*mm. ";
              message.precision(6);
              DumpInfo();
              G4Exception("G4SubtractionSolid::DistanceToIn(p,v)",
                          "GeomSolids1001", JustWarning, message,
                          "Returning candidate distance.");
              return dist;
            }
          }    
        }
        while( Inside(p+dist*v) == kOutside ) ;
      }
    }
    else // p outside A, start in A
    {
      dist = fPtrSolidA->DistanceToIn(p,v) ;

      if( dist == kInfinity ) // past A, hence past A\B
      {
        return kInfinity ;
      }
      else
      {
        G4int count2=0;
        while( Inside(p+dist*v) == kOutside )  // pushing loop
        {
          disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v) ;
          dist += disTmp ;

          if( Inside(p+dist*v) == kOutside )
          { 
            disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;

            if(disTmp == kInfinity) // past A, hence past A\B
            {  
              return kInfinity ;
            }                 
            dist2 = dist+disTmp;
            if (dist == dist2)  { return dist; }   // no progress
            dist = dist2 ;
            count2++;
            if( count2 > 1000 )  // Infinite loop detected
            {
              G4String nameB = fPtrSolidB->GetName();
              if(fPtrSolidB->GetEntityType()=="G4DisplacedSolid")
              {
                nameB = (dynamic_cast<G4DisplacedSolid*>(fPtrSolidB))
                        ->GetConstituentMovedSolid()->GetName();
              }
              std::ostringstream message;
              message << "Illegal condition caused by solids: "
                      << fPtrSolidA->GetName() << " and " << nameB << G4endl;
              message.precision(16);
              message << "Looping detected in point " << p+dist*v
                      << ", from original point " << p
                      << " and direction " << v << G4endl
                      << "Computed candidate distance: " << dist << "*mm. ";
              message.precision(6);
              DumpInfo();
              G4Exception("G4SubtractionSolid::DistanceToIn(p,v)",
                          "GeomSolids1001", JustWarning, message,
                          "Returning candidate distance.");
              return dist;
            }
          }
        }    // Loop checking, 13.08.2015, G.Cosmo
      }
    }
  
  return dist ;
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids. It is usually underestimated from the point of view of
// isotropic safety

G4double 
G4SubtractionSolid::DistanceToIn( const G4ThreeVector& p ) const 
{
  G4double dist=0.0;

#ifdef G4BOOLDEBUG
  if( Inside(p) == kInside )
  {
    G4cout << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToIn(p)" << G4endl
           << "  Point p is inside !" << G4endl;
    G4cerr << "          p = " << p << G4endl;
  }
#endif

  if( ( fPtrSolidA->Inside(p) != kOutside) &&   // case 1
      ( fPtrSolidB->Inside(p) != kOutside)    )
  {
      dist= fPtrSolidB->DistanceToOut(p)  ;
  }
  else
  {
      dist= fPtrSolidA->DistanceToIn(p) ; 
  }
  
  return dist; 
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4SubtractionSolid::DistanceToOut( const G4ThreeVector& p,
                 const G4ThreeVector& v,
                 const G4bool calcNorm,
                       G4bool *validNorm,
                       G4ThreeVector *n ) const 
{
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
             << "G4SubtractionSolid::DistanceToOut(p,v)" << G4endl
             << "  Point p is outside !" << G4endl;
      G4cout << "          p = " << p << G4endl;
      G4cout << "          v = " << v << G4endl;
      G4cerr << "WARNING - Invalid call in "
             << "G4SubtractionSolid::DistanceToOut(p,v)" << G4endl
             << "  Point p is outside !" << G4endl;
      G4cerr << "          p = " << p << G4endl;
      G4cerr << "          v = " << v << G4endl;
    }
#endif

    G4double distout;
    G4double distA = fPtrSolidA->DistanceToOut(p,v,calcNorm,validNorm,n) ;
    G4double distB = fPtrSolidB->DistanceToIn(p,v) ;
    if(distB < distA)
    {
      if(calcNorm)
      {
        *n = -(fPtrSolidB->SurfaceNormal(p+distB*v)) ;
        *validNorm = false ;
      }
      distout= distB ;
    }
    else
    {
      distout= distA ; 
    } 
    return distout;
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4SubtractionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4double dist=0.0;

  if( Inside(p) == kOutside )
  { 
#ifdef G4BOOLDEBUG
    G4cout << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside" << G4endl;
    G4cout << "          p = " << p << G4endl;
    G4cerr << "WARNING - Invalid call in "
           << "G4SubtractionSolid::DistanceToOut(p)" << G4endl
           << "  Point p is outside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
#endif
  }
  else
  {
     dist= std::min(fPtrSolidA->DistanceToOut(p),
                      fPtrSolidB->DistanceToIn(p) ) ; 
  }
  return dist; 
}

//////////////////////////////////////////////////////////////
//
//

G4GeometryType G4SubtractionSolid::GetEntityType() const 
{
  return G4String("G4SubtractionSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4SubtractionSolid::Clone() const
{
  return new G4SubtractionSolid(*this);
}

//////////////////////////////////////////////////////////////
//
//

void 
G4SubtractionSolid::ComputeDimensions(       G4VPVParameterisation*,
                                       const G4int,
                                       const G4VPhysicalVolume* ) 
{
}

/////////////////////////////////////////////////
//
//                    

void 
G4SubtractionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4SubtractionSolid::CreatePolyhedron () const 
{
  HepPolyhedronProcessor processor;
  // Stack components and components of components recursively
  // See G4BooleanSolid::StackPolyhedron
  G4Polyhedron* top = StackPolyhedron(processor, this);
  G4Polyhedron* result = new G4Polyhedron(*top);
  if (processor.execute(*result)) { return result; }
  else { return 0; }
}
