// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4UnionSolid.cc,v 1.13 2001-01-23 08:26:44 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 12.09.98 V.Grichine, implementation
// 28.11.98 V.Grichine, J. Apostolakis, while loops in DistToIn/Out 
// 27.07.99 V.Grichine, modifications in Distance ToOut(p,v,...)
//                      while -> do-while

#include "G4UnionSolid.hh"
// #include "G4PlacedSolid.hh"

#include "G4RotationMatrix.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"
#include "G4AffineTransform.hh"

#include "G4VoxelLimits.hh"
#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"

///////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid

G4UnionSolid:: G4UnionSolid( const G4String& pName,
                                           G4VSolid* pSolidA ,
                                           G4VSolid* pSolidB   ):
G4BooleanSolid(pName,pSolidA,pSolidB)
{
   ;
}

/////////////////////////////////////////////////////////////////////
//
//
 
G4UnionSolid:: G4UnionSolid( const G4String& pName,
                                   G4VSolid* pSolidA ,
                                   G4VSolid* pSolidB ,
                                   G4RotationMatrix* rotMatrix,
                             const G4ThreeVector& transVector    ):
G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)

{
   ;
}

///////////////////////////////////////////////////////////
//
//
 
G4UnionSolid:: G4UnionSolid( const G4String& pName,
                                   G4VSolid* pSolidA ,
                                   G4VSolid* pSolidB ,
                             const G4Transform3D& transform  ):
G4BooleanSolid(pName,pSolidA,pSolidB,transform)

{
   ;
} 



G4UnionSolid::~G4UnionSolid()
{
    ;
}

///////////////////////////////////////////////////////////////
//
//
     
G4bool 
G4UnionSolid::CalculateExtent(const EAxis pAxis,
				     const G4VoxelLimits& pVoxelLimit,
				     const G4AffineTransform& pTransform,
				     G4double& pMin, G4double& pMax) const 
{
  G4bool   touchesA, touchesB;
  G4double minA, minB, maxA, maxB; 

  touchesA= fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, pTransform, minA, maxA);
  touchesB= fPtrSolidB->CalculateExtent( pAxis, pVoxelLimit, pTransform, minB, maxB);

  pMin = G4std::min( minA, minB ); 
  pMax = G4std::max( maxA, maxB ); 

  return touchesA || touchesB ;  // It exists in this slice if either one does.
}
 
/////////////////////////////////////////////////////
//
// Important comment: When solids A and B touch together along flat
// surface the surface points will be considered as kSurface, while points 
// located around will correspond to kInside

EInside G4UnionSolid::Inside(const G4ThreeVector& p) const
{
  EInside positionA = fPtrSolidA->Inside(p) ;
  EInside positionB = fPtrSolidB->Inside(p) ;
  
  if( positionA == kInside  || positionB == kInside )
  {    
    return kInside ;
  }
  else
  {
    if((positionA != kInside && positionB == kSurface) ||
       (positionB != kInside && positionA == kSurface) ||
       (positionA == kSurface && positionB == kSurface)   )
    {
      return kSurface ;
    }
    else
    {
      return kOutside ;
    }
  }
}

//////////////////////////////////////////////////////////////
//
//

G4ThreeVector 
G4UnionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
    G4ThreeVector normal;

    if( Inside(p) == kOutside )
    {
       G4Exception("Invalid call in G4UnionSolid::SurfaceNormal(p),  point p is outside") ;
    }

    if(fPtrSolidA->Inside(p) == kSurface && fPtrSolidB->Inside(p) != kInside) 
    {
       normal= fPtrSolidA->SurfaceNormal(p) ;
    }
    else if(fPtrSolidB->Inside(p) == kSurface && 
            fPtrSolidA->Inside(p) != kInside)
    {
       normal= fPtrSolidB->SurfaceNormal(p) ;
    }
    else 
    {
      G4Exception("Invalid call in G4UnionSolid::SurfaceNormal(p),  point p is inside") ;
    }

    return normal;
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4UnionSolid::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4UnionSolid::DistanceToIn(p,v),  point p is inside") ;
  }
  return G4std::min(fPtrSolidA->DistanceToIn(p,v),
                    fPtrSolidB->DistanceToIn(p,v) ) ;
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the union of
// two solids

G4double 
G4UnionSolid::DistanceToIn( const G4ThreeVector& p) const 
{
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4UnionSolid::DistanceToIn(p),  point p is inside") ;
  }
  G4double distA = fPtrSolidA->DistanceToIn(p) ;
  G4double distB = fPtrSolidB->DistanceToIn(p) ;

  return G4std::min(distA,distB) ;
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
     G4cout << "Position:"  << G4endl << G4endl;
     G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
     G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
     G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
     G4cout << "Direction:" << G4endl << G4endl;
     G4cout << "v.x() = "   << v.x() << G4endl;
     G4cout << "v.y() = "   << v.y() << G4endl;
     G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
    G4Exception("Invalid call in G4UnionSolid::DistanceToOut(p,v),  point p is outside") ;
  }
  else
  {
    EInside positionA = fPtrSolidA->Inside(p) ;
    // EInside positionB = fPtrSolidB->Inside(p) ;

    if( positionA != kOutside )
    { 
      do
      {
        disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp)        ;
        dist += disTmp ;

        if(fPtrSolidB->Inside(p+dist*v) != kOutside)
	{ 
          disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                            validNorm,nTmp)         ;
          dist += disTmp ;
	}
      }
      //     while( Inside(p+dist*v) == kInside ) ;
           while( fPtrSolidA->Inside(p+dist*v) != kOutside && 
                  disTmp > 0.5*kCarTolerance ) ;
      *n = *nTmp ; 
    }
    else // if( positionB != kOutside )
    {
      do
      {
        disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp)        ; 
        dist += disTmp ;

        if(fPtrSolidA->Inside(p+dist*v) != kOutside)
	{ 
          disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                            validNorm,nTmp)         ;
          dist += disTmp ;
	}
      }
      //  while( Inside(p+dist*v) == kInside ) ;
        while( fPtrSolidB->Inside(p+dist*v) != kOutside && 
              disTmp > 0.5*kCarTolerance ) ;
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
  G4double distout = kInfinity;
  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4UnionSolid::DistanceToOut(p),  point p is outside") ;
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
      distout= G4std::max(fPtrSolidA->DistanceToOut(p),
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

//////////////////////////////////////////////////////////////
//
//

void 
G4UnionSolid::ComputeDimensions( G4VPVParameterisation* p,
	                                const G4int n,
                                        const G4VPhysicalVolume* pRep ) 
{
  return ;
}

/////////////////////////////////////////////////
//
//                    

void 
G4UnionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddThis (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4UnionSolid::CreatePolyhedron () const 
{
  G4Polyhedron* pA = fPtrSolidA->CreatePolyhedron();
  G4Polyhedron* pB = fPtrSolidB->CreatePolyhedron();
  G4Polyhedron* resultant = new G4Polyhedron (pA->add(*pB));
  delete pB;
  delete pA;
  return resultant;
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4UnionSolid::CreateNURBS      () const 
{
  // Take into account boolean operation - see CreatePolyhedron.
  // return new G4NURBSbox (1.0, 1.0, 1.0);
  return 0;
}
