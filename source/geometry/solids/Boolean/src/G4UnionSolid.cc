// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 12.09.98 V.Grichine 

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
#include "G4VisExtent.hh"

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

  pMin = min( minA, minB ); 
  pMax = max( maxA, maxB ); 

  return touchesA || touchesB ;  // It exists in this slice if either one does.
}
 
/////////////////////////////////////////////////////
//
// Important comment: When solids A and B touch together along flat
// surface the surface points will be considered as kSurface. To avoid
// infinite loops it is better to combine solids with intersection/gap
// thicker than 3 kCarTolerance

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
       G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p),  point p is outside") ;
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
      G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p),  point p is inside") ;
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
  G4double dist ;
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToIn(p,v),  point p is inside") ;
  }
  return min(fPtrSolidA->DistanceToIn(p,v),
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

  return min(distA,distB) ;
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
  G4double disTmp = 0.0, dist = 0.0 ;
  G4ThreeVector normTmp;
  G4ThreeVector* nTmp= &normTmp;

  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p,v),  point p is outside") ;
  }
  else
  {
    EInside positionA = fPtrSolidA->Inside(p) ;
    EInside positionB = fPtrSolidB->Inside(p) ;
  
    if( positionA != kOutside )
    { 
      do
      {
        disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp)        ;
        dist += disTmp ;
        if( Inside(p+dist*v) == kInside )
	{ 
          disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                            validNorm,nTmp)         ;
          dist += disTmp ;
	}
        else
	{
          break ;
	}
      }
      while( Inside(p+dist*v) == kInside ) ;
      *n = *nTmp ; 
    }
    else
    {
      do
      {
        disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v,calcNorm,
                                           validNorm,nTmp)        ;
        dist += disTmp ;
        if( Inside(p+dist*v) == kInside )
	{ 
          disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v,calcNorm,
                                            validNorm,nTmp)         ;
          dist += disTmp ;
	}
        else
	{
          break ;
	}
      }
      while( Inside(p+dist*v) == kInside ) ;
      *n = *nTmp ;   
    }
  }
  *validNorm = false ;
  return dist ;
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4UnionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4double distout;
  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p),  point p is outside") ;
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
      distout= max(fPtrSolidA->DistanceToOut(p),
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
  return ;
}

/////////////////////////////////////////////////////////////
//
//

G4VisExtent   
G4UnionSolid::GetExtent        () const 
{
  return   G4VisExtent(-1.0,1.0,-1.0,1.0,-1.0,1.0) ;
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4UnionSolid::CreatePolyhedron () const 
{
  return new G4PolyhedronBox (1.0, 1.0, 1.0);
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4UnionSolid::CreateNURBS      () const 
{
  return new G4NURBSbox (1.0, 1.0, 1.0);
}





