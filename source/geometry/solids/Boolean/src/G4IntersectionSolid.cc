// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 12.09.98 V.Grichine 
// 29.07.99 V.Grichine, modifications in DistanceToIn(p,v)

#include "G4IntersectionSolid.hh"
// #include "G4DisplacedSolid.hh"

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

/////////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid
//

G4IntersectionSolid:: G4IntersectionSolid( const G4String& pName,
                                           G4VSolid* pSolidA ,
                                           G4VSolid* pSolidB   ):
G4BooleanSolid(pName,pSolidA,pSolidB)
{
   ;
} 


///////////////////////////////////////////////////////////////////
//

G4IntersectionSolid:: 
G4IntersectionSolid( const G4String& pName,
                           G4VSolid* pSolidA ,
                           G4VSolid* pSolidB,
                           G4RotationMatrix* rotMatrix,
                     const G4ThreeVector& transVector    ):
G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)
{
   ;
}

//////////////////////////////////////////////////////////////////
//
// 
 
G4IntersectionSolid:: 
G4IntersectionSolid( const G4String& pName,
                           G4VSolid* pSolidA ,
                           G4VSolid* pSolidB ,
                     const G4Transform3D& transform  ):
G4BooleanSolid(pName,pSolidA,pSolidB,transform)
{
   ;
} 


G4IntersectionSolid::~G4IntersectionSolid()
{
    ;
}

///////////////////////////////////////////////////////////////
//
//
     
G4bool 
G4IntersectionSolid::CalculateExtent(const EAxis pAxis,
				     const G4VoxelLimits& pVoxelLimit,
				     const G4AffineTransform& pTransform,
				     G4double& pMin, G4double& pMax) const 
{
  G4bool   retA, retB;
  G4double minA, minB, maxA, maxB; 

  retA= fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, pTransform, minA, maxA);
  retB= fPtrSolidB->CalculateExtent( pAxis, pVoxelLimit, pTransform, minB, maxB);

  pMin = max( minA, minB ); 
  pMax = min( maxA, maxB ); 

  return retA && retB ; // It exists in this slice only if both exist in it.
}
 
/////////////////////////////////////////////////////
//
// Touching ? Empty intersection ?

EInside G4IntersectionSolid::Inside(const G4ThreeVector& p) const
{
  EInside positionA = fPtrSolidA->Inside(p) ;
  EInside positionB = fPtrSolidB->Inside(p) ;
  
  if(positionA == kInside && positionB == kInside)
  {
    return kInside ;
  }
  else
  {
    if((positionA == kInside && positionB == kSurface) ||
       (positionB == kInside && positionA == kSurface) ||
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
G4IntersectionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector normal;
  G4bool insideA, insideB;
  
  insideA= fPtrSolidA->Inside(p);
  insideB= fPtrSolidB->Inside(p);

  // if( Inside(p) == kOutside )
  if( (insideA == kOutside) || (insideB == kOutside) )
  {
    G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p), point p is outside") ;
  }

  // OLD: if(fPtrSolidA->DistanceToOut(p) <= fPtrSolidB->DistanceToOut(p) ) 

  // On the surface of both is difficult ... treat it like on A now!
  //
// if( (insideA == kSurface) && (insideB == kSurface) )
//    normal= fPtrSolidA->SurfaceNormal(p) ;
// else 
  if( insideA == kSurface )
    {
      normal= fPtrSolidA->SurfaceNormal(p) ;
    }
  else if( insideB == kSurface )
    {
      normal= fPtrSolidB->SurfaceNormal(p) ;
    } 
  // We are on neither surface, so we should generate an exception
  else
    {
      G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p), point p is not on the surface of the volume.") ;

      // Or else
      if(fPtrSolidA->DistanceToOut(p) <= fPtrSolidB->DistanceToOut(p) ) 
	 normal= fPtrSolidA->SurfaceNormal(p) ;   
      else
	 normal= fPtrSolidB->SurfaceNormal(p) ;   
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
  G4double dist = 0.0, disTmp = 0.0 ;
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToIn(p,v),  point p is inside") ;
  }
  else
  {
    if( fPtrSolidA->Inside(p) != kOutside    )
    {
      do
      {
        disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
   
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
            disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
            if( disTmp != kInfinity  )
	    {
              dist += disTmp ;
	    }
            else
	    {
              return kInfinity ;
	    }
	  }
          else
	  {
            break ;
	  }
        }
        else
        {
          return kInfinity ;
        } 
      }
      while( Inside(p+dist*v) == kOutside ) ;
      // while( fPtrSolidB->Inside(p) != kOutside ) ;
    }
    else  if(  fPtrSolidB->Inside(p) != kOutside    )
    {     
      do 
      {
        disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
   
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
            disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
            if( disTmp != kInfinity  )
	    {
              dist += disTmp ;
	    }
            else
	    {
              return kInfinity ;
	    }
	  }
          else
	  {
            break ;
	  }
        }
        else
        {
          return kInfinity ;
        } 
      }
      while( Inside(p+dist*v) == kOutside ) ;
    }
    else
    {
      do
      {
        disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
   
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
            disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
            if( disTmp != kInfinity  )
	    {
              dist += disTmp ;
	    }
            else
	    {
              return kInfinity ;
	    }
	  }
          else
	  {
            break ;
	  }
        }
        else
        {
          return kInfinity ;
        } 
      }
      while( Inside(p+dist*v) == kOutside ) ;
    }
  }
  return dist ;  
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids

G4double 
G4IntersectionSolid::DistanceToIn( const G4ThreeVector& p) const 
{
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToIn(p),  point p is inside") ;
  }
  EInside sideA = fPtrSolidA->Inside(p) ;
  EInside sideB = fPtrSolidB->Inside(p) ;
  G4double dist ;

  if( sideA != kInside && sideB  != kOutside    )
  {
    dist = fPtrSolidA->DistanceToIn(p) ;
  }
  else
  {
    if( sideB != kInside  && sideA != kOutside    )
    {
      dist = fPtrSolidB->DistanceToIn(p) ;
    }
    else
    {
      dist =  min(fPtrSolidA->DistanceToIn(p),
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
			            G4ThreeVector *n      ) const 
{
  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p,v),  point p is outside") ;
  }
  G4double distA = fPtrSolidA->DistanceToOut(p,v,calcNorm,validNorm,n) ;
  G4double distB = fPtrSolidB->DistanceToOut(p,v,calcNorm,validNorm,n) ;
  G4double dist = min(distA,distB) ; 
  return dist ; 
  //  return min(fPtrSolidA->DistanceToOut(p,v,calcNorm,validNorm,n),
  //	     fPtrSolidB->DistanceToOut(p,v,calcNorm,validNorm,n) ) ; 
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4IntersectionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p),  point p is outside") ;
  }


  return min(fPtrSolidA->DistanceToOut(p),
             fPtrSolidB->DistanceToOut(p) ) ; 

}

//////////////////////////////////////////////////////////////
//
//

void 
G4IntersectionSolid::ComputeDimensions( G4VPVParameterisation* p,
	                                const G4int n,
                                        const G4VPhysicalVolume* pRep ) 
{
  return ;
}

/////////////////////////////////////////////////
//
//                    

void 
G4IntersectionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  return ;
}

/////////////////////////////////////////////////////////////
//
//

G4VisExtent   
G4IntersectionSolid::GetExtent        () const 
{
  return   G4VisExtent(-1.0,1.0,-1.0,1.0,-1.0,1.0) ;
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4IntersectionSolid::CreatePolyhedron () const 
{
  return new G4PolyhedronBox (1.0, 1.0, 1.0);
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4IntersectionSolid::CreateNURBS      () const 
{
  return new G4NURBSbox (1.0, 1.0, 1.0);
}





