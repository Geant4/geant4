// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SubtractionSolid.cc,v 1.5 2000-04-27 09:59:40 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 14.10.98 V.Grichine, implementation of the first version 
// 19.10.98 V.Grichine, new algorithm of DistanceToIn(p,v) according to
//          J.Apostolakis recommendations (while loops)
// 02.08.99 V.Grichine, bugs fixed in DistanceToOut(p,v,...)
//                      while -> do-while & surfaceA limitations
//

#include "G4SubtractionSolid.hh"
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

///////////////////////////////////////////////////////////////////
//
// Transfer all data members to G4BooleanSolid which is responsible
// for them. pName will be in turn sent to G4VSolid

G4SubtractionSolid:: G4SubtractionSolid( const G4String& pName,
                                           G4VSolid* pSolidA ,
                                           G4VSolid* pSolidB   ):
G4BooleanSolid(pName,pSolidA,pSolidB)
{
   ;
}

///////////////////////////////////////////////////////////////
//
//
 
G4SubtractionSolid:: 
G4SubtractionSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB ,
                          G4RotationMatrix* rotMatrix,
                    const G4ThreeVector& transVector    ):
G4BooleanSolid(pName,pSolidA,pSolidB,rotMatrix,transVector)
{
   ;
} 

///////////////////////////////////////////////////////////////
//
//

G4SubtractionSolid:: 
G4SubtractionSolid( const G4String& pName,
                          G4VSolid* pSolidA ,
                          G4VSolid* pSolidB ,
                    const G4Transform3D& transform  ):
G4BooleanSolid(pName,pSolidA,pSolidB,transform)
{
   ;
} 


G4SubtractionSolid::~G4SubtractionSolid()
{
    ;
}

///////////////////////////////////////////////////////////////
//
//
     
G4bool 
G4SubtractionSolid::CalculateExtent(const EAxis pAxis,
				     const G4VoxelLimits& pVoxelLimit,
				     const G4AffineTransform& pTransform,
				     G4double& pMin, G4double& pMax) const 
{
  // Since we cannot be sure how much the second solid subtracts 
  // from the first,    we must use the first solid's extent!

  return fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, 
                                      pTransform, pMin, pMax);
}
 
/////////////////////////////////////////////////////
//
// Touching ? Empty subtraction ?

EInside G4SubtractionSolid::Inside(const G4ThreeVector& p) const
{
  EInside positionA = fPtrSolidA->Inside(p) ;
  EInside positionB = fPtrSolidB->Inside(p) ;
  
  if(positionA == kInside && positionB == kOutside)
  {
    return kInside ;
  }
  else
  {
    if((positionA == kInside && positionB == kSurface) ||
       (positionB == kOutside && positionA == kSurface) ||
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
G4SubtractionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector normal;
  if( Inside(p) == kOutside )
  {
    G4Exception("Invalid call in G4SubtractionSolid::SurfaceNormal(p), point p is outside") ;
  }
  else
  { 
    if( fPtrSolidA->Inside(p) == kSurface && 
        fPtrSolidB->Inside(p) != kInside      ) 
    {
      normal= fPtrSolidA->SurfaceNormal(p) ;
    }
    else if( fPtrSolidA->Inside(p) == kInside && 
             fPtrSolidB->Inside(p) != kOutside    )
    {
      normal= -fPtrSolidB->SurfaceNormal(p) ;
    }
    else 
    {
      G4Exception("Invalid call in G4SubtractionSolid::SurfaceNormal(p),  point p is inside rather than on surface") ;
    }
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
  G4double dist = 0.0,disTmp = 0.0 ;
    
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4SubtractionSolid::DistanceToIn(p,v),  point p is inside") ;
  }

  if( // ( fPtrSolidA->Inside(p) != kOutside) &&     // case1:p in both A&B 

      ( fPtrSolidB->Inside(p) != kOutside)    )   // start: out of B
    {
      dist = fPtrSolidB->DistanceToOut(p,v) ; // ,calcNorm,validNorm,n) ;
      
      if( fPtrSolidA->Inside(p+dist*v) != kInside )
      {
        do
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
            dist += disTmp ;
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
	 while( Inside(p+dist*v) == kOutside )  // pushing loop
	//  do
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
            dist += disTmp ;
	  }
	}
	//  while( Inside(p+dist*v) == kOutside ) ;
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
G4SubtractionSolid::DistanceToIn( const G4ThreeVector& p) const 
{
  G4double dist;
  if( Inside(p) == kInside )
  {
    G4Exception("Invalid call in G4SubtractionSolid::DistanceToIn(p), point p is inside") ;
  }

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
			           G4ThreeVector *n      ) const 
{
    if( Inside(p) == kOutside )
    {
         G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p,v), point p is outside") ;
    }

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
      *validNorm = true ;
    } 
    return distout;
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4SubtractionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4double dist;

  if( Inside(p) == kOutside )
  { 
    G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p),  point p is outside") ;
  }
  else
  {
     dist= G4std::min(fPtrSolidA->DistanceToOut(p),
               fPtrSolidB->DistanceToIn(p) ) ; 
  }
  return dist; 
}

//////////////////////////////////////////////////////////////
//
//

void 
G4SubtractionSolid::ComputeDimensions( G4VPVParameterisation* p,
	                                const G4int n,
                                        const G4VPhysicalVolume* pRep ) 
{
  return ;
}

/////////////////////////////////////////////////
//
//                    

void 
G4SubtractionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddThis (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4SubtractionSolid::CreatePolyhedron () const 
{
  G4Polyhedron* pA = fPtrSolidA->CreatePolyhedron();
  G4Polyhedron* pB = fPtrSolidB->CreatePolyhedron();
  G4Polyhedron* resultant = new G4Polyhedron (pA->subtract(*pB));
  delete pB;
  delete pA;
  return resultant;
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4SubtractionSolid::CreateNURBS      () const 
{
  // Take into account boolean operation - see CreatePolyhedron.
  // return new G4NURBSbox (1.0, 1.0, 1.0);
  return 0;
}
