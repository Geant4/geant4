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
// $Id: G4IntersectionSolid.cc,v 1.16 2001-10-02 08:51:48 grichine Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation of methods for the class G4IntersectionSolid
//
// History:
//
// 29.05.01 V.Grichine, bug was fixed in DistanceToIn(p,v)
// 16.03.01 V.Grichine, modifications in CalculateExtent() and Inside()
//                      based on D.Williams proposal
// 29.07.99 V.Grichine, modifications in DistanceToIn(p,v)
// 12.09.98 V.Grichine, implementation based on discussions with J. Apostolakis and 
//                      S. Giani 

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
  G4bool   retA, retB, out ;
  G4double minA, minB, maxA, maxB ; 

  retA= fPtrSolidA->CalculateExtent( pAxis, pVoxelLimit, pTransform, minA, maxA);
  retB= fPtrSolidB->CalculateExtent( pAxis, pVoxelLimit, pTransform, minB, maxB);
  if(retA && retB)
  {
    pMin = G4std::max( minA, minB ) ; 
    pMax = G4std::min( maxA, maxB ) ;
    out  = true ;
  }
  else out = false ;

  return out ; // It exists in this slice only if both exist in it.
}
 
/////////////////////////////////////////////////////
//
// Touching ? Empty intersection ?

EInside G4IntersectionSolid::Inside(const G4ThreeVector& p) const
{
  EInside positionA = fPtrSolidA->Inside(p) ;

  if( positionA == kOutside ) return kOutside ;

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

G4ThreeVector 
G4IntersectionSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector normal;
  EInside insideA, insideB;
  
  insideA= fPtrSolidA->Inside(p);
  insideB= fPtrSolidB->Inside(p);

  // if( Inside(p) == kOutside )
  if( (insideA == kOutside) || (insideB == kOutside) )
  {
    G4cerr << "WARNING - Invalid call in G4IntersectionSolid::SurfaceNormal(p),"
           << " point p is outside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    #ifdef G4BOOLDEBUG
       G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p), p is outside") ;
    #endif
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
      if(fPtrSolidA->DistanceToOut(p) <= fPtrSolidB->DistanceToOut(p) ) 
	 normal= fPtrSolidA->SurfaceNormal(p) ;   
      else
	 normal= fPtrSolidB->SurfaceNormal(p) ;   
      G4cerr << "WARNING - Invalid call in G4IntersectionSolid::SurfaceNormal(p),"
             << " point p is outsurface" << G4endl;
      G4cerr << "          p = " << p << G4endl;
      #ifdef G4BOOLDEBUG
         G4Exception("Invalid call in G4IntersectionSolid::SurfaceNormal(p), p is outsurface") ;
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
  G4double dist = 0.0, disTmp = 0.0 ;
  if( Inside(p) == kInside )
  {
    G4cerr << "WARNING - Invalid call in G4IntersectionSolid::DistanceToIn(p,v),"
           << " point p is inside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
    #ifdef G4BOOLDEBUG
      G4Exception("Invalid call in G4IntersectionSolid::DistanceToIn(p,v), p is inside") ;
    #endif
  }
  else
  {
    if( fPtrSolidA->Inside(p) != kOutside    )
    {
      do
      {
	if( fPtrSolidB->Inside(p+dist*v) == kInside    )
	{
	  disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v) ; 
	}
	else
	{
          disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
	}
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
	    if( fPtrSolidA->Inside(p+dist*v) == kInside    )
	    {
	      disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v) ; 
	    }
	    else
	    {
              disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
	    }
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
	if( fPtrSolidA->Inside(p+dist*v) == kInside    )
	{
	  disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v) ; 
	}
	else
	{
          disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
	}   
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
	    if( fPtrSolidB->Inside(p+dist*v) == kInside    )
	    {
	      disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v) ; 
	    }
	    else
	    {
              disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
	    }
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
	if( fPtrSolidB->Inside(p+dist*v) == kInside    )
	{
	  disTmp = fPtrSolidB->DistanceToOut(p+dist*v,v) ; 
	}
	else
	{
          disTmp = fPtrSolidB->DistanceToIn(p+dist*v,v) ;
	}
        if( disTmp != kInfinity )
        {
          dist += disTmp ;
          if(Inside(p+dist*v) == kOutside  )
	  {
	    if( fPtrSolidA->Inside(p+dist*v) == kInside    )
	    {
	      disTmp = fPtrSolidA->DistanceToOut(p+dist*v,v) ; 
	    }
	    else
	    {
              disTmp = fPtrSolidA->DistanceToIn(p+dist*v,v) ;
	    }
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
    G4cerr << "WARNING - Invalid call in G4IntersectionSolid::DistanceToIn(p),"
           << " point p is inside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    #ifdef G4BOOLDEBUG
      G4Exception("Invalid call in G4IntersectionSolid::DistanceToIn(p), p is inside") ;
    #endif
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
      dist =  G4std::min(fPtrSolidA->DistanceToIn(p),
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
  G4bool         validNormA, validNormB;
  G4ThreeVector  nA, nB;

  if( Inside(p) == kOutside )
  {
    G4cerr << "WARNING - Invalid call in G4IntersectionSolid::DistanceToOut(p,v),"
           << " point p is outside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    G4cerr << "          v = " << v << G4endl;
    G4cout << "Position:"  << G4endl << G4endl;
    G4cout << "p.x() = "   << p.x()/mm << " mm" << G4endl;
    G4cout << "p.y() = "   << p.y()/mm << " mm" << G4endl;
    G4cout << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl;
    G4cout << "Direction:" << G4endl << G4endl;
    G4cout << "v.x() = "   << v.x() << G4endl;
    G4cout << "v.y() = "   << v.y() << G4endl;
    G4cout << "v.z() = "   << v.z() << G4endl << G4endl;
    #ifdef G4BOOLDEBUG
      G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p,v), p is outside") ;
    #endif
  }
  G4double distA = fPtrSolidA->DistanceToOut(p,v,calcNorm,&validNormA,&nA) ;
  G4double distB = fPtrSolidB->DistanceToOut(p,v,calcNorm,&validNormB,&nB) ;

  G4double dist = G4std::min(distA,distB) ; 

  if( calcNorm ){
     if (distA < distB ) {
        *validNorm = validNormA;
	*n =         nA;
     }else{   
        *validNorm = validNormB;
        *n =         nB;
     }
  }

  return dist ; 
  //  return G4std::min(fPtrSolidA->DistanceToOut(p,v,calcNorm,validNorm,n),
  //	                fPtrSolidB->DistanceToOut(p,v,calcNorm,validNorm,n) ) ; 
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4IntersectionSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  if( Inside(p) == kOutside )
  {
    G4cerr << "WARNING - Invalid call in G4IntersectionSolid::DistanceToOut(p),"
           << " point p is outside" << G4endl;
    G4cerr << "          p = " << p << G4endl;
    #ifdef G4BOOLDEBUG
      G4Exception("Invalid call in G4IntersectionSolid::DistanceToOut(p), p is outside") ;
    #endif
  }


  return G4std::min(fPtrSolidA->DistanceToOut(p),
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

G4GeometryType G4IntersectionSolid::GetEntityType() const 
{
  return G4String("G4IntersectionSolid");
}

/////////////////////////////////////////////////
//
//                    

void 
G4IntersectionSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddThis (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4IntersectionSolid::CreatePolyhedron () const 
{
  G4Polyhedron* pA = fPtrSolidA->CreatePolyhedron();
  G4Polyhedron* pB = fPtrSolidB->CreatePolyhedron();
  G4Polyhedron* resultant = new G4Polyhedron (pA->intersect(*pB));
  delete pB;
  delete pA;
  return resultant;
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4IntersectionSolid::CreateNURBS      () const 
{
  // Take into account boolean operation - see CreatePolyhedron.
  // return new G4NURBSbox (1.0, 1.0, 1.0);
  return 0;
}
