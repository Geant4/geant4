// Implementation for G4DisplacedSolid class for boolean 
// operations between other solids
//
// History:
//
// 28.10.98 V.Grichine, creation according J. Apostolakis's recommendations

#include "G4DisplacedSolid.hh"


#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"
#include "G4VisExtent.hh"

////////////////////////////////////////////////////////////////
//
// Constractor for transformation like rotation of frame then translation 
// in new frame. It is similar to 1st constractor in G4PVPlacement

G4DisplacedSolid::
G4DisplacedSolid( const G4String& pName,
                     G4VSolid* pSolid ,
                     G4RotationMatrix* rotMatrix,
               const G4ThreeVector& transVector    )
   : G4VSolid(pName)
{
  fPtrSolid = pSolid ;
  fPtrTransform = new G4AffineTransform(rotMatrix,transVector) ;
  fPtrTransform->Invert() ;
  fDirectTransform = new G4AffineTransform(rotMatrix,transVector) ;
}

/////////////////////////////////////////////////////////////////////////////////
//
//
G4DisplacedSolid::G4DisplacedSolid( const G4String& pName,
                                    G4VSolid* pSolid ,
			      const G4Transform3D& transform  ) :
  G4VSolid(pName)
{
  fPtrSolid = pSolid ;
  fDirectTransform = new G4AffineTransform(transform.getRotation().inverse(),
                                           transform.getTranslation()) ;

  fPtrTransform    = new G4AffineTransform(transform.getRotation().inverse(),
                                           transform.getTranslation()) ;
  fPtrTransform->Invert() ;
}

/////////////////////////////////////////////////////////////////////////////////
//  Constructor for use with creation of Transient object from Persistent object
//

G4DisplacedSolid::G4DisplacedSolid( const G4String& pName,
                                    G4VSolid* pSolid ,
			      const G4AffineTransform directTransform ) :
  G4VSolid(pName)
{
  fPtrSolid = pSolid ;
  fDirectTransform = new G4AffineTransform( directTransform );
  fPtrTransform    = new G4AffineTransform( directTransform.Inverse() ) ; 
}

///////////////////////////////////////////////////////////////////
//

G4DisplacedSolid::~G4DisplacedSolid() 
{
  if(fPtrTransform)
  {
   delete fPtrTransform ;
   delete fDirectTransform;
  }
}

const G4DisplacedSolid* G4DisplacedSolid::GetDisplacedSolidPtr() const   
{ return this; }

      G4DisplacedSolid* G4DisplacedSolid::GetDisplacedSolidPtr() 
{ return this; }

G4VSolid* G4DisplacedSolid::GetConstituentSolid()
{ return fPtrSolid; } 

const G4AffineTransform  G4DisplacedSolid::GetTransform() const
{
   G4AffineTransform aTransform= *fPtrTransform;
   return aTransform;
}

G4RotationMatrix G4DisplacedSolid::GetFrameRotation() const
{
   G4RotationMatrix InvRotation= fPtrTransform->NetRotation();
   InvRotation.invert();
   return InvRotation;
}

G4ThreeVector  G4DisplacedSolid::GetFrameTranslation() const
{
   return fPtrTransform->NetTranslation();
}

///////////////////////////////////////////////////////////////
G4RotationMatrix G4DisplacedSolid::GetObjectRotation() const
{
   G4RotationMatrix InvRotation= fPtrTransform->NetRotation();
   return InvRotation;
}

G4ThreeVector  G4DisplacedSolid::GetObjectTranslation() const
{
   return fPtrTransform->NetTranslation();
}
///////////////////////////////////////////////////////////////
//
//
     
G4bool 
G4DisplacedSolid::CalculateExtent(const EAxis pAxis,
			       const G4VoxelLimits& pVoxelLimit,
			       const G4AffineTransform& pTransform,
				     G4double& pMin, G4double& pMax) const 
{
  return fPtrSolid->CalculateExtent(pAxis,pVoxelLimit,pTransform,
                                    pMin,pMax) ;
}
 
/////////////////////////////////////////////////////
//
// 

EInside G4DisplacedSolid::Inside(const G4ThreeVector& p) const
{
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  return fPtrSolid->Inside(newPoint) ; 
}

//////////////////////////////////////////////////////////////
//
//

G4ThreeVector 
G4DisplacedSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  G4ThreeVector normal = fPtrSolid->SurfaceNormal(newPoint) ; 
  return fDirectTransform->TransformAxis(normal) ;
    
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4DisplacedSolid::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{    
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  G4ThreeVector newDirection = fPtrTransform->TransformAxis(v) ;
  return fPtrSolid->DistanceToIn(newPoint,newDirection) ;   
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids

G4double 
G4DisplacedSolid::DistanceToIn( const G4ThreeVector& p) const 
{
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  return fPtrSolid->DistanceToIn(newPoint) ;   
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4DisplacedSolid::DistanceToOut( const G4ThreeVector& p,
			            const G4ThreeVector& v,
			            const G4bool calcNorm,
			            G4bool *validNorm,
			            G4ThreeVector *n      ) const 
{
  G4ThreeVector solNorm ; 
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  G4ThreeVector newDirection = fPtrTransform->TransformAxis(v) ;
  G4double dist = fPtrSolid->DistanceToOut(newPoint,newDirection,
                                  calcNorm,validNorm,&solNorm) ;
  if(calcNorm)
  { 
    *n = fDirectTransform->TransformAxis(solNorm) ;
  }
  return dist ;  
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4DisplacedSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4ThreeVector newPoint = fPtrTransform->TransformPoint(p) ;
  return fPtrSolid->DistanceToOut(newPoint) ;   
}

//////////////////////////////////////////////////////////////
//
//

void 
G4DisplacedSolid::ComputeDimensions( G4VPVParameterisation* p,
	                                const G4int n,
                                        const G4VPhysicalVolume* pRep ) 
{
  // fPtrSolid->ComputeDimensions(p,n,pRep);

  G4Exception("ERROR: ComputeDimensions has no meaning for a G4DisplacedSolid. It cannot be called.");
}

/////////////////////////////////////////////////
//
//                    

void 
G4DisplacedSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  fPtrSolid->DescribeYourselfTo(scene) ;
}

/////////////////////////////////////////////////////////////
//
//

G4VisExtent   
G4DisplacedSolid::GetExtent        () const 
{
  return fPtrSolid->GetExtent()  ;
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4DisplacedSolid::CreatePolyhedron () const 
{
  return fPtrSolid->CreatePolyhedron() ;
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4DisplacedSolid::CreateNURBS      () const 
{
  return fPtrSolid->CreateNURBS() ;
}





