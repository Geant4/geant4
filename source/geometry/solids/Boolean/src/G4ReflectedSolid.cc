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
// $Id: G4ReflectedSolid.cc,v 1.7 2002-04-16 10:14:22 grichine Exp $
//
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Implementation for G4ReflectedSolid class for boolean 
// operations between other solids
//
// Author: Vladimir Grichine, 23.07.01  (Vladimir.Grichine@cern.ch)

#include "G4ReflectedSolid.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4VSolid.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBSbox.hh"


/////////////////////////////////////////////////////////////////////////////////
//
// Constructor using HepTransform3D, in fact HepReflect3D

G4ReflectedSolid::G4ReflectedSolid( const G4String& pName,
                                          G4VSolid* pSolid ,
			            const G4Transform3D& transform  ) :
  G4VSolid(pName)
{
  fPtrSolid = pSolid ;
  G4RotationMatrix rotMatrix ;
  
  fDirectTransform = new G4AffineTransform(rotMatrix, transform.getTranslation()) ;  

  fPtrTransform    = new G4AffineTransform(rotMatrix, transform.getTranslation()) ; 
  fPtrTransform->Invert() ;

  fDirectTransform3D = new G4Transform3D(transform) ;
  fPtrTransform3D    = new G4Transform3D(transform.inverse()) ;   
}

/* **************************************************************

////////////////////////////////////////////////////////////////
//
// Constractor for transformation like rotation of frame then translation 
// in new frame. It is similar to 1st constractor in G4PVPlacement

G4ReflectedSolid::
G4ReflectedSolid( const G4String& pName,
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
//  Constructor for use with creation of Transient object from Persistent object
//

G4ReflectedSolid::G4ReflectedSolid( const G4String& pName,
                                    G4VSolid* pSolid ,
			      const G4AffineTransform directTransform ) :
  G4VSolid(pName)
{
  fPtrSolid = pSolid ;
  fDirectTransform = new G4AffineTransform( directTransform );
  fPtrTransform    = new G4AffineTransform( directTransform.Inverse() ) ; 
}

********************************************************** */


///////////////////////////////////////////////////////////////////
//

G4ReflectedSolid::~G4ReflectedSolid() 
{
  if(fPtrTransform)
  {
   delete fPtrTransform ;
   delete fDirectTransform;
  }
}

G4GeometryType G4ReflectedSolid::GetEntityType() const 
{
   return G4String("G4ReflectedSolid");
}

const G4ReflectedSolid* G4ReflectedSolid::GetReflectedSolidPtr() const   
{ return this; }

      G4ReflectedSolid* G4ReflectedSolid::GetReflectedSolidPtr() 
{ return this; }

G4VSolid* G4ReflectedSolid::GetConstituentMovedSolid() const
{ 
  return fPtrSolid; 
} 

/////////////////////////////////////////////////////////////////////////////

G4AffineTransform  G4ReflectedSolid::GetTransform() const
{
   G4AffineTransform aTransform = *fPtrTransform;
   return aTransform;
}

void G4ReflectedSolid::SetTransform(G4AffineTransform& transform) 
{
   fPtrTransform = &transform ;
}

//////////////////////////////////////////////////////////////////////////////

G4AffineTransform  G4ReflectedSolid::GetDirectTransform() const
{
   G4AffineTransform aTransform= *fDirectTransform;
   return aTransform;
}

void G4ReflectedSolid::SetDirectTransform(G4AffineTransform& transform) 
{
   fDirectTransform = &transform ;
}

/////////////////////////////////////////////////////////////////////////////





/////////////////////////////////////////////////////////////////////////////

G4Transform3D  G4ReflectedSolid::GetTransform3D() const
{
   G4Transform3D aTransform = *fPtrTransform3D;
   return aTransform;
}

void G4ReflectedSolid::SetTransform3D(G4Transform3D& transform) 
{
   fPtrTransform3D = &transform ;
}

//////////////////////////////////////////////////////////////////////////////

G4Transform3D  G4ReflectedSolid::GetDirectTransform3D() const
{
   G4Transform3D aTransform= *fDirectTransform3D;
   return aTransform;
}

void G4ReflectedSolid::SetDirectTransform3D(G4Transform3D& transform) 
{
   fDirectTransform3D = &transform ;
}

/////////////////////////////////////////////////////////////////////////////

G4RotationMatrix G4ReflectedSolid::GetFrameRotation() const
{
   G4RotationMatrix InvRotation= fDirectTransform->NetRotation();
   return InvRotation;
}

void G4ReflectedSolid::SetFrameRotation(const G4RotationMatrix& matrix)
{
   fDirectTransform->SetNetRotation(matrix);
}

/////////////////////////////////////////////////////////////////////////////

G4ThreeVector  G4ReflectedSolid::GetFrameTranslation() const
{
   return fPtrTransform->NetTranslation();
}

void G4ReflectedSolid::SetFrameTranslation(const G4ThreeVector& vector)
{
  fPtrTransform->SetNetTranslation(vector);
}

///////////////////////////////////////////////////////////////

G4RotationMatrix G4ReflectedSolid::GetObjectRotation() const
{
   G4RotationMatrix Rotation= fPtrTransform->NetRotation();
   return Rotation;
}

void G4ReflectedSolid::SetObjectRotation(const G4RotationMatrix& matrix)
{
   fPtrTransform->SetNetRotation(matrix);
}

///////////////////////////////////////////////////////////////////////

G4ThreeVector  G4ReflectedSolid::GetObjectTranslation() const
{
   return fDirectTransform->NetTranslation();
}

void G4ReflectedSolid::SetObjectTranslation(const G4ThreeVector& vector)
{
  fDirectTransform->SetNetTranslation(vector);
}

///////////////////////////////////////////////////////////////
//
//
     
G4bool 
G4ReflectedSolid::CalculateExtent( const EAxis pAxis,
			           const G4VoxelLimits& pVoxelLimit,
			           const G4AffineTransform& pTransform,
				         G4double& pMin, 
                                         G4double& pMax           ) const 
{

  G4AffineTransform sumTransform ;
  G4Transform3D sumT, pt3d ;

  G4ReflectX3D tX ;
  G4ReflectY3D tY ;
  G4ReflectZ3D tZ ;

  G4bool extentR ;
  G4double minR, maxR ;

  pt3d = G4Transform3D(pTransform.NetRotation().inverse(),
                      pTransform.NetTranslation()          );
  sumT = pt3d * (*fDirectTransform3D) ;
  sumT = ((sumT*tX)*tY)*tZ;
  sumTransform = G4AffineTransform( sumT.getRotation(),
                                    // sumT.getRotation().inverse(),
                                    sumT.getTranslation()         );


  //  sumTransform.Product(*fDirectTransform,pTransform) ;
  //  extent  = fPtrSolid->CalculateExtent(pAxis,pVoxelLimit,pTransform,min,max) ;
  extentR = fPtrSolid->CalculateExtent(pAxis,pVoxelLimit,sumTransform,
                                       minR,maxR) ;

  // Additional extension, just in case:
  /*  
  if( maxR > 0 ) pMax = 2.0*maxR ;
  else           pMax = 0.5*maxR ;
  if( minR > 0 ) pMin = 0.5*minR ;
  else           pMin = 2.0*minR ;
  maxR =  pMax;
  minR =  pMin;
  */
  // and John's trick:
  pMin = -maxR ;  // minR ; 
  pMax = -minR ;  // maxR ; 
  return extentR ;

  /* ******************************************

  // General rotated case from G4Tube

  G4int i, noEntries, noBetweenSections4 ;
  G4bool existsAfterClip = false ;
  G4Point3D refPoint;
  G4ThreeVector point; 
  G4ThreeVectorList* vertices    = new G4ThreeVectorList() ; 
  G4ThreeVectorList* tmpVertices = fPtrSolid->CreateRotatedVertices(pTransform) ;
    
  pMin      = +kInfinity ;
  pMax      = -kInfinity ;
  noEntries = tmpVertices->size() ;

  for(i = 0 ; i < noEntries ; i++ )
  {
      point = (*tmpVertices)[i];
      //  G4Point3D refPoint = (*fDirectTransform3D)*G4Point3D(point) ;
      G4Point3D refPoint = (*fPtrTransform3D)*G4Point3D(point) ;
      point = G4ThreeVector(refPoint.x(),refPoint.y(),refPoint.z()); 
      vertices->push_back(point);     
  }	    
  noBetweenSections4 = noEntries - 4 ;
  
  G4cout<<G4endl;
  for(i = 0 ; i < noBetweenSections4 ; i++ )
  {
    G4cout<<i<<"\t"<<"v.x = "<<((*vertices)[i]).x()<<"\t"
                   <<"v.y = "<<((*vertices)[i]).y()<<"\t"
                   <<"v.z = "<<((*vertices)[i]).z()<<"\t"
          <<G4endl;
  }	    
  G4cout<<G4endl;
  
  for (i = 0 ; i < noEntries ; i += 4 )
  {
    ClipCrossSection(vertices,i,pVoxelLimit,pAxis,pMin,pMax) ;
  }
  for (i = 0 ; i < noBetweenSections4 ; i += 4 )
  {
    ClipBetweenSections(vertices,i,pVoxelLimit,pAxis,pMin,pMax) ;
  }
  if (pMin != kInfinity || pMax != -kInfinity )
  {
    existsAfterClip = true ;
    pMin -= kCarTolerance ; // Add 2*tolerance to avoid precision troubles
    pMax += kCarTolerance ;
  }
  else
  {
// Check for case where completely enveloping clipping volume
// If point inside then we are confident that the solid completely
// envelopes the clipping volume. Hence set min/max extents according
// to clipping volume extents along the specified axis.

    G4ThreeVector clipCentre(

    (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
    (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
    (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5 ) ;
		    
    if ( Inside(pTransform.Inverse().TransformPoint(clipCentre)) != kOutside )
    {
      existsAfterClip = true ;
      pMin            = pVoxelLimit.GetMinExtent(pAxis) ;
      pMax            = pVoxelLimit.GetMaxExtent(pAxis) ;
    }
  }
  delete vertices;
  delete tmpVertices;

  return existsAfterClip;

  ****************************************** */
}

 
/////////////////////////////////////////////////////
//
// 

EInside G4ReflectedSolid::Inside(const G4ThreeVector& p) const
{

  G4Point3D newPoint = (*fDirectTransform3D)*G4Point3D(p) ;
  // G4Point3D newPoint = (*fPtrTransform3D)*G4Point3D(p) ;

  return fPtrSolid->Inside(G4ThreeVector(newPoint.x(),newPoint.y(),newPoint.z())) ; 
}

//////////////////////////////////////////////////////////////
//
//

G4ThreeVector 
G4ReflectedSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4Point3D newPoint = (*fDirectTransform3D)*G4Point3D(p) ;
  G4ThreeVector normal = fPtrSolid->SurfaceNormal(G4ThreeVector(newPoint.x(),
                                                                newPoint.y(),
                                                                newPoint.z() ) ) ;
  G4Point3D newN = (*fDirectTransform3D)*G4Point3D(normal) ;
  newN.unit() ;

  return G4ThreeVector(newN.x(),newN.y(),newN.z()) ;
    
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4ReflectedSolid::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v  ) const 
{    
  G4Point3D newPoint     = (*fDirectTransform3D)*G4Point3D(p) ;
  G4Point3D newDirection = (*fDirectTransform3D)*G4Point3D(v) ;
  newDirection.unit() ;
  return fPtrSolid->DistanceToIn(
       G4ThreeVector(newPoint.x(),newPoint.y(),newPoint.z()),
       G4ThreeVector(newDirection.x(),newDirection.y(),newDirection.z())) ;   
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids

G4double 
G4ReflectedSolid::DistanceToIn( const G4ThreeVector& p) const 
{
  G4Point3D newPoint = (*fDirectTransform3D)*G4Point3D(p) ;
  return fPtrSolid->DistanceToIn(
                    G4ThreeVector(newPoint.x(),newPoint.y(),newPoint.z())) ;   
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4ReflectedSolid::DistanceToOut( const G4ThreeVector& p,
			            const G4ThreeVector& v,
			            const G4bool calcNorm,
			            G4bool *validNorm,
			            G4ThreeVector *n      ) const 
{
  G4ThreeVector solNorm ; 

  G4Point3D newPoint     = (*fDirectTransform3D)*G4Point3D(p) ;
  G4Point3D newDirection = (*fDirectTransform3D)*G4Point3D(v) ;
  newDirection.unit() ;

  G4double dist = fPtrSolid->DistanceToOut(
                  G4ThreeVector(newPoint.x(),newPoint.y(),newPoint.z()),
                  G4ThreeVector(newDirection.x(),newDirection.y(),newDirection.z()),
                                           calcNorm,validNorm,&solNorm) ;
  if(calcNorm)
  { 
    G4Point3D newN = (*fDirectTransform3D)*G4Point3D(solNorm) ;
    newN.unit() ;
    *n = G4ThreeVector(newN.x(),newN.y(),newN.z()) ;
  }
  return dist ;  
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4ReflectedSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4Point3D newPoint = (*fDirectTransform3D)*G4Point3D(p) ;
  return fPtrSolid->DistanceToOut(
                    G4ThreeVector(newPoint.x(),newPoint.y(),newPoint.z())) ;   
}

//////////////////////////////////////////////////////////////
//
//

void 
G4ReflectedSolid::ComputeDimensions( G4VPVParameterisation* p,
	                                const G4int n,
                                        const G4VPhysicalVolume* pRep ) 
{
  // fPtrSolid->ComputeDimensions(p,n,pRep);

  G4Exception("ERROR: ComputeDimensions has no meaning for a G4ReflectedSolid. It cannot be called.");
}

/////////////////////////////////////////////////
//
//                    

void 
G4ReflectedSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddThis (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4ReflectedSolid::CreatePolyhedron () const 
{
  G4Polyhedron* polyhedron = fPtrSolid->CreatePolyhedron();
  polyhedron->Transform(*fDirectTransform3D);

  return polyhedron;
}

/////////////////////////////////////////////////////////
//
//

G4NURBS*      
G4ReflectedSolid::CreateNURBS      () const 
{
  // Take into account local transformation - see CreatePolyhedron.
  // return fPtrSolid->CreateNURBS() ;
  return 0;
}
