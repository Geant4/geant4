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
// $Id: G4DisplacedSolid.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// Implementation for G4DisplacedSolid class for boolean 
// operations between other solids
//
// History:
//
// 28.10.98 V.Grichine: created
// 14.11.99 V.Grichine: modifications in CalculateExtent(...) method
// 22.11.00 V.Grichine: new set methods for matrix/vectors
//
// --------------------------------------------------------------------

#include "G4DisplacedSolid.hh"

#include "G4VoxelLimits.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

////////////////////////////////////////////////////////////////
//
// Constructor for transformation like rotation of frame then translation 
// in new frame. It is similar to 1st constractor in G4PVPlacement

G4DisplacedSolid::G4DisplacedSolid( const G4String& pName,
                                          G4VSolid* pSolid ,
                                          G4RotationMatrix* rotMatrix,
                                    const G4ThreeVector& transVector    )
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fPtrSolid = pSolid ;
  fPtrTransform = new G4AffineTransform(rotMatrix,transVector) ;
  fPtrTransform->Invert() ;
  fDirectTransform = new G4AffineTransform(rotMatrix,transVector) ;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Constructor

G4DisplacedSolid::G4DisplacedSolid( const G4String& pName,
                                          G4VSolid* pSolid ,
                                    const G4Transform3D& transform  )
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fPtrSolid = pSolid ;
  fDirectTransform = new G4AffineTransform(transform.getRotation().inverse(),
                                           transform.getTranslation()) ;

  fPtrTransform    = new G4AffineTransform(transform.getRotation().inverse(),
                                           transform.getTranslation()) ;
  fPtrTransform->Invert() ;
}

///////////////////////////////////////////////////////////////////
//
// Constructor for use with creation of Transient object
// from Persistent object

G4DisplacedSolid::G4DisplacedSolid( const G4String& pName,
                                          G4VSolid* pSolid ,
                                    const G4AffineTransform directTransform )
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fPtrSolid = pSolid ;
  fDirectTransform = new G4AffineTransform( directTransform );
  fPtrTransform    = new G4AffineTransform( directTransform.Inverse() ) ; 
}

///////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.

G4DisplacedSolid::G4DisplacedSolid( __void__& a )
  : G4VSolid(a), fPtrSolid(0), fPtrTransform(0),
    fDirectTransform(0), fRebuildPolyhedron(false), fpPolyhedron(0)
{
}

///////////////////////////////////////////////////////////////////
//
// Destructor

G4DisplacedSolid::~G4DisplacedSolid() 
{
  CleanTransformations();
  delete fpPolyhedron; fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////
//
// Copy constructor

G4DisplacedSolid::G4DisplacedSolid(const G4DisplacedSolid& rhs)
  : G4VSolid (rhs), fPtrSolid(rhs.fPtrSolid),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fPtrTransform = new G4AffineTransform(*(rhs.fPtrTransform));
  fDirectTransform = new G4AffineTransform(*(rhs.fDirectTransform));
}

///////////////////////////////////////////////////////////////
//
// Assignment operator

G4DisplacedSolid& G4DisplacedSolid::operator = (const G4DisplacedSolid& rhs) 
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fPtrSolid = rhs.fPtrSolid;
  delete fPtrTransform; delete fDirectTransform;
  fPtrTransform = new G4AffineTransform(*(rhs.fPtrTransform));
  fDirectTransform = new G4AffineTransform(*(rhs.fDirectTransform));
  fRebuildPolyhedron = false;
  delete fpPolyhedron; fpPolyhedron= 0;

  return *this;
}  

void G4DisplacedSolid::CleanTransformations()
{
  if(fPtrTransform)
  {
    delete fPtrTransform;  fPtrTransform=0;
    delete fDirectTransform;  fDirectTransform=0;
  }
}

const G4DisplacedSolid* G4DisplacedSolid::GetDisplacedSolidPtr() const   
{
  return this;
}

G4DisplacedSolid* G4DisplacedSolid::GetDisplacedSolidPtr() 
{
  return this;
}

G4VSolid* G4DisplacedSolid::GetConstituentMovedSolid() const
{ 
  return fPtrSolid; 
} 

/////////////////////////////////////////////////////////////////////////////

G4AffineTransform  G4DisplacedSolid::GetTransform() const
{
  G4AffineTransform aTransform = *fPtrTransform;
  return aTransform;
}

void G4DisplacedSolid::SetTransform(G4AffineTransform& transform) 
{
  fPtrTransform = &transform ;
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////////

G4AffineTransform  G4DisplacedSolid::GetDirectTransform() const
{
  G4AffineTransform aTransform= *fDirectTransform;
  return aTransform;
}

void G4DisplacedSolid::SetDirectTransform(G4AffineTransform& transform) 
{
  fDirectTransform = &transform ;
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////////

G4RotationMatrix G4DisplacedSolid::GetFrameRotation() const
{
  G4RotationMatrix InvRotation= fDirectTransform->NetRotation();
  return InvRotation;
}

void G4DisplacedSolid::SetFrameRotation(const G4RotationMatrix& matrix)
{
  fDirectTransform->SetNetRotation(matrix);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////////

G4ThreeVector  G4DisplacedSolid::GetFrameTranslation() const
{
  return fPtrTransform->NetTranslation();
}

void G4DisplacedSolid::SetFrameTranslation(const G4ThreeVector& vector)
{
  fPtrTransform->SetNetTranslation(vector);
  fRebuildPolyhedron = true;
}

///////////////////////////////////////////////////////////////

G4RotationMatrix G4DisplacedSolid::GetObjectRotation() const
{
  G4RotationMatrix Rotation= fPtrTransform->NetRotation();
  return Rotation;
}

void G4DisplacedSolid::SetObjectRotation(const G4RotationMatrix& matrix)
{
  fPtrTransform->SetNetRotation(matrix);
  fRebuildPolyhedron = true;
}

///////////////////////////////////////////////////////////////////////

G4ThreeVector  G4DisplacedSolid::GetObjectTranslation() const
{
  return fDirectTransform->NetTranslation();
}

void G4DisplacedSolid::SetObjectTranslation(const G4ThreeVector& vector)
{
  fDirectTransform->SetNetTranslation(vector);
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4DisplacedSolid::BoundingLimits(G4ThreeVector& pMin,
                                      G4ThreeVector& pMax) const
{
  if (!fDirectTransform->IsRotated())
  {
    // Special case of pure translation
    //
    fPtrSolid->BoundingLimits(pMin,pMax);
    G4ThreeVector offset = fDirectTransform->NetTranslation();
    pMin += offset;
    pMax += offset;
  }
  else
  {
    // General case, use CalculateExtent() to find bounding box
    //
    G4VoxelLimits unLimit;
    G4double xmin,xmax,ymin,ymax,zmin,zmax;
    fPtrSolid->CalculateExtent(kXAxis,unLimit,*fDirectTransform,xmin,xmax);
    fPtrSolid->CalculateExtent(kYAxis,unLimit,*fDirectTransform,ymin,ymax);
    fPtrSolid->CalculateExtent(kZAxis,unLimit,*fDirectTransform,zmin,zmax);
    pMin.set(xmin,ymin,zmin);
    pMax.set(xmax,ymax,zmax);
  }
  
  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4DisplacedSolid::BoundingLimits()", "GeomMgt0001",
               JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
     
G4bool 
G4DisplacedSolid::CalculateExtent( const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, 
                                         G4double& pMax           ) const 
{
  G4AffineTransform sumTransform ;
  sumTransform.Product(*fDirectTransform,pTransform) ;
  return fPtrSolid->CalculateExtent(pAxis,pVoxelLimit,sumTransform,pMin,pMax) ;
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
G4DisplacedSolid::DistanceToIn( const G4ThreeVector& p ) const 
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
                                       G4ThreeVector *n   ) const 
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
G4DisplacedSolid::ComputeDimensions(       G4VPVParameterisation*,
                                     const G4int,
                                     const G4VPhysicalVolume* ) 
{
  DumpInfo();
  G4Exception("G4DisplacedSolid::ComputeDimensions()",
              "GeomSolids0001", FatalException,
              "Method not applicable in this context!");
}

//////////////////////////////////////////////////////////////////////////
//
// Returns a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface
//

G4ThreeVector G4DisplacedSolid::GetPointOnSurface() const
{
  G4ThreeVector p =  fPtrSolid->GetPointOnSurface();
  return fDirectTransform->TransformPoint(p);
}

//////////////////////////////////////////////////////////////////////////
//
// Return object type name

G4GeometryType G4DisplacedSolid::GetEntityType() const 
{
  return G4String("G4DisplacedSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4DisplacedSolid::Clone() const
{
  return new G4DisplacedSolid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4DisplacedSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for Displaced solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters of constituent solid: \n"
     << "===========================================================\n";
  fPtrSolid->StreamInfo(os);
  os << "===========================================================\n"
     << " Transformations: \n"
     << "    Direct transformation - translation : \n"
     << "           " << fDirectTransform->NetTranslation() << "\n"
     << "                          - rotation    : \n"
     << "           ";
  fDirectTransform->NetRotation().print(os);
  os << "\n"
     << "===========================================================\n";

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
//                    

void 
G4DisplacedSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

//////////////////////////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4DisplacedSolid::CreatePolyhedron () const 
{
  G4Polyhedron* polyhedron = fPtrSolid->CreatePolyhedron();
  polyhedron
    ->Transform(G4Transform3D(GetObjectRotation(),GetObjectTranslation()));
  return polyhedron;
}

//////////////////////////////////////////////////////////////////////////
//
//

G4Polyhedron* G4DisplacedSolid::GetPolyhedron () const
{
  if (!fpPolyhedron ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
    {
      fpPolyhedron = CreatePolyhedron();
      fRebuildPolyhedron = false;
    }
  return fpPolyhedron;
}
