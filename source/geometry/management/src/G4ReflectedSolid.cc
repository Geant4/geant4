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
// $Id: G4ReflectedSolid.cc 104317 2017-05-24 13:08:38Z gcosmo $
//
//
// Implementation for G4ReflectedSolid class
//
// Author: Vladimir Grichine, 23.07.01  (Vladimir.Grichine@cern.ch)
//
// --------------------------------------------------------------------

#include "G4ReflectedSolid.hh"

#include <sstream>

#include "G4Point3D.hh"
#include "G4Vector3D.hh"

#include "G4AffineTransform.hh"
#include "G4Transform3D.hh"
#include "G4VoxelLimits.hh"

#include "G4VPVParameterisation.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

/////////////////////////////////////////////////////////////////
//
// Constructor using HepTransform3D, in fact HepReflect3D

G4ReflectedSolid::G4ReflectedSolid( const G4String& pName,
                                          G4VSolid* pSolid ,
                                    const G4Transform3D& transform )
  : G4VSolid(pName), fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fPtrSolid = pSolid;
  fDirectTransform3D = new G4Transform3D(transform);
}

///////////////////////////////////////////////////////////////////
//

G4ReflectedSolid::~G4ReflectedSolid() 
{
  delete fDirectTransform3D; fDirectTransform3D=0;
  delete fpPolyhedron; fpPolyhedron = 0;
}

///////////////////////////////////////////////////////////////////
//

G4ReflectedSolid::G4ReflectedSolid(const G4ReflectedSolid& rhs)
  : G4VSolid(rhs), fPtrSolid(rhs.fPtrSolid),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fDirectTransform3D = new G4Transform3D(*rhs.fDirectTransform3D);
}

///////////////////////////////////////////////////////////////////
//

G4ReflectedSolid& G4ReflectedSolid::operator=(const G4ReflectedSolid& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)  { return *this; }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fPtrSolid= rhs.fPtrSolid;
  delete fDirectTransform3D;
  fDirectTransform3D= new G4Transform3D(*rhs.fDirectTransform3D);
  fRebuildPolyhedron = false;
  delete fpPolyhedron; fpPolyhedron= 0;

  return *this;
}

///////////////////////////////////////////////////////////////////
//

G4GeometryType G4ReflectedSolid::GetEntityType() const 
{
  return G4String("G4ReflectedSolid");
}

const G4ReflectedSolid* G4ReflectedSolid::GetReflectedSolidPtr() const   
{
  return this;
}

G4ReflectedSolid* G4ReflectedSolid::GetReflectedSolidPtr() 
{
  return this;
}

G4VSolid* G4ReflectedSolid::GetConstituentMovedSolid() const
{ 
  return fPtrSolid; 
} 

/////////////////////////////////////////////////////////////////////////////
//

G4Transform3D  G4ReflectedSolid::GetTransform3D() const
{
  return fDirectTransform3D->inverse();
}

G4Transform3D  G4ReflectedSolid::GetDirectTransform3D() const
{
  G4Transform3D aTransform= *fDirectTransform3D;
  return aTransform;
}

void G4ReflectedSolid::SetDirectTransform3D(G4Transform3D& transform) 
{
  fDirectTransform3D = &transform;
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4ReflectedSolid::BoundingLimits(G4ThreeVector& pMin,
                                      G4ThreeVector& pMax) const
{
  fPtrSolid->BoundingLimits(pMin,pMax);
  G4double xmin = pMin.x(), ymin = pMin.y(), zmin = pMin.z();
  G4double xmax = pMax.x(), ymax = pMax.y(), zmax = pMax.z();
  G4double xx = fDirectTransform3D->xx();
  G4double yy = fDirectTransform3D->yy();
  G4double zz = fDirectTransform3D->zz();

  if (std::abs(xx) == 1 && std::abs(yy) == 1 && std::abs(zz) == 1)
  {
    // Special case of reflection in axis and pure translation
    //
    if (xx == -1) { G4double tmp = -xmin; xmin = -xmax; xmax = tmp; }
    if (yy == -1) { G4double tmp = -ymin; ymin = -ymax; ymax = tmp; }
    if (zz == -1) { G4double tmp = -zmin; zmin = -zmax; zmax = tmp; }
    xmin += fDirectTransform3D->dx();
    xmax += fDirectTransform3D->dx();
    ymin += fDirectTransform3D->dy();
    ymax += fDirectTransform3D->dy();
    zmin += fDirectTransform3D->dz();
    zmax += fDirectTransform3D->dz();
  }
  else
  {
    // Use additional reflection in Z to set up affine transformation
    //
    G4Transform3D transform3D = G4ReflectZ3D()*(*fDirectTransform3D);
    G4AffineTransform transform(transform3D.getRotation().inverse(),
                                transform3D.getTranslation());
  
    // Find bounding box
    //
    G4VoxelLimits unLimit;
    fPtrSolid->CalculateExtent(kXAxis,unLimit,transform,xmin,xmax);
    fPtrSolid->CalculateExtent(kYAxis,unLimit,transform,ymin,ymax);
    fPtrSolid->CalculateExtent(kZAxis,unLimit,transform,zmin,zmax); 
  }

  pMin.set(xmin,ymin,-zmax);
  pMax.set(xmax,ymax,-zmin);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4ReflectedSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool 
G4ReflectedSolid::CalculateExtent( const EAxis pAxis,
                                   const G4VoxelLimits& pVoxelLimits,
                                   const G4AffineTransform& pTransform,
                                         G4double& pMin, 
                                         G4double& pMax ) const 
{
  // Separation of transformations. Calculation of the extent is done
  // in a reflection of the global space. In such way, the voxel is
  // reflected, but the solid is transformed just by G4AffineTransform.
  // It allows to use CalculateExtent() of the solid.

  // Reflect voxel limits in Z
  //
  G4VoxelLimits limits;
  limits.AddLimit(kXAxis, pVoxelLimits.GetMinXExtent(),
                          pVoxelLimits.GetMaxXExtent());
  limits.AddLimit(kYAxis, pVoxelLimits.GetMinYExtent(),
                          pVoxelLimits.GetMaxYExtent());
  limits.AddLimit(kZAxis,-pVoxelLimits.GetMaxZExtent(),
                         -pVoxelLimits.GetMinZExtent());

  // Set affine transformation
  //
  G4Transform3D transform3D = G4ReflectZ3D()*pTransform*(*fDirectTransform3D);
  G4AffineTransform transform(transform3D.getRotation().inverse(),
                              transform3D.getTranslation());

  // Find extent
  //
  if (!fPtrSolid->CalculateExtent(pAxis, limits, transform, pMin, pMax))
  {
    return false;
  }
  if (pAxis == kZAxis)
  {
    G4double tmp= -pMin; pMin= -pMax; pMax= tmp;
  }

  return true;
}

//////////////////////////////////////////////////////////////
//
// 

EInside G4ReflectedSolid::Inside(const G4ThreeVector& p ) const
{
  G4ThreeVector newPoint = (*fDirectTransform3D)*G4Point3D(p);
  return fPtrSolid->Inside(newPoint); 
}

//////////////////////////////////////////////////////////////
//
//

G4ThreeVector 
G4ReflectedSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  G4ThreeVector newPoint = (*fDirectTransform3D)*G4Point3D(p);
  G4Vector3D normal = fPtrSolid->SurfaceNormal(newPoint);
  return (*fDirectTransform3D)*normal;    
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)

G4double 
G4ReflectedSolid::DistanceToIn( const G4ThreeVector& p,
                                const G4ThreeVector& v ) const 
{    
  G4ThreeVector newPoint     = (*fDirectTransform3D)*G4Point3D(p);
  G4ThreeVector newDirection = (*fDirectTransform3D)*G4Vector3D(v);
  return fPtrSolid->DistanceToIn(newPoint,newDirection);   
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the intersection of
// two solids

G4double 
G4ReflectedSolid::DistanceToIn( const G4ThreeVector& p ) const 
{
  G4ThreeVector newPoint = (*fDirectTransform3D)*G4Point3D(p);
  return fPtrSolid->DistanceToIn(newPoint);   
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)

G4double 
G4ReflectedSolid::DistanceToOut( const G4ThreeVector& p,
                                 const G4ThreeVector& v,
                                 const G4bool calcNorm,
                                       G4bool *validNorm,
                                       G4ThreeVector *n ) const 
{
  G4ThreeVector solNorm; 

  G4ThreeVector newPoint     = (*fDirectTransform3D)*G4Point3D(p);
  G4ThreeVector newDirection = (*fDirectTransform3D)*G4Vector3D(v);

  G4double dist = fPtrSolid->DistanceToOut(newPoint, newDirection,
                                           calcNorm, validNorm, &solNorm);
  if(calcNorm)
  { 
    *n = (*fDirectTransform3D)*G4Vector3D(solNorm);
  }
  return dist;  
}

//////////////////////////////////////////////////////////////
//
// Inverted algorithm of DistanceToIn(p)

G4double 
G4ReflectedSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  G4ThreeVector newPoint = (*fDirectTransform3D)*G4Point3D(p);
  return fPtrSolid->DistanceToOut(newPoint);   
}

//////////////////////////////////////////////////////////////
//
//

void 
G4ReflectedSolid::ComputeDimensions(       G4VPVParameterisation*,
                                     const G4int,
                                     const G4VPhysicalVolume* ) 
{
  DumpInfo();
  G4Exception("G4ReflectedSolid::ComputeDimensions()",
               "GeomMgt0001", FatalException,
               "Method not applicable in this context!");
}

//////////////////////////////////////////////////////////////
//
// Return a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface

G4ThreeVector G4ReflectedSolid::GetPointOnSurface() const
{
  G4ThreeVector p        =  fPtrSolid->GetPointOnSurface();
  return (*fDirectTransform3D)*G4Point3D(p);
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of this object

G4VSolid* G4ReflectedSolid::Clone() const
{
  return new G4ReflectedSolid(*this);
}


//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4ReflectedSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for Reflected solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters of constituent solid: \n"
     << "===========================================================\n";
  fPtrSolid->StreamInfo(os);
  os << "===========================================================\n"
     << " Transformations: \n"
     << "    Direct transformation - translation : \n"
     << "           " << fDirectTransform3D->getTranslation() << "\n"
     << "                          - rotation    : \n"
     << "           ";
  fDirectTransform3D->getRotation().print(os);
  os << "\n"
     << "===========================================================\n";

  return os;
}

/////////////////////////////////////////////////
//
//                    

void 
G4ReflectedSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////
//
//

G4Polyhedron* 
G4ReflectedSolid::CreatePolyhedron () const 
{
  G4Polyhedron* polyhedron = fPtrSolid->CreatePolyhedron();
  if (polyhedron)
  {
    polyhedron->Transform(*fDirectTransform3D);
    return polyhedron;
  }
  else
  {
    std::ostringstream message;
    message << "Solid - " << GetName()
            << " - original solid has no" << G4endl
            << "corresponding polyhedron. Returning NULL!";
    G4Exception("G4ReflectedSolid::CreatePolyhedron()",
                "GeomMgt1001", JustWarning, message);
    return 0;
  }
}

/////////////////////////////////////////////////////////
//
//

G4Polyhedron*
G4ReflectedSolid::GetPolyhedron () const
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
