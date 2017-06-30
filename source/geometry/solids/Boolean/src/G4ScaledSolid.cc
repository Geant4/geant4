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
// $Id:$
//
// Implementation for G4ScaledSolid class
//
// History:
//
// 27.10.15 G.Cosmo: created, based on implementation also provided in Root
//
// --------------------------------------------------------------------

#include "G4ScaledSolid.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "G4ScaleTransform.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

///////////////////////////////////////////////////////////////////
//
// Constructor
//
G4ScaledSolid::G4ScaledSolid( const G4String& pName,
                                    G4VSolid* pSolid ,
                              const G4Scale3D& pScale  )
  : G4VSolid(pName), fPtrSolid(pSolid),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fScale = new G4ScaleTransform(pScale);
}

///////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4ScaledSolid::G4ScaledSolid( __void__& a )
  : G4VSolid(a), fPtrSolid(0), fScale(0),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
}

///////////////////////////////////////////////////////////////////
//
// Destructor
//
G4ScaledSolid::~G4ScaledSolid() 
{
  delete fpPolyhedron; fpPolyhedron= 0;
  delete fScale; fScale= 0;
}

///////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4ScaledSolid::G4ScaledSolid(const G4ScaledSolid& rhs)
  : G4VSolid (rhs), fPtrSolid(rhs.fPtrSolid),
    fRebuildPolyhedron(false), fpPolyhedron(0)
{
  fScale = new G4ScaleTransform(*(rhs.fScale));
}

///////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4ScaledSolid& G4ScaledSolid::operator = (const G4ScaledSolid& rhs) 
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
  delete fScale;
  fScale = new G4ScaleTransform(*(rhs.fScale));
  fRebuildPolyhedron = false;
  delete fpPolyhedron; fpPolyhedron= 0;

  return *this;
}  

//////////////////////////////////////////////////////////////////////////
//
// Return original solid not scaled
//
G4VSolid* G4ScaledSolid::GetUnscaledSolid() const
{ 
  return fPtrSolid; 
} 

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4ScaledSolid::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
  G4ThreeVector bmin,bmax;
  G4ThreeVector scale = fScale->GetScale();
 
  fPtrSolid->BoundingLimits(bmin,bmax);
  pMin.set(bmin.x()*scale.x(),bmin.y()*scale.y(),bmin.z()*scale.z());
  pMax.set(bmax.x()*scale.x(),bmax.y()*scale.y(),bmax.z()*scale.z());

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
           << "\npMax = " << pMax;
    G4Exception("G4ScaledSolid::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
//
G4bool 
G4ScaledSolid::CalculateExtent( const EAxis pAxis,
                                const G4VoxelLimits& pVoxelLimit,
                                const G4AffineTransform& pTransform,
                                      G4double& pMin,
                                      G4double& pMax ) const
{
  // Find bounding box of unscaled solid
  G4ThreeVector bmin,bmax;
  fPtrSolid->BoundingLimits(bmin,bmax);

  // Set combined transformation
  G4Transform3D transform3D =
    G4Transform3D(pTransform.NetRotation().inverse(),
                  pTransform.NetTranslation())*GetScaleTransform();

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,transform3D,pMin,pMax);
}
 
/////////////////////////////////////////////////////
//
// Inside
//
EInside G4ScaledSolid::Inside(const G4ThreeVector& p) const
{
  return fPtrSolid->Inside(fScale->Transform(p));
}

//////////////////////////////////////////////////////////////
//
// SurfaceNormal
//
G4ThreeVector 
G4ScaledSolid::SurfaceNormal( const G4ThreeVector& p ) const 
{
  // Transform point to unscaled shape frame
  G4ThreeVector newPoint;
  fScale->Transform(p, newPoint);

  // Compute normal in unscaled frame
  G4ThreeVector newNormal = fPtrSolid->SurfaceNormal(newPoint); 
  G4ThreeVector normal;

  // Convert normal to scaled frame
  fScale->InverseTransformNormal(newNormal, normal);
  return normal/normal.mag();
}

/////////////////////////////////////////////////////////////
//
// The same algorithm as in DistanceToIn(p)
//
G4double 
G4ScaledSolid::DistanceToIn( const G4ThreeVector& p,
                             const G4ThreeVector& v  ) const 
{    
  // Transform point and direction to unscaled shape frame
  G4ThreeVector newPoint;
  fScale->Transform(p, newPoint);

  // Direction is un-normalized after scale transformation
  G4ThreeVector newDirection;
  fScale->Transform(v, newDirection);
  newDirection = newDirection/newDirection.mag();

  // Compute distance in unscaled system
  G4double dist = fPtrSolid->DistanceToIn(newPoint,newDirection);

  // Return converted distance to global
  return fScale->InverseTransformDistance(dist, newDirection);
}

////////////////////////////////////////////////////////
//
// Approximate nearest distance from the point p to the solid from outside
//
G4double 
G4ScaledSolid::DistanceToIn( const G4ThreeVector& p ) const 
{
  // Transform point to unscaled shape frame
  G4ThreeVector newPoint;
  fScale->Transform(p, newPoint);

  // Compute unscaled safety, then scale it.
  G4double dist = fPtrSolid->DistanceToIn(newPoint);
  return fScale->InverseTransformDistance(dist);
}

//////////////////////////////////////////////////////////
//
// The same algorithm as DistanceToOut(p)
//
G4double 
G4ScaledSolid::DistanceToOut( const G4ThreeVector& p,
                              const G4ThreeVector& v,
                              const G4bool calcNorm,
                                    G4bool *validNorm,
                                    G4ThreeVector *n   ) const 
{
  // Transform point and direction to unscaled shape frame
  G4ThreeVector newPoint;
  fScale->Transform(p, newPoint);

  // Direction is un-normalized after scale transformation
  G4ThreeVector newDirection;
  fScale->Transform(v, newDirection);
  newDirection = newDirection/newDirection.mag();

  // Compute distance in unscaled system
  G4ThreeVector solNorm; 
  G4double dist = fPtrSolid->DistanceToOut(newPoint,newDirection,
                                           calcNorm,validNorm,&solNorm);
  if(calcNorm)
  { 
    G4ThreeVector normal;
    fScale->TransformNormal(solNorm, normal);
    *n = normal/normal.mag();
  }

  // Return distance converted to global
  return fScale->InverseTransformDistance(dist, newDirection);
}

//////////////////////////////////////////////////////////////
//
//  Approximate nearest distance from the point p to the solid from inside
//
G4double 
G4ScaledSolid::DistanceToOut( const G4ThreeVector& p ) const 
{
  // Transform point to unscaled shape frame
  G4ThreeVector newPoint;
  fScale->Transform(p, newPoint);

  // Compute unscaled safety, then scale it.
  G4double dist = fPtrSolid->DistanceToOut(newPoint);
  return fScale->InverseTransformDistance(dist);
}

//////////////////////////////////////////////////////////////
//
// ComputeDimensions
//
void 
G4ScaledSolid::ComputeDimensions( G4VPVParameterisation*,
                                  const G4int,
                                  const G4VPhysicalVolume* ) 
{
  DumpInfo();
  G4Exception("G4ScaledSolid::ComputeDimensions()",
              "GeomSolids0001", FatalException,
              "Method not applicable in this context!");
}

//////////////////////////////////////////////////////////////////////////
//
// Returns a point (G4ThreeVector) randomly and uniformly selected
// on the solid surface
//
G4ThreeVector G4ScaledSolid::GetPointOnSurface() const
{
  return fScale->InverseTransform(fPtrSolid->GetPointOnSurface());
}

//////////////////////////////////////////////////////////////////////////
//
// Return object type name
//
G4GeometryType G4ScaledSolid::GetEntityType() const 
{
  return G4String("G4ScaledSolid");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4ScaledSolid::Clone() const
{
  return new G4ScaledSolid(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Returning the scaling transformation
//
G4Scale3D G4ScaledSolid::GetScaleTransform() const
{
  return G4Scale3D(fScale->GetScale().x(),
                   fScale->GetScale().y(),
                   fScale->GetScale().z());
}

//////////////////////////////////////////////////////////////////////////
//
// Setting the scaling transformation
//
void G4ScaledSolid::SetScaleTransform(const G4Scale3D& scale)
{
  if (fScale) { delete fScale; }
  fScale = new G4ScaleTransform(scale);
  fRebuildPolyhedron = true;
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream
//
std::ostream& G4ScaledSolid::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for Scaled solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters of constituent solid: \n"
     << "===========================================================\n";
  fPtrSolid->StreamInfo(os);
  os << "===========================================================\n"
     << " Scaling: \n"
     << "    Scale transformation : \n"
     << "           " << fScale->GetScale().x() << ", "
                      << fScale->GetScale().y() << ", "
                      << fScale->GetScale().z() << "\n"
     << "===========================================================\n";

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// DescribeYourselfTo
//
void 
G4ScaledSolid::DescribeYourselfTo ( G4VGraphicsScene& scene ) const 
{
  scene.AddSolid (*this);
}

//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* 
G4ScaledSolid::CreatePolyhedron () const 
{
  G4Polyhedron* polyhedron = fPtrSolid->CreatePolyhedron();
  polyhedron->Transform(GetScaleTransform());
  return polyhedron;
}

//////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4ScaledSolid::GetPolyhedron () const
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
