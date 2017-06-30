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
// GEANT4 tag $Name:$
//
//
// G4USolid implementation
//
// --------------------------------------------------------------------

#include "G4USolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VoxelLimits.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"
#include "G4PhysicalConstants.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

G4USolid::G4USolid(const G4String& name, VUSolid* s) :
  G4VSolid(name), fShape(s), fRebuildPolyhedron(false), fPolyhedron(0)
{
}

G4USolid::G4USolid(__void__& a)
  : G4VSolid(a), fShape(0), fRebuildPolyhedron(false), fPolyhedron(0)
{
}

G4USolid::~G4USolid()
{
  delete fPolyhedron; fPolyhedron = 0;
}

G4bool G4USolid::operator==(const G4USolid& s) const
{
  return (this == &s) ? true : false;
}

void G4USolid::ComputeDimensions(G4VPVParameterisation*,
                                 const G4int,
                                 const G4VPhysicalVolume*)
{
    std::ostringstream message;
    message << "Illegal call to G4USolid::ComputeDimensions()" << G4endl
            << "Method not overloaded by derived class !";
    G4Exception("G4USolid::ComputeDimensions()", "GeomSolids0003",
                FatalException, message);
}

void G4USolid::DescribeYourselfTo(G4VGraphicsScene& scene) const
{
  scene.AddSolid(*this);
}

G4GeometryType G4USolid::GetEntityType() const
{

  G4String string = fShape->GetEntityType();
  return "G4" + string;
}

std::ostream& G4USolid::StreamInfo(std::ostream& os) const
{
  return fShape->StreamInfo(os);
}

G4USolid::G4USolid(const G4USolid& rhs)
  : G4VSolid(rhs), fRebuildPolyhedron(false), fPolyhedron(0)
{
  fShape = rhs.fShape->Clone();
}

G4USolid& G4USolid::operator=(const G4USolid& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  G4VSolid::operator=(rhs);

  // Copy data
  //
  fShape = rhs.fShape->Clone();
  fRebuildPolyhedron = false;
  delete fPolyhedron; fPolyhedron = 0;

  return *this;
}

G4VSolid* G4USolid::Clone() const
{
  std::ostringstream message;
  message << "Clone() method not implemented for type: "
          << GetEntityType() << "!" << G4endl
          << "Returning NULL pointer!";
  G4Exception("G4USolid::Clone()", "GeomSolids1001", JustWarning, message);
  return 0;
}

G4bool G4USolid::CalculateExtent(const EAxis pAxis,
                                 const G4VoxelLimits& pVoxelLimit,
                                 const G4AffineTransform& pTransform,
                                 G4double& pMin, G4double& pMax) const
{
  UVector3 vmin, vmax;
  fShape->Extent(vmin,vmax);
  G4ThreeVector bmin(vmin.x(),vmin.y(),vmin.z());
  G4ThreeVector bmax(vmax.x(),vmax.y(),vmax.z());

  // Check correctness of the bounding box
  //
  if (bmin.x() >= bmax.x() || bmin.y() >= bmax.y() || bmin.z() >= bmax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " - " << GetEntityType() << " !"
            << "\nmin = " << bmin
            << "\nmax = " << bmax;
    G4Exception("G4USolid::CalculateExtent()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

G4Polyhedron* G4USolid::CreatePolyhedron() const
{
  // Must be implemented in concrete wrappers...

  std::ostringstream message;
  message << "Visualization not supported for USolid shape "
          << GetEntityType() << "... Sorry!" << G4endl;
  G4Exception("G4USolid::CreatePolyhedron()", "GeomSolids0003",
              FatalException, message);
  return 0;
}

G4Polyhedron* G4USolid::GetPolyhedron() const
{
  if (!fPolyhedron ||
      fRebuildPolyhedron ||
      fPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&polyhedronMutex);
    delete fPolyhedron;
    fPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fPolyhedron;
}

G4VisExtent G4USolid::GetExtent() const
{
  UVector3 vmin, vmax;
  fShape->Extent(vmin,vmax);
  return G4VisExtent(vmin.x(),vmax.x(),
                     vmin.y(),vmax.y(),
                     vmin.z(),vmax.z());
}

#endif  // G4GEOM_USE_USOLIDS
