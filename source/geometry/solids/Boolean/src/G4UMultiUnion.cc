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
// 
// Implementation of G4UMultiUnion wrapper class
// --------------------------------------------------------------------

#include "G4UMultiUnion.hh"

#if defined(G4GEOM_USE_USOLIDS)

#include "G4VoxelLimits.hh"
#include "G4BoundingEnvelope.hh"
#include "G4Polyhedron.hh"
#include "G4DisplacedSolid.hh"
#include "G4RotationMatrix.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UMultiUnion::G4UMultiUnion(const G4String& name)
  : G4USolid(name, new UMultiUnion(name))
{ 
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UMultiUnion::G4UMultiUnion(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UMultiUnion::~G4UMultiUnion()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UMultiUnion::G4UMultiUnion(const G4UMultiUnion &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UMultiUnion& G4UMultiUnion::operator=(const G4UMultiUnion &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}


//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers
//
void G4UMultiUnion::AddNode(G4VSolid& solid, G4Transform3D& trans)
{
  HepGeom::Rotate3D rot;
  HepGeom::Translate3D transl ;
  HepGeom::Scale3D scale;

  trans.getDecomposition(scale,rot,transl); 
  G4ThreeVector pos = transl.getTranslation();
    
  UTransform3D tr;
  tr.fRot[0] = rot.xx(); tr.fRot[1] = rot.xy(); tr.fRot[2] = rot.xz();
  tr.fRot[3] = rot.yx(); tr.fRot[4] = rot.yy(); tr.fRot[5] = rot.yz();
  tr.fRot[6] = rot.zx(); tr.fRot[7] = rot.zy(); tr.fRot[8] = rot.zz();
  tr.fTr = UVector3(pos.x(), pos.y(), pos.z());
 
  GetShape()->AddNode(*(static_cast<G4USolid&>(solid).GetSolid()), tr);
}

G4Transform3D* G4UMultiUnion::GetTransformation(G4int index) const
{
  UTransform3D tr = GetShape()->GetTransformation(index);

  G4RotationMatrix
    rot(CLHEP::HepRep3x3(tr.fRot[0], tr.fRot[1], tr.fRot[2],
                         tr.fRot[3], tr.fRot[4], tr.fRot[5],
                         tr.fRot[6], tr.fRot[7], tr.fRot[8]));
  G4ThreeVector transl(tr.fTr.x(), tr.fTr.y(), tr.fTr.z());

  return new G4Transform3D(rot, transl);
}

G4VSolid* G4UMultiUnion::GetSolid(G4int index) const
{
  VUSolid* solid = GetShape()->GetSolid(index);
  return new G4USolid(solid->GetName(), solid);
}

G4int G4UMultiUnion::GetNumberOfSolids()const
{
  return GetShape()->GetNumberOfSolids();
}

void G4UMultiUnion::Voxelize()
{
  GetShape()->Voxelize();
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UMultiUnion::Extent(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  UVector3 vmin, vmax;
  GetShape()->Extent(vmin,vmax);
  pMin.set(vmin.x(),vmin.y(),vmin.z());
  pMax.set(vmax.x(),vmax.y(),vmax.z());

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UMultiUnion::Extent()", "GeomMgt0001", JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UMultiUnion::CalculateExtent(const EAxis pAxis,
                               const G4VoxelLimits& pVoxelLimit,
                               const G4AffineTransform& pTransform,
                                     G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Get bounding box
  Extent(bmin,bmax);

  // Find extent
  G4BoundingEnvelope bbox(bmin,bmax);
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
}

//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UMultiUnion::CreatePolyhedron() const
{

  HepPolyhedronProcessor processor;
  HepPolyhedronProcessor::Operation operation = HepPolyhedronProcessor::UNION;

  G4VSolid* solidA = GetSolid(0);
  const G4Transform3D* transform0=GetTransformation(0);
  G4RotationMatrix rot0=(*transform0).getRotation();
  const  G4ThreeVector transl0 = (*transform0).getTranslation();
  G4DisplacedSolid dispSolidA("placedA",solidA,&rot0,transl0);
  delete transform0;

  G4Polyhedron* top = new G4Polyhedron(*dispSolidA.GetPolyhedron());
    
  for(G4int i=1; i<GetNumberOfSolids(); ++i)
  {
    G4VSolid* solidB = GetSolid(i);
    const G4Transform3D* transform=GetTransformation(i);
    G4RotationMatrix rot=(*transform).getRotation();
    const  G4ThreeVector transl = (*transform).getTranslation();
    G4DisplacedSolid dispSolidB("placedB",solidB,&rot,transl);
    G4Polyhedron* operand = dispSolidB.GetPolyhedron();
    processor.push_back (operation, *operand);
    delete transform;
  }
   
  if (processor.execute(*top)) { return top; }
  else { return 0; } 
}

#endif  // G4GEOM_USE_USOLIDS
