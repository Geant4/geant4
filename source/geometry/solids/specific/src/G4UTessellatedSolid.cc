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
// Implementation of G4UTessellatedSolid wrapper class
//
// 11.01.18 G.Cosmo, CERN
// --------------------------------------------------------------------

#include "G4TessellatedSolid.hh"
#include "G4UTessellatedSolid.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4TriangularFacet.hh"
#include "G4QuadrangularFacet.hh"

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

////////////////////////////////////////////////////////////////////////
//
// Constructors
//
G4UTessellatedSolid::G4UTessellatedSolid()
 : Base_t("")
{
}

G4UTessellatedSolid::G4UTessellatedSolid(const G4String& name)
 : Base_t(name)
{
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTessellatedSolid::G4UTessellatedSolid(__void__& a)
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTessellatedSolid::~G4UTessellatedSolid()
{
  std::size_t size = fFacets.size();
  for (std::size_t i = 0; i < size; ++i)  { delete fFacets[i]; }
  fFacets.clear();
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTessellatedSolid::G4UTessellatedSolid(const G4UTessellatedSolid& source)
  : Base_t(source)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTessellatedSolid&
G4UTessellatedSolid::operator=(const G4UTessellatedSolid& source)
{
  if (this == &source) return *this;
  
  Base_t::operator=( source );
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors

G4bool G4UTessellatedSolid::AddFacet(G4VFacet* aFacet)
{
  // Add a facet to the structure, checking validity.
  //
  if (GetSolidClosed())
  {
    G4Exception("G4UTessellatedSolid::AddFacet()", "GeomSolids1002",
                JustWarning, "Attempt to add facets when solid is closed.");
    return false;
  }
  if (!aFacet->IsDefined())
  {
    G4Exception("G4UTessellatedSolid::AddFacet()", "GeomSolids1002",
                JustWarning, "Attempt to add facet not properly defined.");    
    aFacet->StreamInfo(G4cout);
    return false;
  }
  if (aFacet->GetNumberOfVertices() == 3)
  {
    G4TriangularFacet* a3Facet = dynamic_cast<G4TriangularFacet*>(aFacet);
    return Base_t::AddTriangularFacet(U3Vector(a3Facet->GetVertex(0).x(),
                                               a3Facet->GetVertex(0).y(),
                                               a3Facet->GetVertex(0).z()),
                                      U3Vector(a3Facet->GetVertex(1).x(),
                                               a3Facet->GetVertex(1).y(),
                                               a3Facet->GetVertex(1).z()),
                                      U3Vector(a3Facet->GetVertex(2).x(),
                                               a3Facet->GetVertex(2).y(),
                                               a3Facet->GetVertex(2).z()),
                                      true);
  }
  else if (aFacet->GetNumberOfVertices() == 4)
  {
    G4QuadrangularFacet* a4Facet = dynamic_cast<G4QuadrangularFacet*>(aFacet);
    return Base_t::AddQuadrilateralFacet(U3Vector(a4Facet->GetVertex(0).x(),
                                                  a4Facet->GetVertex(0).y(),
                                                  a4Facet->GetVertex(0).z()),
                                         U3Vector(a4Facet->GetVertex(1).x(),
                                                  a4Facet->GetVertex(1).y(),
                                                  a4Facet->GetVertex(1).z()),
                                         U3Vector(a4Facet->GetVertex(2).x(),
                                                  a4Facet->GetVertex(2).y(),
                                                  a4Facet->GetVertex(2).z()),
                                         U3Vector(a4Facet->GetVertex(3).x(),
                                                  a4Facet->GetVertex(3).y(),
                                                  a4Facet->GetVertex(3).z()),
                                         true);
  }
  else
  {
    G4Exception("G4UTessellatedSolid::AddFacet()", "GeomSolids1002",
                JustWarning, "Attempt to add facet not properly defined.");    
    aFacet->StreamInfo(G4cout);
    return false;
  }
}

G4VFacet* G4UTessellatedSolid::GetFacet(G4int i) const
{
  return fFacets[i];
}

G4int G4UTessellatedSolid::GetNumberOfFacets() const
{
  return GetNFacets();
}

void G4UTessellatedSolid::SetSolidClosed(const G4bool t)
{
  if (t && !Base_t::IsClosed())
  {
    Base_t::Close();
    std::size_t nVertices = fTessellated.fVertices.size();
    std::size_t nFacets   = fTessellated.fFacets.size();
    for (std::size_t j = 0; j < nVertices; ++j)
    {
      U3Vector vt = fTessellated.fVertices[j];
      fVertexList.push_back(G4ThreeVector(vt.x(), vt.y(), vt.z()));
    }
    for (std::size_t i = 0; i < nFacets; ++i)
    {
      vecgeom::TriangleFacet<G4double>* afacet = Base_t::GetFacet(i);
      std::vector<G4ThreeVector> v;
      for (G4int k=0; k<3; ++k)
      {
        v.push_back(G4ThreeVector(afacet->fVertices[k].x(),
                                  afacet->fVertices[k].y(),
                                  afacet->fVertices[k].z()));
      }
      G4VFacet* facet = new G4TriangularFacet(v[0], v[1], v[2],
                                              G4FacetVertexType::ABSOLUTE);
      facet->SetVertices(&fVertexList);
      for (G4int k=0; k<3; ++k)
      {
        facet->SetVertexIndex(k, afacet->fIndices[k]);
      }
      fFacets.push_back(facet);
    }
  }
}

G4bool G4UTessellatedSolid::GetSolidClosed() const
{
  return Base_t::IsClosed();
}

void G4UTessellatedSolid::SetMaxVoxels(G4int)
{
  // Not yet implemented !
}

G4double G4UTessellatedSolid::GetMinXExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMin.x();
}
G4double G4UTessellatedSolid::GetMaxXExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMax.x();
}
G4double G4UTessellatedSolid::GetMinYExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMin.y();
}
G4double G4UTessellatedSolid::GetMaxYExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMax.y();
}
G4double G4UTessellatedSolid::GetMinZExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMin.z();
}
G4double G4UTessellatedSolid::GetMaxZExtent() const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  return aMax.z();
}

G4int G4UTessellatedSolid::AllocatedMemoryWithoutVoxels()
{
  G4int base = sizeof(*this);
  base += fVertexList.capacity() * sizeof(G4ThreeVector);

  std::size_t limit = fFacets.size();
  for (std::size_t i = 0; i < limit; ++i)
  {
    G4VFacet &facet = *fFacets[i];
    base += facet.AllocatedMemory();
  }
  return base;
}
G4int G4UTessellatedSolid::AllocatedMemory()
{
  return AllocatedMemoryWithoutVoxels();
}
void G4UTessellatedSolid::DisplayAllocatedMemory()
{
  G4int without = AllocatedMemoryWithoutVoxels();
  //  G4int with = AllocatedMemory();
  //  G4double ratio = (G4double) with / without;
  //  G4cout << "G4TessellatedSolid - Allocated memory without voxel overhead "
  //         << without << "; with " << with << "; ratio: " << ratio << G4endl; 
  G4cout << "G4TessellatedSolid - Allocated memory without voxel overhead "
         << without << G4endl; 
}


///////////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTessellatedSolid::BoundingLimits(G4ThreeVector& pMin,
                                         G4ThreeVector& pMax) const
{
  U3Vector aMin, aMax;
  Base_t::Extent(aMin, aMax);
  pMin = G4ThreeVector(aMin.x(), aMin.y(), aMin.z());
  pMax = G4ThreeVector(aMax.x(), aMax.y(), aMax.z());

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UTessellatedSolid::BoundingLimits()",
                "GeomMgt0001", JustWarning, message);
    StreamInfo(G4cout);
  }
}


//////////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTessellatedSolid::CalculateExtent(const EAxis pAxis,
                                     const G4VoxelLimits& pVoxelLimit,
                                     const G4AffineTransform& pTransform,
                                           G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);

  // Use simple bounding-box to help in the case of complex meshes
  //
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);

#if 0
  // Precise extent computation (disabled by default for this shape)
  //
  G4double kCarToleranceHalf = 0.5*kCarTolerance;
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return (pMin < pMax) ? true : false;
  }

  // The extent is calculated as cumulative extent of the pyramids
  // formed by facets and the center of the bounding box.
  //
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  G4ThreeVectorList base;
  G4ThreeVectorList apex(1);
  std::vector<const G4ThreeVectorList *> pyramid(2);
  pyramid[0] = &base;
  pyramid[1] = &apex;
  apex[0] = (bmin+bmax)*0.5;

  // main loop along facets
  pMin =  kInfinity;
  pMax = -kInfinity;
  for (G4int i=0; i<GetNumberOfFacets(); ++i)
  {
    G4VFacet* facet = GetFacet(i);
    if (std::abs((facet->GetSurfaceNormal()).dot(facet->GetVertex(0)-apex[0]))
        < kCarToleranceHalf) continue;

    base.resize(3);
    for (G4int k=0; k<3; ++k) { base[k] = facet->GetVertex(k); }
    G4double emin,emax;
    G4BoundingEnvelope benv(pyramid);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) break; // max possible extent
  }
  return (pMin < pMax);
#endif
}


///////////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron()
//
G4Polyhedron* G4UTessellatedSolid::CreatePolyhedron () const
{
  G4int nVertices = (G4int)fVertexList.size();
  G4int nFacets = (G4int)fFacets.size();
  G4Polyhedron* polyhedron = new G4Polyhedron(nVertices, nFacets);
  for (auto i = 0; i < nVertices; ++i)
  {
    polyhedron->SetVertex(i+1, fVertexList[i]);
  }

  for (auto i = 0; i < nFacets; ++i)
  {
    G4int v[3];  // Only facets with 3 vertices are defined in VecGeom
    G4VFacet* facet = GetFacet(i);
    for (auto j = 0; j < 3; ++j) // Retrieve indexing directly from VecGeom
    {
      v[j] = facet->GetVertexIndex(j) + 1;
    }
    polyhedron->SetFacet(i+1, v[0], v[1], v[2]);
  }
  polyhedron->SetReferences();  

  return polyhedron;
}

#endif  // G4GEOM_USE_USOLIDS
