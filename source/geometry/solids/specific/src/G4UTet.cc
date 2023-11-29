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
// Implementation for G4UTet wrapper class
//
// 1.11.13 G.Cosmo, CERN
// --------------------------------------------------------------------

#include "G4Tet.hh"
#include "G4UTet.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// This class is implemented separately from general polyhedra,
// because the simplex geometry can be computed very quickly,
// which may become important in situations imported from mesh generators,
// in which a very large number of G4Tets are created.
// A Tet has all of its geometrical information precomputed
//
G4UTet::G4UTet(const G4String& pName,
               const G4ThreeVector& anchor,
               const G4ThreeVector& p1,
               const G4ThreeVector& p2,
               const G4ThreeVector& p3, G4bool* degeneracyFlag)
  : Base_t(pName, U3Vector(anchor.x(),anchor.y(),anchor.z()),
                  U3Vector(p1.x(), p1.y(), p1.z()),
                  U3Vector(p2.x(), p2.y(), p2.z()),
                  U3Vector(p3.x(), p3.y(), p3.z()))
{
  // Check for degeneracy
  G4bool degenerate = CheckDegeneracy(anchor, p1, p2, p3);
  if(degeneracyFlag != nullptr) *degeneracyFlag = degenerate;
  else if (degenerate)
  {
    G4Exception("G4UTet::G4UTet()", "GeomSolids0002", FatalException,
                "Degenerate tetrahedron not allowed.");
  }

  // Set bounding box
  for (G4int i = 0; i < 3; ++i)
  {
    fBmin[i] = std::min(std::min(std::min(anchor[i], p1[i]), p2[i]), p3[i]);
    fBmax[i] = std::max(std::max(std::max(anchor[i], p1[i]), p2[i]), p3[i]);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UTet::G4UTet( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UTet::~G4UTet() = default;

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTet::G4UTet(const G4UTet& rhs)
  : Base_t(rhs)
{
  fBmin = rhs.fBmin;
  fBmax = rhs.fBmax;
}


///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTet& G4UTet::operator = (const G4UTet& rhs)
{
  // Check assignment to self
  if (this == &rhs)  { return *this; }

  // Copy base class data
  Base_t::operator=(rhs);

  // Copy bounding box
  fBmin = rhs.fBmin;
  fBmax = rhs.fBmax;

  return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// Return true if tetrahedron is degenerate
// Tetrahedron is concidered as degenerate in case if its minimal
// height is less than the degeneracy tolerance
//
G4bool G4UTet::CheckDegeneracy(const G4ThreeVector& p0,
                               const G4ThreeVector& p1,
                               const G4ThreeVector& p2,
                               const G4ThreeVector& p3) const
{
  G4double hmin = 4. * kCarTolerance; // degeneracy tolerance

  // Calculate volume
  G4double vol = std::abs((p1 - p0).cross(p2 - p0).dot(p3 - p0));

  // Calculate face areas squared
  G4double ss[4];
  ss[0] = ((p1 - p0).cross(p2 - p0)).mag2();
  ss[1] = ((p2 - p0).cross(p3 - p0)).mag2();
  ss[2] = ((p3 - p0).cross(p1 - p0)).mag2();
  ss[3] = ((p2 - p1).cross(p3 - p1)).mag2();

  // Find face with max area
  G4int k = 0;
  for (G4int i = 1; i < 4; ++i) { if (ss[i] > ss[k]) k = i; }

  // Check: vol^2 / s^2 <= hmin^2
  return (vol*vol <= ss[k]*hmin*hmin);
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UTet::ComputeDimensions(G4VPVParameterisation*,
                               const G4int,
                               const G4VPhysicalVolume*)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4UTet::Clone() const
{
  return new G4UTet(*this);
}

///////////////////////////////////////////////////////////////////////////////
//
// Modifier
//
void G4UTet::SetVertices(const G4ThreeVector& anchor,
                         const G4ThreeVector& p1,
                         const G4ThreeVector& p2,
                         const G4ThreeVector& p3,
                         G4bool* degeneracyFlag)
{
  // Check for degeneracy
  G4bool degenerate = CheckDegeneracy(anchor, p1, p2, p3);
  if(degeneracyFlag != nullptr) *degeneracyFlag = degenerate;
  else if (degenerate)
  {
    G4Exception("G4UTet::SetVertices()", "GeomSolids0002", FatalException,
                "Degenerate tetrahedron not allowed.");
  }

  // Change tetrahedron
  *this = G4UTet(GetName(), anchor, p1, p2, p3, &degenerate);
}

///////////////////////////////////////////////////////////////////////////////
//
// Accessors
//
void G4UTet::GetVertices(G4ThreeVector& anchor,
                         G4ThreeVector& p1,
                         G4ThreeVector& p2,
                         G4ThreeVector& p3) const
{
  std::vector<U3Vector> vec(4);
  Base_t::GetVertices(vec[0], vec[1], vec[2], vec[3]);
  anchor = G4ThreeVector(vec[0].x(), vec[0].y(), vec[0].z());
  p1 = G4ThreeVector(vec[1].x(), vec[1].y(), vec[1].z());
  p2 = G4ThreeVector(vec[2].x(), vec[2].y(), vec[2].z());
  p3 = G4ThreeVector(vec[3].x(), vec[3].y(), vec[3].z());
}

std::vector<G4ThreeVector> G4UTet::GetVertices() const
{
  std::vector<U3Vector> vec(4);
  Base_t::GetVertices(vec[0], vec[1], vec[2], vec[3]);
  std::vector<G4ThreeVector> vertices;
  for (unsigned int i=0; i<4; ++i)
  {
    G4ThreeVector v(vec[i].x(), vec[i].y(), vec[i].z());
    vertices.push_back(v);
  }
  return vertices;
}

////////////////////////////////////////////////////////////////////////
//
// Set bounding box
//
void G4UTet::SetBoundingLimits(const G4ThreeVector& pMin,
                               const G4ThreeVector& pMax)
{
  G4ThreeVector fVertex[4];
  GetVertices(fVertex[0], fVertex[1], fVertex[2], fVertex[3]);

  G4int iout[4] = { 0, 0, 0, 0 };
  for (G4int i = 0; i < 4; ++i)
  {
    iout[i] = (G4int)(fVertex[i].x() < pMin.x() ||
                      fVertex[i].y() < pMin.y() ||
                      fVertex[i].z() < pMin.z() ||
                      fVertex[i].x() > pMax.x() ||
                      fVertex[i].y() > pMax.y() ||
                      fVertex[i].z() > pMax.z());
  }
  if (iout[0] + iout[1] + iout[2] + iout[3] != 0)
  {
    std::ostringstream message;
    message << "Attempt to set bounding box that does not encapsulate solid: "
            << GetName() << " !\n"
            << "  Specified bounding box limits:\n"
            << "    pmin: " << pMin << "\n"
            << "    pmax: " << pMax << "\n"
            << "  Tetrahedron vertices:\n"
            << "    anchor " << fVertex[0] << ((iout[0]) != 0 ? " is outside\n" : "\n")
            << "    p1 "     << fVertex[1] << ((iout[1]) != 0 ? " is outside\n" : "\n")
            << "    p2 "     << fVertex[2] << ((iout[2]) != 0 ? " is outside\n" : "\n")
            << "    p3 "     << fVertex[3] << ((iout[3]) != 0 ? " is outside"   : "");
    G4Exception("G4UTet::SetBoundingLimits()", "GeomSolids0002",
                FatalException, message);
  }
  fBmin = pMin;
  fBmax = pMax;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTet::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin = fBmin;
  pMax = fBmax;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UTet::CalculateExtent(const EAxis pAxis,
                        const G4VoxelLimits& pVoxelLimit,
                        const G4AffineTransform& pTransform,
                              G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);

  // Use simple bounding-box to help in the case of complex 3D meshes
  //
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);

#if 0
  // Precise extent computation (disabled by default for this shape)
  //
  G4bool exist;
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Set bounding envelope (benv) and calculate extent
  //
  std::vector<G4ThreeVector> vec = GetVertices();

  G4ThreeVectorList anchor(1);
  anchor[0] = vec[0];

  G4ThreeVectorList base(3);
  base[0] = vec[1];
  base[1] = vec[2];
  base[2] = vec[3];

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &anchor;
  polygons[1] = &base;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  return exists = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UTet::CreatePolyhedron() const
{
  std::vector<U3Vector> vec(4);
  Base_t::GetVertices(vec[0], vec[1], vec[2], vec[3]);

  G4double xyz[4][3];
  const G4int faces[4][4] = {{1,3,2,0},{1,4,3,0},{1,2,4,0},{2,3,4,0}};
  for (unsigned int i=0; i<4; ++i)
  {
    xyz[i][0] = vec[i].x();
    xyz[i][1] = vec[i].y();
    xyz[i][2] = vec[i].z();
  }

  auto ph = new G4Polyhedron;
  ph->createPolyhedron(4,4,xyz,faces);
  return ph;
}

#endif  // G4GEOM_USE_USOLIDS
