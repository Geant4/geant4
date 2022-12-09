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
// * This  code  implementation is the  intellectual property  of the *
// * Vanderbilt University Free Electron Laser Center                 *
// * Vanderbilt University, Nashville, TN, USA                        *
// * Development supported by:                                        *
// * United States MFEL program  under grant FA9550-04-1-0045         *
// * and NASA under contract number NNG04CT05P                        *
// * Written by Marcus H. Mendenhall and Robert A. Weller.            *
// *                                                                  *
// * Contributed to the Geant4 Core, January, 2005.                   *
// *                                                                  *
// ********************************************************************
//
// Implementation for G4Tet class
//
// 03.09.2004 - Marcus Mendenhall, created
// 08.01.2020 - Evgueni Tcherniaev, complete revision, speed up
// --------------------------------------------------------------------

#include "G4Tet.hh"

#if !defined(G4GEOM_USE_UTET)

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"

#include "G4QuickRand.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor - create a tetrahedron
// A Tet has all of its geometrical information precomputed
//
G4Tet::G4Tet(const G4String& pName,
             const G4ThreeVector& p0,
             const G4ThreeVector& p1,
             const G4ThreeVector& p2,
             const G4ThreeVector& p3, G4bool* degeneracyFlag)
  : G4VSolid(pName)
{
  // Check for degeneracy
  G4bool degenerate = CheckDegeneracy(p0, p1, p2, p3);
  if (degeneracyFlag)
  {
    *degeneracyFlag = degenerate;
  }
  else if (degenerate)
  {
    std::ostringstream message;
    message << "Degenerate tetrahedron: " << GetName() << " !\n"
            << "  anchor: " << p0 << "\n"
            << "  p1    : " << p1 << "\n"
            << "  p2    : " << p2 << "\n"
            << "  p3    : " << p3 << "\n"
            << "  volume: "
            << std::abs((p1 - p0).cross(p2 - p0).dot(p3 - p0))/6.;
    G4Exception("G4Tet::G4Tet()", "GeomSolids0002", FatalException, message);
  }

  // Define surface thickness
  halfTolerance = 0.5 * kCarTolerance;

  // Set data members
  Initialize(p0, p1, p2, p3);
}

////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Tet::G4Tet( __void__& a )
  : G4VSolid(a)
{
}

////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4Tet::~G4Tet()
{
  delete fpPolyhedron; fpPolyhedron = nullptr;
}

////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4Tet::G4Tet(const G4Tet& rhs)
  : G4VSolid(rhs)
{
   halfTolerance = rhs.halfTolerance;
   fCubicVolume = rhs.fCubicVolume;
   fSurfaceArea = rhs.fSurfaceArea;
   for (G4int i = 0; i < 4; ++i) { fVertex[i] = rhs.fVertex[i]; }
   for (G4int i = 0; i < 4; ++i) { fNormal[i] = rhs.fNormal[i]; }
   for (G4int i = 0; i < 4; ++i) { fDist[i] = rhs.fDist[i]; }
   for (G4int i = 0; i < 4; ++i) { fArea[i] = rhs.fArea[i]; }
   fBmin = rhs.fBmin;
   fBmax = rhs.fBmax;
}

////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4Tet& G4Tet::operator = (const G4Tet& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4VSolid::operator=(rhs);

   // Copy data
   //
   halfTolerance = rhs.halfTolerance;
   fCubicVolume = rhs.fCubicVolume;
   fSurfaceArea = rhs.fSurfaceArea;
   for (G4int i = 0; i < 4; ++i) { fVertex[i] = rhs.fVertex[i]; }
   for (G4int i = 0; i < 4; ++i) { fNormal[i] = rhs.fNormal[i]; }
   for (G4int i = 0; i < 4; ++i) { fDist[i] = rhs.fDist[i]; }
   for (G4int i = 0; i < 4; ++i) { fArea[i] = rhs.fArea[i]; }
   fBmin = rhs.fBmin;
   fBmax = rhs.fBmax;
   fRebuildPolyhedron = false;
   delete fpPolyhedron; fpPolyhedron = nullptr;

   return *this;
}

////////////////////////////////////////////////////////////////////////
//
// Return true if tetrahedron is degenerate
// Tetrahedron is concidered as degenerate in case if its minimal
// height is less than degeneracy tolerance
//
G4bool G4Tet::CheckDegeneracy(const G4ThreeVector& p0,
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
// Set data members
//
void G4Tet::Initialize(const G4ThreeVector& p0,
                       const G4ThreeVector& p1,
                       const G4ThreeVector& p2,
                       const G4ThreeVector& p3)
{
  // Set vertices
  fVertex[0] = p0;
  fVertex[1] = p1;
  fVertex[2] = p2;
  fVertex[3] = p3;

  G4ThreeVector norm[4];
  norm[0] = (p2 - p0).cross(p1 - p0);
  norm[1] = (p3 - p0).cross(p2 - p0);
  norm[2] = (p1 - p0).cross(p3 - p0);
  norm[3] = (p2 - p1).cross(p3 - p1);
  G4double volume = norm[0].dot(p3 - p0);
  if (volume > 0.)
  {
    for (G4int i = 0; i < 4; ++i) { norm[i] = -norm[i]; }
  }

  // Set normals to face planes
  for (G4int i = 0; i < 4; ++i) { fNormal[i] = norm[i].unit(); }

  // Set distances to planes
  for (G4int i = 0; i < 3; ++i) { fDist[i] = fNormal[i].dot(p0); }
  fDist[3] = fNormal[3].dot(p1);

  // Set face areas
  for (G4int i = 0; i < 4; ++i) { fArea[i] = 0.5*norm[i].mag(); }

  // Set bounding box
  for (G4int i = 0; i < 3; ++i)
  {
    fBmin[i] = std::min(std::min(std::min(p0[i], p1[i]), p2[i]), p3[i]);
    fBmax[i] = std::max(std::max(std::max(p0[i], p1[i]), p2[i]), p3[i]);
  }

  // Set volume and surface area
  fCubicVolume = std::abs(volume)/6.;
  fSurfaceArea = fArea[0] + fArea[1] + fArea[2] + fArea[3];
}

////////////////////////////////////////////////////////////////////////
//
// Set vertices
//
void G4Tet::SetVertices(const G4ThreeVector& p0,
                        const G4ThreeVector& p1,
                        const G4ThreeVector& p2,
                        const G4ThreeVector& p3, G4bool* degeneracyFlag)
{
  // Check for degeneracy
  G4bool degenerate = CheckDegeneracy(p0, p1, p2, p3);
  if (degeneracyFlag)
  {
    *degeneracyFlag = degenerate;
  }
  else if (degenerate)
  {
    std::ostringstream message;
    message << "Degenerate tetrahedron is not permitted: " << GetName() << " !\n"
            << "  anchor: " << p0 << "\n"
            << "  p1    : " << p1 << "\n"
            << "  p2    : " << p2 << "\n"
            << "  p3    : " << p3 << "\n"
            << "  volume: "
            << std::abs((p1 - p0).cross(p2 - p0).dot(p3 - p0))/6.;
    G4Exception("G4Tet::SetVertices()", "GeomSolids0002",
                FatalException, message);
  }

  // Set data members
  Initialize(p0, p1, p2, p3);

  // Set flag to rebuild polyhedron
  fRebuildPolyhedron = true;
}

////////////////////////////////////////////////////////////////////////
//
// Return four vertices
//
void G4Tet::GetVertices(G4ThreeVector& p0,
                        G4ThreeVector& p1,
                        G4ThreeVector& p2,
                        G4ThreeVector& p3) const
{
  p0 = fVertex[0];
  p1 = fVertex[1];
  p2 = fVertex[2];
  p3 = fVertex[3];
}

////////////////////////////////////////////////////////////////////////
//
// Return std::vector of vertices
//
std::vector<G4ThreeVector> G4Tet::GetVertices() const
{
  std::vector<G4ThreeVector> vertices(4);
  for (G4int i = 0; i < 4; ++i) { vertices[i] = fVertex[i]; }
  return vertices;
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4Tet::ComputeDimensions(G4VPVParameterisation* ,
                              const G4int ,
                              const G4VPhysicalVolume* )
{
}

////////////////////////////////////////////////////////////////////////
//
// Set bounding box
//
void G4Tet::SetBoundingLimits(const G4ThreeVector& pMin,
                              const G4ThreeVector& pMax)
{
  G4int iout[4] = { 0, 0, 0, 0 };
  for (G4int i = 0; i < 4; ++i)
  {
    iout[i] = (fVertex[i].x() < pMin.x() ||
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
            << "    anchor " << fVertex[0] << ((iout[0]) ? " is outside\n" : "\n")
            << "    p1 "     << fVertex[1] << ((iout[1]) ? " is outside\n" : "\n")
            << "    p2 "     << fVertex[2] << ((iout[2]) ? " is outside\n" : "\n")
            << "    p3 "     << fVertex[3] << ((iout[3]) ? " is outside"   : "");
    G4Exception("G4Tet::SetBoundingLimits()", "GeomSolids0002",
                FatalException, message);
  }
  fBmin = pMin;
  fBmax = pMax;
}

////////////////////////////////////////////////////////////////////////
//
// Return bounding box
//
void G4Tet::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  pMin = fBmin;
  pMax = fBmax;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit
//
G4bool G4Tet::CalculateExtent(const EAxis pAxis,
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
  anchor[0].set(vec[0].x(),vec[0].y(),vec[0].z());

  G4ThreeVectorList base(3);
  base[0].set(vec[1].x(),vec[1].y(),vec[1].z());
  base[1].set(vec[2].x(),vec[2].y(),vec[2].z());
  base[2].set(vec[3].x(),vec[3].y(),vec[3].z());

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &anchor;
  polygons[1] = &base;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  return exists = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
}

////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface
//
EInside G4Tet::Inside(const G4ThreeVector& p) const
{
  G4double dd[4];
  for (G4int i = 0; i < 4; ++i) { dd[i] = fNormal[i].dot(p) - fDist[i]; }

  G4double dist = std::max(std::max(std::max(dd[0], dd[1]), dd[2]), dd[3]);
  return (dist <= -halfTolerance) ?
    kInside : ((dist <= halfTolerance) ? kSurface : kOutside);
}

////////////////////////////////////////////////////////////////////////
//
// Return unit normal to surface at p
//
G4ThreeVector G4Tet::SurfaceNormal( const G4ThreeVector& p) const
{
  G4double k[4];
  for (G4int i = 0; i < 4; ++i)
  {
    k[i] = std::abs(fNormal[i].dot(p) - fDist[i]) <= halfTolerance;
  }
  G4double nsurf = k[0] + k[1] + k[2] + k[3];
  G4ThreeVector norm =
    k[0]*fNormal[0] + k[1]*fNormal[1] + k[2]*fNormal[2] + k[3]*fNormal[3];

  if (nsurf == 1.) return norm;
  else if (nsurf > 1.) return norm.unit(); // edge or vertex
  {
#ifdef G4SPECSDEBUG
    std::ostringstream message;
    G4long oldprc = message.precision(16);
    message << "Point p is not on surface (!?) of solid: "
            << GetName() << "\n";
    message << "Position:\n";
    message << "   p.x() = " << p.x()/mm << " mm\n";
    message << "   p.y() = " << p.y()/mm << " mm\n";
    message << "   p.z() = " << p.z()/mm << " mm";
    G4cout.precision(oldprc);
    G4Exception("G4Tet::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, message );
    DumpInfo();
#endif
    return ApproxSurfaceNormal(p);
  }
}

////////////////////////////////////////////////////////////////////////
//
// Find surface nearest to point and return corresponding normal
// This method normally should not be called
//
G4ThreeVector G4Tet::ApproxSurfaceNormal(const G4ThreeVector& p) const
{
  G4double dist = -DBL_MAX;
  G4int iside = 0;
  for (G4int i = 0; i < 4; ++i)
  {
    G4double d = fNormal[i].dot(p) - fDist[i];
    if (d > dist) { dist = d; iside = i; }
  }
  return fNormal[iside];
}

////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface from outside,
// return kInfinity if no intersection
//
G4double G4Tet::DistanceToIn(const G4ThreeVector& p,
                             const G4ThreeVector& v) const
{
  G4double tin = -DBL_MAX, tout = DBL_MAX;
  for (G4int i = 0; i < 4; ++i)
  {
    G4double cosa = fNormal[i].dot(v);
    G4double dist = fNormal[i].dot(p) - fDist[i];
    if (dist >= -halfTolerance)
    {
      if (cosa >= 0.) { return kInfinity; }
      tin = std::max(tin, -dist/cosa);
    }
    else if (cosa > 0.)
    {
      tout = std::min(tout, -dist/cosa);
    }
  }

  return (tout - tin <= halfTolerance) ?
    kInfinity : ((tin < halfTolerance) ? 0. : tin);
}

////////////////////////////////////////////////////////////////////////
//
// Estimate safety distance to surface from outside
//
G4double G4Tet::DistanceToIn(const G4ThreeVector& p) const
{
  G4double dd[4];
  for (G4int i = 0; i < 4; ++i) { dd[i] = fNormal[i].dot(p) - fDist[i]; }

  G4double dist = std::max(std::max(std::max(dd[0], dd[1]), dd[2]), dd[3]);
  return (dist > 0.) ? dist : 0.;
}

////////////////////////////////////////////////////////////////////////
//
// Calcluate distance to surface from inside
//
G4double G4Tet::DistanceToOut(const G4ThreeVector& p,
                              const G4ThreeVector& v,
                              const G4bool calcNorm,
                                    G4bool* validNorm,
                                    G4ThreeVector* n) const
{
  // Calculate distances and cosines
  G4double cosa[4], dist[4];
  G4int ind[4] = {0}, nside = 0;
  for (G4int i = 0; i < 4; ++i)
  {
    G4double tmp = fNormal[i].dot(v);
    cosa[i] = tmp;
    ind[nside] = (tmp > 0) * i;
    nside += (tmp > 0);
    dist[i] = fNormal[i].dot(p) - fDist[i];
  }

  // Find intersection (in most of cases nside == 1)
  G4double tout = DBL_MAX;
  G4int iside = 0;
  for (G4int i = 0; i < nside; ++i)
  {
    G4int k = ind[i];
    // Check: leaving the surface
    if (dist[k] >= -halfTolerance) { tout = 0.; iside = k; break; }
    // Compute distance to intersection
    G4double tmp = -dist[k]/cosa[k];
    if (tmp < tout) { tout = tmp; iside = k; }
  }

  // Set normal, if required, and return distance to out
  if (calcNorm)
  {
    *validNorm = true;
    *n = fNormal[iside];
  }
  return tout;
}

////////////////////////////////////////////////////////////////////////
//
// Calculate safety distance to surface from inside
//
G4double G4Tet::DistanceToOut(const G4ThreeVector& p) const
{
  G4double dd[4];
  for (G4int i = 0; i < 4; ++i) { dd[i] = fDist[i] - fNormal[i].dot(p); }

  G4double dist = std::min(std::min(std::min(dd[0], dd[1]), dd[2]), dd[3]);
  return (dist > 0.) ? dist : 0.;
}

////////////////////////////////////////////////////////////////////////
//
// GetEntityType
//
G4GeometryType G4Tet::GetEntityType() const
{
  return G4String("G4Tet");
}

////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4Tet::Clone() const
{
  return new G4Tet(*this);
}

////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream
//
std::ostream& G4Tet::StreamInfo(std::ostream& os) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: " << GetEntityType() << "\n"
     << " Parameters: \n"
     << "    anchor: " << fVertex[0]/mm << " mm\n"
     << "    p1    : " << fVertex[1]/mm << " mm\n"
     << "    p2    : " << fVertex[2]/mm << " mm\n"
     << "    p3    : " << fVertex[3]/mm << " mm\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);
  return os;
}

////////////////////////////////////////////////////////////////////////
//
// Return random point on the surface
//
G4ThreeVector G4Tet::GetPointOnSurface() const
{
  constexpr G4int iface[4][3] = { {0,1,2}, {0,2,3}, {0,3,1}, {1,2,3} };

  // Select face
  G4double select = fSurfaceArea*G4QuickRand();
  G4int i = 0;
  for ( ; i < 4; ++i) { if ((select -= fArea[i]) <= 0.) break; }

  // Set selected triangle
  G4ThreeVector p0 = fVertex[iface[i][0]];
  G4ThreeVector e1 = fVertex[iface[i][1]] - p0;
  G4ThreeVector e2 = fVertex[iface[i][2]] - p0;

  // Return random point
  G4double r1 = G4QuickRand();
  G4double r2 = G4QuickRand();
  return (r1 + r2 > 1.) ?
    p0 + e1*(1. - r1) + e2*(1. - r2) : p0 + e1*r1 + e2*r2;
}

////////////////////////////////////////////////////////////////////////
//
// Return volume of the tetrahedron
//
G4double G4Tet::GetCubicVolume()
{
  return fCubicVolume;
}

////////////////////////////////////////////////////////////////////////
//
// Return surface area of the tetrahedron
//
G4double G4Tet::GetSurfaceArea()
{
  return fSurfaceArea;
}

////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation
//
void G4Tet::DescribeYourselfTo (G4VGraphicsScene& scene) const
{
  scene.AddSolid (*this);
}

////////////////////////////////////////////////////////////////////////
//
// Return VisExtent
//
G4VisExtent G4Tet::GetExtent() const
{
  return G4VisExtent(fBmin.x(), fBmax.x(),
                     fBmin.y(), fBmax.y(),
                     fBmin.z(), fBmax.z());
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4Tet::CreatePolyhedron() const
{
  // Check orientation of vertices
  G4ThreeVector v1 = fVertex[1] - fVertex[0];
  G4ThreeVector v2 = fVertex[2] - fVertex[0];
  G4ThreeVector v3 = fVertex[3] - fVertex[0];
  G4bool invert = v1.cross(v2).dot(v3) < 0.;
  G4int k2 = (invert) ? 3 : 2;
  G4int k3 = (invert) ? 2 : 3;

  // Set coordinates of vertices
  G4double xyz[4][3];
  for (G4int i = 0; i < 3; ++i)
  {
    xyz[0][i] = fVertex[0][i];
    xyz[1][i] = fVertex[1][i];
    xyz[2][i] = fVertex[k2][i];
    xyz[3][i] = fVertex[k3][i];
  }

  // Create polyhedron
  G4int faces[4][4] = { {1,3,2,0}, {1,4,3,0}, {1,2,4,0}, {2,3,4,0} };
  G4Polyhedron* ph = new G4Polyhedron;
  ph->createPolyhedron(4,4,xyz,faces);

  return ph;
}

////////////////////////////////////////////////////////////////////////
//
// GetPolyhedron
//
G4Polyhedron* G4Tet::GetPolyhedron() const
{
  if (fpPolyhedron == nullptr ||
      fRebuildPolyhedron ||
      fpPolyhedron->GetNumberOfRotationStepsAtTimeOfCreation() !=
      fpPolyhedron->GetNumberOfRotationSteps())
  {
    G4AutoLock l(&polyhedronMutex);
    delete fpPolyhedron;
    fpPolyhedron = CreatePolyhedron();
    fRebuildPolyhedron = false;
    l.unlock();
  }
  return fpPolyhedron;
}

#endif
