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
// Implementation of G4UPolycone wrapper class
// --------------------------------------------------------------------

#include "G4Polycone.hh"
#include "G4UPolycone.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor (GEANT3 style parameters)
//  
G4UPolycone::G4UPolycone( const G4String& name, 
                              G4double phiStart,
                              G4double phiTotal,
                              G4int numZPlanes,
                        const G4double zPlane[],
                        const G4double rInner[],
                        const G4double rOuter[]  )
  : Base_t(name, phiStart, phiTotal, numZPlanes, zPlane, rInner, rOuter)
{
  fGenericPcon = false;
  SetOriginalParameters();
  wrStart = phiStart;
  while (wrStart < 0)
  {
    wrStart += twopi;
  }
  wrDelta = phiTotal;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrStart = 0;
    wrDelta = twopi;
  }
  rzcorners.resize(0);
  for (G4int i=0; i<numZPlanes; ++i)
  {
    G4double z = zPlane[i];
    G4double r = rOuter[i];
    rzcorners.push_back(G4TwoVector(r,z));
  }
  for (G4int i=numZPlanes-1; i>=0; --i)
  {
    G4double z = zPlane[i];
    G4double r = rInner[i];
    rzcorners.push_back(G4TwoVector(r,z));
  }
  std::vector<G4int> iout;
  G4GeomTools::RemoveRedundantVertices(rzcorners,iout,2*kCarTolerance);
}


////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UPolycone::G4UPolycone(const G4String& name, 
                               G4double phiStart,
                               G4double phiTotal,
                               G4int    numRZ,
                         const G4double r[],
                         const G4double z[]   )
  : Base_t(name, phiStart, phiTotal, numRZ, r, z)
{ 
  fGenericPcon = true;
  SetOriginalParameters();
  wrStart = phiStart; while (wrStart < 0) wrStart += twopi;
  wrDelta = phiTotal;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrStart = 0;
    wrDelta = twopi;
  }
  rzcorners.resize(0);
  for (G4int i=0; i<numRZ; ++i)
  {
    rzcorners.push_back(G4TwoVector(r[i],z[i]));
  }
  std::vector<G4int> iout;
  G4GeomTools::RemoveRedundantVertices(rzcorners,iout,2*kCarTolerance);
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UPolycone::G4UPolycone( __void__& a )
  : Base_t(a)
{
}


////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UPolycone::~G4UPolycone()
{
}


////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UPolycone::G4UPolycone( const G4UPolycone &source )
  : Base_t( source )
{
  fGenericPcon = source.fGenericPcon;
  fOriginalParameters = source.fOriginalParameters;
  wrStart   = source.wrStart;
  wrDelta   = source.wrDelta;
  rzcorners = source.rzcorners;
}


////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UPolycone &G4UPolycone::operator=( const G4UPolycone &source )
{
  if (this == &source) return *this;
  
  Base_t::operator=( source );
  fGenericPcon = source.fGenericPcon;
  fOriginalParameters = source.fOriginalParameters;
  wrStart   = source.wrStart;
  wrDelta   = source.wrDelta;
  rzcorners = source.rzcorners;

  return *this;
}


////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers
//
G4double G4UPolycone::GetStartPhi() const
{
  return wrStart;
}
G4double G4UPolycone::GetDeltaPhi() const
{
  return wrDelta;
}
G4double G4UPolycone::GetEndPhi() const
{
  return (wrStart + wrDelta);
}
G4double G4UPolycone::GetSinStartPhi() const
{
  if (!IsOpen()) return 0;
  G4double phi = GetStartPhi();
  return std::sin(phi);
}
G4double G4UPolycone::GetCosStartPhi() const
{
  if (!IsOpen()) return 1;
  G4double phi = GetStartPhi();
  return std::cos(phi);
}
G4double G4UPolycone::GetSinEndPhi() const
{
  if (!IsOpen()) return 0;
  G4double phi = GetEndPhi();
  return std::sin(phi);
}
G4double G4UPolycone::GetCosEndPhi() const
{
  if (!IsOpen()) return 1;
  G4double phi = GetEndPhi();
  return std::cos(phi);
}
G4bool G4UPolycone::IsOpen() const
{
  return (wrDelta < twopi);
}
G4int G4UPolycone::GetNumRZCorner() const
{
  return rzcorners.size();
}
G4PolyconeSideRZ G4UPolycone::GetCorner(G4int index) const
{
  G4TwoVector rz = rzcorners.at(index);
  G4PolyconeSideRZ psiderz = { rz.x(), rz.y() };

  return psiderz;
}
G4PolyconeHistorical* G4UPolycone::GetOriginalParameters() const
{
  return new G4PolyconeHistorical(fOriginalParameters);
}
void G4UPolycone::SetOriginalParameters()
{
  G4int numPlanes = GetNSections()+1;
  delete [] fOriginalParameters.Z_values;
  delete [] fOriginalParameters.Rmin;
  delete [] fOriginalParameters.Rmax;
  fOriginalParameters.Z_values = new G4double[numPlanes];
  fOriginalParameters.Rmin = new G4double[numPlanes];
  fOriginalParameters.Rmax = new G4double[numPlanes];

  for (G4int j=0; j<numPlanes; ++j)
  {
    fOriginalParameters.Z_values[j] = GetZAtPlane(j);
    fOriginalParameters.Rmax[j]     = GetRmaxAtPlane(j);
    fOriginalParameters.Rmin[j]     = GetRminAtPlane(j);
  }

  fOriginalParameters.Start_angle   = Base_t::GetStartPhi();
  fOriginalParameters.Opening_angle = Base_t::GetDeltaPhi();
  fOriginalParameters.Num_z_planes  = numPlanes;
}
void G4UPolycone::SetOriginalParameters(G4PolyconeHistorical* pars)
{
  fOriginalParameters = *pars;
  fRebuildPolyhedron = true;
  Reset();
}

G4bool G4UPolycone::Reset()
{
  if (fGenericPcon)
  {
    std::ostringstream message;
    message << "Solid " << GetName() << " built using generic construct."
            << G4endl << "Not applicable to the generic construct !";
    G4Exception("G4UPolycone::Reset()", "GeomSolids1001",
                JustWarning, message, "Parameters NOT resetted.");
    return true;  // error code set
  }

  //
  // Rebuild polycone based on original parameters
  //
  wrStart = fOriginalParameters.Start_angle;
  while (wrStart < 0)
  {
    wrStart += twopi;
  }
  wrDelta = fOriginalParameters.Opening_angle;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrStart = 0;
    wrDelta = twopi;
  }
  rzcorners.resize(0);
  for (G4int i=0; i<fOriginalParameters.Num_z_planes; ++i)
    {
      G4double z = fOriginalParameters.Z_values[i];
      G4double r = fOriginalParameters.Rmax[i];
      rzcorners.push_back(G4TwoVector(r,z));
    }
  for (G4int i=fOriginalParameters.Num_z_planes-1; i>=0; --i)
    {
      G4double z = fOriginalParameters.Z_values[i];
      G4double r = fOriginalParameters.Rmin[i];
      rzcorners.push_back(G4TwoVector(r,z));
    }
  std::vector<G4int> iout;
  G4GeomTools::RemoveRedundantVertices(rzcorners,iout,2*kCarTolerance);

  return false;  // error code unset
}

////////////////////////////////////////////////////////////////////////
//
// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
//
void G4UPolycone::ComputeDimensions(G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Polycone*)this,n,pRep);
}


//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UPolycone::Clone() const
{
  return new G4UPolycone(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UPolycone::BoundingLimits(G4ThreeVector& pMin,
                                 G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;
  static G4bool checkPhi  = true;

  G4double rmin = kInfinity, rmax = -kInfinity;
  G4double zmin = kInfinity, zmax = -kInfinity;

  for (G4int i=0; i<GetNumRZCorner(); ++i)
  {
    G4PolyconeSideRZ corner = GetCorner(i);
    if (corner.r < rmin) rmin = corner.r;
    if (corner.r > rmax) rmax = corner.r;
    if (corner.z < zmin) zmin = corner.z;
    if (corner.z > zmax) zmax = corner.z;
  }

  if (IsOpen())
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rmin,rmax,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            vmin,vmax);
    pMin.set(vmin.x(),vmin.y(),zmin);
    pMax.set(vmax.x(),vmax.y(),zmax);
  }
  else
  {
    pMin.set(-rmax,-rmax, zmin);
    pMax.set( rmax, rmax, zmax);
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
    G4Exception("G4UPolycone::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    U3Vector vmin, vmax;
    Extent(vmin,vmax);
    if (std::abs(pMin.x()-vmin.x()) > kCarTolerance ||
        std::abs(pMin.y()-vmin.y()) > kCarTolerance ||
        std::abs(pMin.z()-vmin.z()) > kCarTolerance ||
        std::abs(pMax.x()-vmax.x()) > kCarTolerance ||
        std::abs(pMax.y()-vmax.y()) > kCarTolerance ||
        std::abs(pMax.z()-vmax.z()) > kCarTolerance)
    {
      std::ostringstream message;
      message << "Inconsistency in bounding boxes for solid: "
              << GetName() << " !"
              << "\nBBox min: wrapper = " << pMin << " solid = " << vmin
              << "\nBBox max: wrapper = " << pMax << " solid = " << vmax;
      G4Exception("G4UPolycone::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }

  // Check consistency of angles
  //
  if (checkPhi)
  {
    if (GetStartPhi() != Base_t::GetStartPhi() ||
        GetEndPhi()   != Base_t::GetEndPhi()   ||
        IsOpen()      != (Base_t::GetDeltaPhi() < twopi))
    {
      std::ostringstream message;
      message << "Inconsistency in Phi angles or # of sides for solid: "
              << GetName() << " !"
              << "\nPhi start  : wrapper = " << GetStartPhi()
              << " solid = " <<     Base_t::GetStartPhi()
              << "\nPhi end    : wrapper = " << GetEndPhi()
              << " solid = " <<     Base_t::GetEndPhi()
              << "\nPhi is open: wrapper = " << (IsOpen() ? "true" : "false")
              << " solid = "
              << ((Base_t::GetDeltaPhi() < twopi) ? "true" : "false");
      G4Exception("G4UPolycone::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkPhi = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4UPolycone::CalculateExtent(const EAxis pAxis,
                                    const G4VoxelLimits& pVoxelLimit,
                                    const G4AffineTransform& pTransform,
                                          G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Check bounding box (bbox)
  //
  BoundingLimits(bmin,bmax);
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // To find the extent, RZ contour of the polycone is subdivided
  // in triangles. The extent is calculated as cumulative extent of
  // all sub-polycones formed by rotation of triangles around Z
  //
  G4TwoVectorList contourRZ;
  G4TwoVectorList triangles;
  std::vector<G4int> iout;
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  // get RZ contour, ensure anticlockwise order of corners
  for (G4int i=0; i<GetNumRZCorner(); ++i)
  {
    G4PolyconeSideRZ corner = GetCorner(i);
    contourRZ.push_back(G4TwoVector(corner.r,corner.z));
  }
  G4GeomTools::RemoveRedundantVertices(contourRZ,iout,2*kCarTolerance);
  G4double area = G4GeomTools::PolygonArea(contourRZ);
  if (area < 0.) std::reverse(contourRZ.begin(),contourRZ.end());

  // triangulate RZ countour
  if (!G4GeomTools::TriangulatePolygon(contourRZ,triangles))
  {
    std::ostringstream message;
    message << "Triangulation of RZ contour has failed for solid: "
            << GetName() << " !"
            << "\nExtent has been calculated using boundary box";
    G4Exception("G4UPolycone::CalculateExtent()",
                "GeomMgt1002", JustWarning, message);
    return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }

  // set trigonometric values
  const G4int NSTEPS = 24;            // number of steps for whole circle
  G4double astep  = twopi/NSTEPS;     // max angle for one step

  G4double sphi   = GetStartPhi();
  G4double ephi   = GetEndPhi();
  G4double dphi   = IsOpen() ? ephi-sphi : twopi;
  G4int    ksteps = (dphi <= astep) ? 1 : (G4int)((dphi-deg)/astep) + 1;
  G4double ang    = dphi/ksteps;

  G4double sinHalf = std::sin(0.5*ang);
  G4double cosHalf = std::cos(0.5*ang);
  G4double sinStep = 2.*sinHalf*cosHalf;
  G4double cosStep = 1. - 2.*sinHalf*sinHalf;

  G4double sinStart = GetSinStartPhi();
  G4double cosStart = GetCosStartPhi();
  G4double sinEnd   = GetSinEndPhi();
  G4double cosEnd   = GetCosEndPhi();

  // define vectors and arrays
  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(ksteps+2);
  G4ThreeVectorList pols[NSTEPS+2];
  for (G4int k=0; k<ksteps+2; ++k) pols[k].resize(6);
  for (G4int k=0; k<ksteps+2; ++k) polygons[k] = &pols[k];
  G4double r0[6],z0[6]; // contour with original edges of triangle
  G4double r1[6];       // shifted radii of external edges of triangle

  // main loop along triangles
  pMin = kInfinity;
  pMax =-kInfinity;
  G4int ntria = triangles.size()/3;
  for (G4int i=0; i<ntria; ++i)
  {
    G4int i3 = i*3;
    for (G4int k=0; k<3; ++k)
    {
      G4int e0 = i3+k, e1 = (k<2) ? e0+1 : i3;
      G4int k2 = k*2;
      // set contour with original edges of triangle
      r0[k2+0] = triangles[e0].x(); z0[k2+0] = triangles[e0].y();
      r0[k2+1] = triangles[e1].x(); z0[k2+1] = triangles[e1].y();
      // set shifted radii
      r1[k2+0] = r0[k2+0];
      r1[k2+1] = r0[k2+1];
      if (z0[k2+1] - z0[k2+0] <= 0) continue;
      r1[k2+0] /= cosHalf;
      r1[k2+1] /= cosHalf;
    }

    // rotate countour, set sequence of 6-sided polygons
    G4double sinCur = sinStart*cosHalf + cosStart*sinHalf;
    G4double cosCur = cosStart*cosHalf - sinStart*sinHalf;
    for (G4int j=0; j<6; ++j) pols[0][j].set(r0[j]*cosStart,r0[j]*sinStart,z0[j]);
    for (G4int k=1; k<ksteps+1; ++k)
    {
      for (G4int j=0; j<6; ++j) pols[k][j].set(r1[j]*cosCur,r1[j]*sinCur,z0[j]);
      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    for (G4int j=0; j<6; ++j) pols[ksteps+1][j].set(r0[j]*cosEnd,r0[j]*sinEnd,z0[j]);

    // set sub-envelope and adjust extent
    G4double emin,emax;
    G4BoundingEnvelope benv(polygons);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) return true; // max possible extent
  }
  return (pMin < pMax);
}

////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UPolycone::CreatePolyhedron() const
{
  G4PolyhedronPcon*
  polyhedron = new G4PolyhedronPcon( fOriginalParameters.Start_angle,
                                     fOriginalParameters.Opening_angle,
                                     fOriginalParameters.Num_z_planes,
                                     fOriginalParameters.Z_values,
                                     fOriginalParameters.Rmin,
                                     fOriginalParameters.Rmax );
  return polyhedron;
}

#endif  // G4GEOM_USE_USOLIDS
