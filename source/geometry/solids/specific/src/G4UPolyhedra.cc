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

#include "G4Polyhedra.hh"
#include "G4UPolyhedra.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4GeometryTolerance.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor (GEANT3 style parameters)
//
// GEANT3 PGON radii are specified in the distance to the norm of each face.
//  
G4UPolyhedra::G4UPolyhedra(const G4String& name, 
                                 G4double phiStart,
                                 G4double phiTotal,
                                 G4int numSide,  
                                 G4int numZPlanes,
                           const G4double zPlane[],
                           const G4double rInner[],
                           const G4double rOuter[]  )
  : Base_t(name, phiStart, phiTotal, numSide,
           numZPlanes, zPlane, rInner, rOuter)
{
  fGenericPgon = false;
  SetOriginalParameters();
  wrStart = phiStart;
  while (wrStart < 0)
  {
    wrStart += twopi;
  }
  wrDelta = phiTotal;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrDelta = twopi;
  }
  wrNumSide = numSide;
  G4double convertRad = 1./std::cos(0.5*wrDelta/wrNumSide);
  rzcorners.resize(0);
  for (G4int i=0; i<numZPlanes; ++i)
  {
    G4double z = zPlane[i];
    G4double r = rOuter[i]*convertRad;
    rzcorners.push_back(G4TwoVector(r,z));
  }
  for (G4int i=numZPlanes-1; i>=0; --i)
  {
    G4double z = zPlane[i];
    G4double r = rInner[i]*convertRad;
    rzcorners.push_back(G4TwoVector(r,z));
  }
  std::vector<G4int> iout;
  G4GeomTools::RemoveRedundantVertices(rzcorners,iout,2*kCarTolerance);
}


////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UPolyhedra::G4UPolyhedra(const G4String& name, 
                                 G4double phiStart,
                                 G4double phiTotal,
                                 G4int    numSide,  
                                 G4int    numRZ,
                           const G4double r[],
                           const G4double z[]   )
  : Base_t(name, phiStart, phiTotal, numSide, numRZ/2, r, z)
{
  fGenericPgon = true;
  SetOriginalParameters();
  wrStart = phiStart;
  while (wrStart < 0)
  {
    wrStart += twopi;
  }
  wrDelta = phiTotal;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrDelta = twopi;
  }
  wrNumSide = numSide;
  G4double convertRad = 1./std::cos(0.5*wrDelta/wrNumSide);
  rzcorners.resize(0);
  for (G4int i=0; i<numRZ; ++i)
  {
    rzcorners.push_back(G4TwoVector(r[i]*convertRad,z[i]));
  }
  std::vector<G4int> iout;
  G4GeomTools::RemoveRedundantVertices(rzcorners,iout,2*kCarTolerance);
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UPolyhedra::G4UPolyhedra( __void__& a )
  : Base_t(a)
{
}


////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UPolyhedra::~G4UPolyhedra()
{
}


////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UPolyhedra::G4UPolyhedra( const G4UPolyhedra &source )
  : Base_t( source )
{
  fGenericPgon = source.fGenericPgon;
  fOriginalParameters = source.fOriginalParameters;
  wrStart   = source.wrStart;
  wrDelta   = source.wrDelta;
  wrNumSide = source.wrNumSide;
  rzcorners = source.rzcorners;
}


////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UPolyhedra& G4UPolyhedra::operator=( const G4UPolyhedra &source )
{
  if (this == &source) return *this;

  Base_t::operator=( source );
  fGenericPgon = source.fGenericPgon;
  fOriginalParameters = source.fOriginalParameters;
  wrStart   = source.wrStart;
  wrDelta   = source.wrDelta;
  wrNumSide = source.wrNumSide;
  rzcorners = source.rzcorners;

  return *this;
}


////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers
//
G4int G4UPolyhedra::GetNumSide() const
{
  return wrNumSide;
}
G4double G4UPolyhedra::GetStartPhi() const
{
  return wrStart;
}
G4double G4UPolyhedra::GetEndPhi() const
{
  return (wrStart + wrDelta);
}
G4double G4UPolyhedra::GetSinStartPhi() const
{
  G4double phi = GetStartPhi();
  return std::sin(phi);
}
G4double G4UPolyhedra::GetCosStartPhi() const
{
  G4double phi = GetStartPhi();
  return std::cos(phi);
}
G4double G4UPolyhedra::GetSinEndPhi() const
{
  G4double phi = GetEndPhi();
  return std::sin(phi);
}
G4double G4UPolyhedra::GetCosEndPhi() const
{
  G4double phi = GetEndPhi();
  return std::cos(phi);
}
G4bool G4UPolyhedra::IsOpen() const
{
  return (wrDelta < twopi);
}
G4bool G4UPolyhedra::IsGeneric() const
{
  return fGenericPgon;
}
G4int G4UPolyhedra::GetNumRZCorner() const
{
  return rzcorners.size();
}
G4PolyhedraSideRZ G4UPolyhedra::GetCorner(G4int index) const
{
  G4TwoVector rz = rzcorners.at(index);
  G4PolyhedraSideRZ psiderz = { rz.x(), rz.y() };

  return psiderz;
}
G4PolyhedraHistorical* G4UPolyhedra::GetOriginalParameters() const
{
  return new G4PolyhedraHistorical(fOriginalParameters);
}
void G4UPolyhedra::SetOriginalParameters()
{
  G4int numPlanes = GetZSegmentCount() + 1;
  delete [] fOriginalParameters.Z_values;
  delete [] fOriginalParameters.Rmin;
  delete [] fOriginalParameters.Rmax;
  fOriginalParameters.Z_values = new G4double[numPlanes];
  fOriginalParameters.Rmin = new G4double[numPlanes];
  fOriginalParameters.Rmax = new G4double[numPlanes];

  for (G4int j=0; j<numPlanes; ++j)
  {
    fOriginalParameters.Z_values[j] = GetZPlanes()[j];
    fOriginalParameters.Rmax[j]     = GetRMax()[j];
    fOriginalParameters.Rmin[j]     = GetRMin()[j];
  }

  fOriginalParameters.Start_angle   = GetPhiStart();
  fOriginalParameters.Opening_angle = GetPhiDelta();
  fOriginalParameters.Num_z_planes  = numPlanes;
  fOriginalParameters.numSide       = GetSideCount();
}
void G4UPolyhedra::SetOriginalParameters(G4PolyhedraHistorical* pars)
{
  fOriginalParameters = *pars;
  fRebuildPolyhedron = true;
  Reset();
}

G4bool G4UPolyhedra::Reset()
{
  if (fGenericPgon)
  {
    std::ostringstream message;
    message << "Solid " << GetName() << " built using generic construct."
            << G4endl << "Not applicable to the generic construct !";
    G4Exception("G4UPolyhedra::Reset()", "GeomSolids1001",
                JustWarning, message, "Parameters NOT resetted.");
    return true;  // error code set
  }

  //
  // Rebuild polyhedra based on original parameters
  //
  wrStart = fOriginalParameters.Start_angle;
  while (wrStart < 0)
  {
    wrStart += twopi;
  }
  wrDelta = fOriginalParameters.Opening_angle;
  if (wrDelta <= 0 || wrDelta >= twopi*(1-DBL_EPSILON))
  {
    wrDelta = twopi;
  }
  wrNumSide = fOriginalParameters.numSide;
  G4double convertRad = 1./std::cos(0.5*wrDelta/wrNumSide);
  rzcorners.resize(0);
  for (G4int i=0; i<fOriginalParameters.Num_z_planes; ++i)
  {
    G4double z = fOriginalParameters.Z_values[i];
    G4double r = fOriginalParameters.Rmax[i]*convertRad;
    rzcorners.push_back(G4TwoVector(r,z));
  }
  for (G4int i=fOriginalParameters.Num_z_planes-1; i>=0; --i)
  {
    G4double z = fOriginalParameters.Z_values[i];
    G4double r = fOriginalParameters.Rmin[i]*convertRad;
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
void G4UPolyhedra::ComputeDimensions(G4VPVParameterisation* p,
                                     const G4int n,
                                     const G4VPhysicalVolume* pRep)
{
  p->ComputeDimensions(*(G4Polyhedra*)this,n,pRep);
}


//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UPolyhedra::Clone() const
{
  return new G4UPolyhedra(*this);
}


//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UPolyhedra::BoundingLimits(G4ThreeVector& pMin,
                                  G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;
  static G4bool checkPhi  = true;

  G4double rmin = kInfinity, rmax = -kInfinity;
  G4double zmin = kInfinity, zmax = -kInfinity;
  for (G4int i=0; i<GetNumRZCorner(); ++i)
  {
    G4PolyhedraSideRZ corner = GetCorner(i);
    if (corner.r < rmin) rmin = corner.r;
    if (corner.r > rmax) rmax = corner.r;
    if (corner.z < zmin) zmin = corner.z;
    if (corner.z > zmax) zmax = corner.z;
  }

  G4double sphi    = GetStartPhi();
  G4double ephi    = GetEndPhi();
  G4double dphi    = IsOpen() ? ephi-sphi : twopi;
  G4int    ksteps  = GetNumSide();
  G4double astep   = dphi/ksteps;
  G4double sinStep = std::sin(astep);
  G4double cosStep = std::cos(astep);

  G4double sinCur = GetSinStartPhi();
  G4double cosCur = GetCosStartPhi();
  if (!IsOpen()) rmin = 0;
  G4double xmin = rmin*cosCur, xmax = xmin;
  G4double ymin = rmin*sinCur, ymax = ymin;
  for (G4int k=0; k<ksteps+1; ++k)
  {
    G4double x = rmax*cosCur;
    if (x < xmin) xmin = x;
    if (x > xmax) xmax = x;
    G4double y = rmax*sinCur;
    if (y < ymin) ymin = y;
    if (y > ymax) ymax = y;
    if (rmin > 0)
    {
      G4double xx = rmin*cosCur;
      if (xx < xmin) xmin = xx;
      if (xx > xmax) xmax = xx;
      G4double yy = rmin*sinCur;
      if (yy < ymin) ymin = yy;
      if (yy > ymax) ymax = yy;
    }
    G4double sinTmp = sinCur;
    sinCur = sinCur*cosStep + cosCur*sinStep;
    cosCur = cosCur*cosStep - sinTmp*sinStep;
  }
  pMin.set(xmin,ymin,zmin);
  pMax.set(xmax,ymax,zmax);

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UPolyhedra::BoundingLimits()", "GeomMgt0001",
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
      G4Exception("G4UPolyhedra::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }

  // Check consistency of angles
  //
  if (checkPhi)
  {
    if (GetStartPhi() != GetPhiStart() ||
	GetEndPhi()   != GetPhiEnd()   ||
	GetNumSide()  != GetSideCount() ||
        IsOpen()      != (Base_t::GetPhiDelta() < twopi))
    {
      std::ostringstream message;
      message << "Inconsistency in Phi angles or # of sides for solid: "
              << GetName() << " !"
              << "\nPhi start  : wrapper = " << GetStartPhi()
              << " solid = " <<     GetPhiStart()
              << "\nPhi end    : wrapper = " << GetEndPhi()
              << " solid = " <<     GetPhiEnd()
              << "\nPhi # sides: wrapper = " << GetNumSide()
              << " solid = " <<     GetSideCount()
              << "\nPhi is open: wrapper = " << (IsOpen() ? "true" : "false")
              << " solid = "
              << ((Base_t::GetPhiDelta() < twopi) ? "true" : "false");
      G4Exception("G4UPolyhedra::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkPhi = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UPolyhedra::CalculateExtent(const EAxis pAxis,
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
    G4PolyhedraSideRZ corner = GetCorner(i);
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
    G4Exception("G4UPolyhedra::CalculateExtent()",
                "GeomMgt1002",JustWarning,message);
    return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }

  // set trigonometric values
  G4double sphi     = GetStartPhi();
  G4double ephi     = GetEndPhi();
  G4double dphi     = IsOpen() ? ephi-sphi : twopi;
  G4int    ksteps   = GetNumSide();
  G4double astep    = dphi/ksteps;
  G4double sinStep  = std::sin(astep);
  G4double cosStep  = std::cos(astep);
  G4double sinStart = GetSinStartPhi();
  G4double cosStart = GetCosStartPhi();

  // allocate vector lists
  std::vector<const G4ThreeVectorList *> polygons;
  polygons.resize(ksteps+1);
  for (G4int k=0; k<ksteps+1; ++k) {
    polygons[k] = new G4ThreeVectorList(3);
  }

  // main loop along triangles
  pMin =  kInfinity;
  pMax = -kInfinity;
  G4int ntria = triangles.size()/3;
  for (G4int i=0; i<ntria; ++i)
  {
    G4double sinCur = sinStart;
    G4double cosCur = cosStart;
    G4int i3 = i*3;
    for (G4int k=0; k<ksteps+1; ++k) // rotate triangle
    {
      G4ThreeVectorList* ptr = const_cast<G4ThreeVectorList*>(polygons[k]);
      G4ThreeVectorList::iterator iter = ptr->begin();
      iter->set(triangles[i3+0].x()*cosCur,
                triangles[i3+0].x()*sinCur,
                triangles[i3+0].y());
      iter++;
      iter->set(triangles[i3+1].x()*cosCur,
                triangles[i3+1].x()*sinCur,
                triangles[i3+1].y());
      iter++;
      iter->set(triangles[i3+2].x()*cosCur,
                triangles[i3+2].x()*sinCur,
                triangles[i3+2].y());

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }

    // set sub-envelope and adjust extent
    G4double emin,emax;
    G4BoundingEnvelope benv(polygons);
    if (!benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,emin,emax)) continue;
    if (emin < pMin) pMin = emin;
    if (emax > pMax) pMax = emax;
    if (eminlim > pMin && emaxlim < pMax) break; // max possible extent
  }
  // free memory
  for (G4int k=0; k<ksteps+1; ++k) { delete polygons[k]; polygons[k]=0;}
  return (pMin < pMax);
}


////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron
//
G4Polyhedron* G4UPolyhedra::CreatePolyhedron() const
{
  if (!IsGeneric())
  {
    return new G4PolyhedronPgon( fOriginalParameters.Start_angle,
                                 fOriginalParameters.Opening_angle,
                                 fOriginalParameters.numSide,
                                 fOriginalParameters.Num_z_planes,
                                 fOriginalParameters.Z_values,
                                 fOriginalParameters.Rmin,
                                 fOriginalParameters.Rmax);
  }
  else
  {
    // The following code prepares for:
    // HepPolyhedron::createPolyhedron(int Nnodes, int Nfaces,
    //                                 const double xyz[][3],
    //                                 const int faces_vec[][4])
    // Here is an extract from the header file HepPolyhedron.h:
    /**
     * Creates user defined polyhedron.
     * This function allows to the user to define arbitrary polyhedron.
     * The faces of the polyhedron should be either triangles or planar
     * quadrilateral. Nodes of a face are defined by indexes pointing to
     * the elements in the xyz array. Numeration of the elements in the
     * array starts from 1 (like in fortran). The indexes can be positive
     * or negative. Negative sign means that the corresponding edge is
     * invisible. The normal of the face should be directed to exterior
     * of the polyhedron. 
     * 
     * @param  Nnodes number of nodes
     * @param  Nfaces number of faces
     * @param  xyz    nodes
     * @param  faces_vec  faces (quadrilaterals or triangles)
     * @return status of the operation - is non-zero in case of problem
     */
    G4int nNodes;
    G4int nFaces;
    typedef G4double double3[3];
    double3* xyz;
    typedef G4int int4[4];
    int4* faces_vec;
    if (IsOpen())
    {
      // Triangulate open ends.  Simple ear-chopping algorithm...
      // I'm not sure how robust this algorithm is (J.Allison).
      //
      std::vector<G4bool> chopped(GetNumRZCorner(), false);
      std::vector<G4int*> triQuads;
      G4int remaining = GetNumRZCorner();
      G4int iStarter = 0;
      while (remaining >= 3)    // Loop checking, 13.08.2015, G.Cosmo
      {
        // Find unchopped corners...
        //
        G4int A = -1, B = -1, C = -1;
        G4int iStepper = iStarter;
        do    // Loop checking, 13.08.2015, G.Cosmo
        {
          if (A < 0)      { A = iStepper; }
          else if (B < 0) { B = iStepper; }
          else if (C < 0) { C = iStepper; }
          do    // Loop checking, 13.08.2015, G.Cosmo
          {
            if (++iStepper >= GetNumRZCorner()) iStepper = 0;
          }
          while (chopped[iStepper]);
        }
        while (C < 0 && iStepper != iStarter);

        // Check triangle at B is pointing outward (an "ear").
        // Sign of z cross product determines...

        G4double BAr = GetCorner(A).r - GetCorner(B).r;
        G4double BAz = GetCorner(A).z - GetCorner(B).z;
        G4double BCr = GetCorner(C).r - GetCorner(B).r;
        G4double BCz = GetCorner(C).z - GetCorner(B).z;
        if (BAr * BCz - BAz * BCr < kCarTolerance)
        {
          G4int* tq = new G4int[3];
          tq[0] = A + 1;
          tq[1] = B + 1;
          tq[2] = C + 1;
          triQuads.push_back(tq);
          chopped[B] = true;
          --remaining;
        }
        else
        {
          do    // Loop checking, 13.08.2015, G.Cosmo
          {
            if (++iStarter >= GetNumRZCorner()) { iStarter = 0; }
          }
          while (chopped[iStarter]);
        }
      }

      // Transfer to faces...
      G4int numSide=GetNumSide();
      nNodes = (numSide + 1) * GetNumRZCorner();
      nFaces = numSide * GetNumRZCorner() + 2 * triQuads.size();
      faces_vec = new int4[nFaces];
      G4int iface = 0;
      G4int addition = GetNumRZCorner() * numSide;
      G4int d = GetNumRZCorner() - 1;
      for (G4int iEnd = 0; iEnd < 2; ++iEnd)
      {
        for (size_t i = 0; i < triQuads.size(); ++i)
        {
          // Negative for soft/auxiliary/normally invisible edges...
          //
          G4int a, b, c;
          if (iEnd == 0)
          {
            a = triQuads[i][0];
            b = triQuads[i][1];
            c = triQuads[i][2];
          }
          else
          {
            a = triQuads[i][0] + addition;
            b = triQuads[i][2] + addition;
            c = triQuads[i][1] + addition;
          }
          G4int ab = std::abs(b - a);
          G4int bc = std::abs(c - b);
          G4int ca = std::abs(a - c);
          faces_vec[iface][0] = (ab == 1 || ab == d)? a: -a;
          faces_vec[iface][1] = (bc == 1 || bc == d)? b: -b;
          faces_vec[iface][2] = (ca == 1 || ca == d)? c: -c;
          faces_vec[iface][3] = 0;
          ++iface;
        }
      }

      // Continue with sides...

      xyz = new double3[nNodes];
      const G4double dPhi = (GetEndPhi() - GetStartPhi()) / numSide;
      G4double phi = GetStartPhi();
      G4int ixyz = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
        {
          xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
          xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
          xyz[ixyz][2] = GetCorner(iCorner).z;
          if (iCorner < GetNumRZCorner() - 1)
          {
            faces_vec[iface][0] = ixyz + 1;
            faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
            faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
            faces_vec[iface][3] = ixyz + 2;
          }
          else
          {
            faces_vec[iface][0] = ixyz + 1;
            faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
            faces_vec[iface][2] = ixyz + 2;
            faces_vec[iface][3] = ixyz - GetNumRZCorner() + 2;
          }
          ++iface;
          ++ixyz;
        }
        phi += dPhi;
      }

      // Last GetCorner...

      for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
      {
        xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
        xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
        xyz[ixyz][2] = GetCorner(iCorner).z;
        ++ixyz;
      }
    }
    else  // !phiIsOpen - i.e., a complete 360 degrees.
    {
      nNodes = GetNumSide() * GetNumRZCorner();
      nFaces = GetNumSide() * GetNumRZCorner();;
      xyz = new double3[nNodes];
      faces_vec = new int4[nFaces];
      // const G4double dPhi = (endPhi - startPhi) / numSide;
      const G4double dPhi = twopi / GetNumSide();
      // !phiIsOpen endPhi-startPhi = 360 degrees.
      G4double phi = GetStartPhi();
      G4int ixyz = 0, iface = 0;
      for (G4int iSide = 0; iSide < GetNumSide(); ++iSide)
      {
        for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
        {
          xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
          xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
          xyz[ixyz][2] = GetCorner(iCorner).z;
          if (iSide < GetNumSide() - 1)
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
              faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = ixyz - GetNumRZCorner() + 2;
            }
          }
          else   // Last side joins ends...
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + GetNumRZCorner() - nFaces + 1;
              faces_vec[iface][2] = ixyz + GetNumRZCorner() - nFaces + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz - nFaces + GetNumRZCorner() + 1;
              faces_vec[iface][2] = ixyz - nFaces + 2;
              faces_vec[iface][3] = ixyz - GetNumRZCorner() + 2;
            }
          }
          ++ixyz;
          ++iface;
        }
        phi += dPhi;
      }
    }
    G4Polyhedron* polyhedron = new G4Polyhedron;
    G4int prob = polyhedron->createPolyhedron(nNodes, nFaces, xyz, faces_vec);
    delete [] faces_vec;
    delete [] xyz;
    if (prob)
    {
      std::ostringstream message;
      message << "Problem creating G4Polyhedron for: " << GetName();
      G4Exception("G4Polyhedra::CreatePolyhedron()", "GeomSolids1002",
                  JustWarning, message);
      delete polyhedron;
      return 0;
    }
    else
    {
      return polyhedron;
    }
  }
}

#endif  // G4GEOM_USE_USOLIDS
