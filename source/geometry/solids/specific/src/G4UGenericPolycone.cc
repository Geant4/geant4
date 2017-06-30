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
// Implementation of G4UGenericPolycone wrapper class
// --------------------------------------------------------------------

#include "G4GenericPolycone.hh"
#include "G4UGenericPolycone.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

#include "G4Polyhedron.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UGenericPolycone::G4UGenericPolycone(const G4String& name, 
                                             G4double phiStart,
                                             G4double phiTotal,
                                             G4int    numRZ,
                                       const G4double r[],
                                       const G4double z[]   )
  : G4USolid(name, new UGenericPolycone(name, phiStart, phiTotal, numRZ, r, z))
{ 
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UGenericPolycone::G4UGenericPolycone(__void__& a)
  : G4USolid(a)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UGenericPolycone::~G4UGenericPolycone()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UGenericPolycone::G4UGenericPolycone(const G4UGenericPolycone &source)
  : G4USolid(source)
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UGenericPolycone&
G4UGenericPolycone::operator=(const G4UGenericPolycone &source)
{
  if (this == &source) return *this;
  
  G4USolid::operator=( source );
  
  return *this;
}

G4double G4UGenericPolycone::GetStartPhi() const
{
  return GetShape()->GetStartPhi();
}
G4double G4UGenericPolycone::GetEndPhi() const
{
  return GetShape()->GetEndPhi();
}
G4double G4UGenericPolycone::GetSinStartPhi() const
{
  if (!GetShape()->IsOpen()) return 0;
  G4double phi = GetShape()->GetStartPhi();
  return std::sin(phi);
}
G4double G4UGenericPolycone::GetCosStartPhi() const
{
  if (!GetShape()->IsOpen()) return 1;
  G4double phi = GetShape()->GetStartPhi();
  return std::cos(phi);
}
G4double G4UGenericPolycone::GetSinEndPhi() const
{
  if (!GetShape()->IsOpen()) return 0;
  G4double phi = GetShape()->GetEndPhi();
  return std::sin(phi);
}
G4double G4UGenericPolycone::GetCosEndPhi() const
{
  if (!GetShape()->IsOpen()) return 1;
  G4double phi = GetShape()->GetEndPhi();
  return std::cos(phi);
}
G4bool G4UGenericPolycone::IsOpen() const
{
  return GetShape()->IsOpen();
}
G4int G4UGenericPolycone::GetNumRZCorner() const
{
  return GetShape()->GetNumRZCorner();
}
G4PolyconeSideRZ G4UGenericPolycone::GetCorner(G4int index) const
{
  UPolyconeSideRZ pside = GetShape()->GetCorner(index);
  G4PolyconeSideRZ psiderz = { pside.r, pside.z };

  return psiderz;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void
G4UGenericPolycone::BoundingLimits(G4ThreeVector& pMin,
                                   G4ThreeVector& pMax) const
{
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
    G4Exception("G4UGenericPolycone::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UGenericPolycone::CalculateExtent(const EAxis pAxis,
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
  G4double eminlim = pVoxelLimit.GetMinExtent(pAxis);
  G4double emaxlim = pVoxelLimit.GetMaxExtent(pAxis);

  // get RZ contour, ensure anticlockwise order of corners
  for (G4int i=0; i<GetNumRZCorner(); ++i)
  {
    G4PolyconeSideRZ corner = GetCorner(i);
    contourRZ.push_back(G4TwoVector(corner.r,corner.z));
  }
  G4double area = G4GeomTools::PolygonArea(contourRZ);
  if (area < 0.) std::reverse(contourRZ.begin(),contourRZ.end());

  // triangulate RZ countour
  if (!G4GeomTools::TriangulatePolygon(contourRZ,triangles))
  {
    std::ostringstream message;
    message << "Triangulation of RZ contour has failed for solid: "
            << GetName() << " !"
            << "\nExtent has been calculated using boundary box";
    G4Exception("G4UGenericPolycone::CalculateExtent()",
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
    for (G4int j=0; j<6; ++j)
    {
      pols[0][j].set(r0[j]*cosStart,r0[j]*sinStart,z0[j]);
    }
    for (G4int k=1; k<ksteps+1; ++k)
    {
      for (G4int j=0; j<6; ++j)
      {
        pols[k][j].set(r1[j]*cosCur,r1[j]*sinCur,z0[j]);
      }
      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    for (G4int j=0; j<6; ++j)
    {
      pols[ksteps+1][j].set(r0[j]*cosEnd,r0[j]*sinEnd,z0[j]);
    }

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

G4Polyhedron* G4UGenericPolycone::CreatePolyhedron() const
{


 // The following code prepares for:
    // HepPolyhedron::createPolyhedron(int Nnodes, int Nfaces,
    //                                  const double xyz[][3],
    //                                  const int faces_vec[][4])
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
    const G4int numSide =
          G4int(G4Polyhedron::GetNumberOfRotationSteps()
                * (GetEndPhi() - GetStartPhi()) / twopi) + 1;
    G4int nNodes;
    G4int nFaces;
    typedef G4double double3[3];
    double3* xyz;
    typedef G4int int4[4];
    int4* faces_vec;
    if (IsOpen())
    {
      // Triangulate open ends. Simple ear-chopping algorithm...
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
            if (++iStepper >= GetNumRZCorner()) { iStepper = 0; }
          }
          while (chopped[iStepper]);
        }
        while (C < 0 && iStepper != iStarter);

        // Check triangle at B is pointing outward (an "ear").
        // Sign of z cross product determines...
        //
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
      //
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
          if (iSide == 0)   // startPhi
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = ixyz - GetNumRZCorner() + 2;
            }
          }
          else if (iSide == numSide - 1)   // endPhi
          {
            if (iCorner < GetNumRZCorner() - 1)
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
                faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
                faces_vec[iface][3] = -(ixyz + 2);
              }
            else
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = ixyz + GetNumRZCorner() + 1;
                faces_vec[iface][2] = ixyz + 2;
                faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
              }
          }
          else
          {
            if (iCorner < GetNumRZCorner() - 1)
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
                faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
                faces_vec[iface][3] = -(ixyz + 2);
              }
              else
              {
                faces_vec[iface][0] = ixyz + 1;
                faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
                faces_vec[iface][2] = ixyz + 2;
                faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
              }
            }
            ++iface;
            ++ixyz;
        }
        phi += dPhi;
      }

      // Last corners...

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
      nNodes = numSide * GetNumRZCorner();
      nFaces = numSide * GetNumRZCorner();;
      xyz = new double3[nNodes];
      faces_vec = new int4[nFaces];
      const G4double dPhi = (GetEndPhi() - GetStartPhi()) / numSide;
      G4double phi = GetStartPhi();
      G4int ixyz = 0, iface = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < GetNumRZCorner(); ++iCorner)
        {
          xyz[ixyz][0] = GetCorner(iCorner).r * std::cos(phi);
          xyz[ixyz][1] = GetCorner(iCorner).r * std::sin(phi);
          xyz[ixyz][2] = GetCorner(iCorner).z;

          if (iSide < numSide - 1)
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() + 2;
              faces_vec[iface][3] = -(ixyz + 2);
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
            }
          }
          else   // Last side joins ends...
          {
            if (iCorner < GetNumRZCorner() - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz + GetNumRZCorner() - nFaces + 1);
              faces_vec[iface][2] = ixyz + GetNumRZCorner() - nFaces + 2;
              faces_vec[iface][3] = -(ixyz + 2);
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = -(ixyz - nFaces + GetNumRZCorner() + 1);
              faces_vec[iface][2] = ixyz - nFaces + 2;
              faces_vec[iface][3] = -(ixyz - GetNumRZCorner() + 2);
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
      G4Exception("G4GenericPolycone::CreatePolyhedron()", "GeomSolids1002",
                  JustWarning, message);
      delete polyhedron;
      return 0;
    }
    else
    {
      return polyhedron;
    }
}

#endif  // G4GEOM_USE_USOLIDS
