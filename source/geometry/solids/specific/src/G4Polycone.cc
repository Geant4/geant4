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
// Implementation of G4Polycone, a CSG polycone
//
// Author: David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4Polycone.hh"

#if !defined(G4GEOM_USE_UPOLYCONE)

#include "G4PolyconeSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "G4QuickRand.hh"

#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"
#include "G4VPVParameterisation.hh"

namespace
{
  G4Mutex surface_elementsMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

// Constructor (GEANT3 style parameters)
//
G4Polycone::G4Polycone( const G4String& name,
                              G4double phiStart,
                              G4double phiTotal,
                              G4int numZPlanes,
                        const G4double zPlane[],
                        const G4double rInner[],
                        const G4double rOuter[]  )
  : G4VCSGfaceted( name )
{
  //
  // Some historical ugliness
  //
  original_parameters = new G4PolyconeHistorical();

  original_parameters->Start_angle = phiStart;
  original_parameters->Opening_angle = phiTotal;
  original_parameters->Num_z_planes = numZPlanes;
  original_parameters->Z_values = new G4double[numZPlanes];
  original_parameters->Rmin = new G4double[numZPlanes];
  original_parameters->Rmax = new G4double[numZPlanes];

  for (G4int i=0; i<numZPlanes; ++i)
  {
    if(rInner[i]>rOuter[i])
    {
      DumpInfo();
      std::ostringstream message;
      message << "Cannot create a Polycone with rInner > rOuter for the same Z"
              << G4endl
              << "        rInner > rOuter for the same Z !" << G4endl
              << "        rMin[" << i << "] = " << rInner[i]
              << " -- rMax[" << i << "] = " << rOuter[i];
      G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002",
                  FatalErrorInArgument, message);
    }
    if (( i < numZPlanes-1) && ( zPlane[i] == zPlane[i+1] ))
    {
      if( (rInner[i]   > rOuter[i+1])
        ||(rInner[i+1] > rOuter[i])   )
      {
        DumpInfo();
        std::ostringstream message;
        message << "Cannot create a Polycone with no contiguous segments."
                << G4endl
                << "        Segments are not contiguous !" << G4endl
                << "        rMin[" << i << "] = " << rInner[i]
                << " -- rMax[" << i+1 << "] = " << rOuter[i+1] << G4endl
                << "        rMin[" << i+1 << "] = " << rInner[i+1]
                << " -- rMax[" << i << "] = " << rOuter[i];
        G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002",
                    FatalErrorInArgument, message);
      }
    }
    original_parameters->Z_values[i] = zPlane[i];
    original_parameters->Rmin[i] = rInner[i];
    original_parameters->Rmax[i] = rOuter[i];
  }

  //
  // Build RZ polygon using special PCON/PGON GEANT3 constructor
  //
  G4ReduciblePolygon *rz =
    new G4ReduciblePolygon( rInner, rOuter, zPlane, numZPlanes );

  //
  // Do the real work
  //
  Create( phiStart, phiTotal, rz );

  delete rz;
}

// Constructor (generic parameters)
//
G4Polycone::G4Polycone( const G4String& name,
                              G4double phiStart,
                              G4double phiTotal,
                              G4int    numRZ,
                        const G4double r[],
                        const G4double z[]   )
  : G4VCSGfaceted( name )
{

  G4ReduciblePolygon* rz = new G4ReduciblePolygon( r, z, numRZ );

  Create( phiStart, phiTotal, rz );

  // Set original_parameters struct for consistency
  //

  G4bool convertible = SetOriginalParameters(rz);

  if(!convertible)
  {
    std::ostringstream message;
    message << "Polycone " << GetName() << "cannot be converted" << G4endl
            << "to Polycone with (Rmin,Rmaz,Z) parameters!";
    G4Exception("G4Polycone::G4Polycone()", "GeomSolids0002",
                FatalException, message, "Use G4GenericPolycone instead!");
  }
  else
  {
    G4cout << "INFO: Converting polycone " << GetName() << G4endl
           << "to optimized polycone with (Rmin,Rmaz,Z) parameters !"
           << G4endl;
  }
  delete rz;
}

// Create
//
// Generic create routine, called by each constructor after
// conversion of arguments
//
void G4Polycone::Create( G4double phiStart,
                         G4double phiTotal,
                         G4ReduciblePolygon* rz )
{
  //
  // Perform checks of rz values
  //
  if (rz->Amin() < 0.0)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        All R values must be >= 0 !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  G4double rzArea = rz->Area();
  if (rzArea < -kCarTolerance)
  {
    rz->ReverseOrder();
  }
  else if (rzArea < kCarTolerance)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        R/Z cross section is zero or near zero: " << rzArea;
    G4Exception("G4Polycone::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  if ( (!rz->RemoveDuplicateVertices( kCarTolerance ))
    || (!rz->RemoveRedundantVertices( kCarTolerance )) )
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        Too few unique R/Z values !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  if (rz->CrossesItself(1/kInfinity))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        R/Z segments cross !";
    G4Exception("G4Polycone::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  numCorner = rz->NumVertices();

  startPhi = phiStart;
  while( startPhi < 0. )    // Loop checking, 13.08.2015, G.Cosmo
    startPhi += twopi;
  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if ( (phiTotal <= 0) || (phiTotal > twopi*(1-DBL_EPSILON)) )
  {
    phiIsOpen = false;
    startPhi = 0.;
    endPhi = twopi;
  }
  else
  {
    phiIsOpen = true;
    endPhi = startPhi + phiTotal;
  }

  //
  // Allocate corner array.
  //
  corners = new G4PolyconeSideRZ[numCorner];

  //
  // Copy corners
  //
  G4ReduciblePolygonIterator iterRZ(rz);

  G4PolyconeSideRZ *next = corners;
  iterRZ.Begin();
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    next->r = iterRZ.GetA();
    next->z = iterRZ.GetB();
  } while( ++next, iterRZ.Next() );

  //
  // Allocate face pointer array
  //
  numFace = phiIsOpen ? numCorner+2 : numCorner;
  faces = new G4VCSGface*[numFace];

  //
  // Construct conical faces
  //
  // But! Don't construct a face if both points are at zero radius!
  //
  G4PolyconeSideRZ* corner = corners,
                  * prev = corners + numCorner-1,
                  * nextNext;
  G4VCSGface  **face = faces;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    next = corner+1;
    if (next >= corners+numCorner) next = corners;
    nextNext = next+1;
    if (nextNext >= corners+numCorner) nextNext = corners;

    if (corner->r < 1/kInfinity && next->r < 1/kInfinity) continue;

    //
    // We must decide here if we can dare declare one of our faces
    // as having a "valid" normal (i.e. allBehind = true). This
    // is never possible if the face faces "inward" in r.
    //
    G4bool allBehind;
    if (corner->z > next->z)
    {
      allBehind = false;
    }
    else
    {
      //
      // Otherwise, it is only true if the line passing
      // through the two points of the segment do not
      // split the r/z cross section
      //
      allBehind = !rz->BisectedBy( corner->r, corner->z,
                 next->r, next->z, kCarTolerance );
    }

    *face++ = new G4PolyconeSide( prev, corner, next, nextNext,
                startPhi, endPhi-startPhi, phiIsOpen, allBehind );
  } while( prev=corner, corner=next, corner > corners );

  if (phiIsOpen)
  {
    //
    // Construct phi open edges
    //
    *face++ = new G4PolyPhiFace( rz, startPhi, 0, endPhi  );
    *face++ = new G4PolyPhiFace( rz, endPhi,   0, startPhi );
  }

  //
  // We might have dropped a face or two: recalculate numFace
  //
  numFace = (G4int)(face-faces);

  //
  // Make enclosingCylinder
  //
  enclosingCylinder =
    new G4EnclosingCylinder( rz, phiIsOpen, phiStart, phiTotal );
}

// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Polycone::G4Polycone( __void__& a )
  : G4VCSGfaceted(a), startPhi(0.),  endPhi(0.), numCorner(0)
{
}


//
// Destructor
//
G4Polycone::~G4Polycone()
{
  delete [] corners;
  delete original_parameters;
  delete enclosingCylinder;
  delete fElements;
  delete fpPolyhedron;
  corners = nullptr;
  original_parameters = nullptr;
  enclosingCylinder = nullptr;
  fElements = nullptr;
  fpPolyhedron = nullptr;
}

// Copy constructor
//
G4Polycone::G4Polycone( const G4Polycone& source )
  : G4VCSGfaceted( source )
{
  CopyStuff( source );
}

// Assignment operator
//
G4Polycone &G4Polycone::operator=( const G4Polycone& source )
{
  if (this == &source) return *this;

  G4VCSGfaceted::operator=( source );

  delete [] corners;
  if (original_parameters) delete original_parameters;

  delete enclosingCylinder;

  CopyStuff( source );

  return *this;
}

// CopyStuff
//
void G4Polycone::CopyStuff( const G4Polycone& source )
{
  //
  // Simple stuff
  //
  startPhi  = source.startPhi;
  endPhi    = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  numCorner = source.numCorner;

  //
  // The corner array
  //
  corners = new G4PolyconeSideRZ[numCorner];

  G4PolyconeSideRZ* corn = corners,
                  * sourceCorn = source.corners;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    *corn = *sourceCorn;
  } while( ++sourceCorn, ++corn < corners+numCorner );

  //
  // Original parameters
  //
  if (source.original_parameters)
  {
    original_parameters =
      new G4PolyconeHistorical( *source.original_parameters );
  }

  //
  // Enclosing cylinder
  //
  enclosingCylinder = new G4EnclosingCylinder( *source.enclosingCylinder );

  //
  // Surface elements
  //
  delete fElements;
  fElements = nullptr;

  //
  // Polyhedron
  //
  fRebuildPolyhedron = false;
  delete fpPolyhedron;
  fpPolyhedron = nullptr;
}

// Reset
//
G4bool G4Polycone::Reset()
{
  //
  // Clear old setup
  //
  G4VCSGfaceted::DeleteStuff();
  delete [] corners;
  delete enclosingCylinder;
  delete fElements;
  corners = nullptr;
  fElements = nullptr;
  enclosingCylinder = nullptr;

  //
  // Rebuild polycone
  //
  G4ReduciblePolygon *rz =
    new G4ReduciblePolygon( original_parameters->Rmin,
                            original_parameters->Rmax,
                            original_parameters->Z_values,
                            original_parameters->Num_z_planes );
  Create( original_parameters->Start_angle,
          original_parameters->Opening_angle, rz );
  delete rz;

  return false;
}

// Inside
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
EInside G4Polycone::Inside( const G4ThreeVector& p ) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->MustBeOutside(p)) return kOutside;

  //
  // Long answer
  //
  return G4VCSGfaceted::Inside(p);
}

// DistanceToIn
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
G4double G4Polycone::DistanceToIn( const G4ThreeVector& p,
                                   const G4ThreeVector& v ) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->ShouldMiss(p,v))
    return kInfinity;

  //
  // Long answer
  //
  return G4VCSGfaceted::DistanceToIn( p, v );
}

// DistanceToIn
//
G4double G4Polycone::DistanceToIn( const G4ThreeVector& p ) const
{
  return G4VCSGfaceted::DistanceToIn(p);
}

// Get bounding box
//
void G4Polycone::BoundingLimits(G4ThreeVector& pMin,
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
    G4Exception("G4Polycone::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

// Calculate extent under transform and specified limit
//
G4bool G4Polycone::CalculateExtent(const EAxis pAxis,
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
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
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
    G4Exception("G4Polycone::CalculateExtent()",
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
  G4int ntria = (G4int)triangles.size()/3;
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

// ComputeDimensions
//
void G4Polycone::ComputeDimensions(       G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

// GetEntityType
//
G4GeometryType  G4Polycone::GetEntityType() const
{
  return G4String("G4Polycone");
}

// Make a clone of the object
//
G4VSolid* G4Polycone::Clone() const
{
  return new G4Polycone(*this);
}

//
// Stream object contents to an output stream
//
std::ostream& G4Polycone::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Polycone\n"
     << " Parameters: \n"
     << "    starting phi angle : " << startPhi/degree << " degrees \n"
     << "    ending phi angle   : " << endPhi/degree << " degrees \n";
  G4int i=0;

    G4int numPlanes = original_parameters->Num_z_planes;
    os << "    number of Z planes: " << numPlanes << "\n"
       << "              Z values: \n";
    for (i=0; i<numPlanes; ++i)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Z_values[i] << "\n";
    }
    os << "              Tangent distances to inner surface (Rmin): \n";
    for (i=0; i<numPlanes; ++i)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Rmin[i] << "\n";
    }
    os << "              Tangent distances to outer surface (Rmax): \n";
    for (i=0; i<numPlanes; ++i)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Rmax[i] << "\n";
    }

  os << "    number of RZ points: " << numCorner << "\n"
     << "              RZ values (corners): \n";
     for (i=0; i<numCorner; ++i)
     {
       os << "                         "
          << corners[i].r << ", " << corners[i].z << "\n";
     }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//////////////////////////////////////////////////////////////////////////
//
// Return volume

G4double G4Polycone::GetCubicVolume()
{
  if (fCubicVolume == 0.)
  {
    G4double total = 0.;
    G4int nrz = GetNumRZCorner();
    G4PolyconeSideRZ a = GetCorner(nrz - 1);
    for (G4int i=0; i<nrz; ++i)
    {
      G4PolyconeSideRZ b = GetCorner(i);
      total += (b.r*b.r + b.r*a.r + a.r*a.r)*(b.z - a.z);
      a = b;
    }
    fCubicVolume = std::abs(total)*(GetEndPhi() - GetStartPhi())/6.;
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Return surface area

G4double G4Polycone::GetSurfaceArea()
{
  if (fSurfaceArea == 0.)
  {
    // phi cut area
    G4int nrz = GetNumRZCorner();
    G4double scut = 0.;
    if (IsOpen())
    {
      G4PolyconeSideRZ a = GetCorner(nrz - 1);
      for (G4int i=0; i<nrz; ++i)
      {
        G4PolyconeSideRZ b = GetCorner(i);
        scut += a.r*b.z - a.z*b.r;
        a = b;
      }
      scut = std::abs(scut);
    }
    // lateral surface area
    G4double slat = 0;
    G4PolyconeSideRZ a = GetCorner(nrz - 1);
    for (G4int i=0; i<nrz; ++i)
    {
      G4PolyconeSideRZ b = GetCorner(i);
      G4double h = std::sqrt((b.r - a.r)*(b.r - a.r) + (b.z - a.z)*(b.z - a.z));
      slat += (b.r + a.r)*h;
      a = b;
    }
    slat *= (GetEndPhi() - GetStartPhi())/2.;
    fSurfaceArea = scut + slat;
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Set vector of surface elements, auxiliary method for sampling
// random points on surface

void G4Polycone::SetSurfaceElements() const
{
  fElements = new std::vector<G4Polycone::surface_element>;
  G4double total = 0.;
  G4int nrz = GetNumRZCorner();

  // set lateral surface elements
  G4double dphi = GetEndPhi() - GetStartPhi();
  G4int ia = nrz - 1;
  for (G4int ib=0; ib<nrz; ++ib)
  {
    G4PolyconeSideRZ a = GetCorner(ia);
    G4PolyconeSideRZ b = GetCorner(ib);
    G4Polycone::surface_element selem;
    selem.i0 = ia;
    selem.i1 = ib;
    selem.i2 = -1;
    ia = ib;
    if (a.r == 0. && b.r == 0.) continue;
    G4double h = std::sqrt((b.r - a.r)*(b.r - a.r) + (b.z - a.z)*(b.z - a.z));
    total += 0.5*dphi*(b.r + a.r)*h;
    selem.area = total;
    fElements->push_back(selem);
  }

  // set elements for phi cuts
  if (IsOpen())
  {
    G4TwoVectorList contourRZ;
    std::vector<G4int> triangles;
    for (G4int i=0; i<nrz; ++i)
    {
      G4PolyconeSideRZ corner = GetCorner(i);
      contourRZ.push_back(G4TwoVector(corner.r, corner.z));
    }
    G4GeomTools::TriangulatePolygon(contourRZ, triangles);
    G4int ntria = (G4int)triangles.size();
    for (G4int i=0; i<ntria; i+=3)
    {
      G4Polycone::surface_element selem;
      selem.i0 = triangles[i];
      selem.i1 = triangles[i+1];
      selem.i2 = triangles[i+2];
      G4PolyconeSideRZ a = GetCorner(selem.i0);
      G4PolyconeSideRZ b = GetCorner(selem.i1);
      G4PolyconeSideRZ c = GetCorner(selem.i2);
      G4double stria =
        std::abs(G4GeomTools::TriangleArea(a.r, a.z, b.r, b.z, c.r, c.z));
      total += stria;
      selem.area = total;
      fElements->push_back(selem); // start phi
      total += stria;
      selem.area = total;
      selem.i0 += nrz;
      fElements->push_back(selem); // end phi
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Generate random point on surface

G4ThreeVector G4Polycone::GetPointOnSurface() const
{
  // Set surface elements
  if (!fElements)
  {
    G4AutoLock l(&surface_elementsMutex);
    SetSurfaceElements();
    l.unlock();
  }

  // Select surface element
  G4Polycone::surface_element selem;
  selem = fElements->back();
  G4double select = selem.area*G4QuickRand();
  auto it = std::lower_bound(fElements->begin(), fElements->end(), select,
                             [](const G4Polycone::surface_element& x, G4double val)
                             -> G4bool { return x.area < val; });

  // Generate random point
  G4double r = 0, z = 0, phi = 0;
  G4double u = G4QuickRand();
  G4double v = G4QuickRand();
  G4int i0 = (*it).i0;
  G4int i1 = (*it).i1;
  G4int i2 = (*it).i2;
  if (i2 < 0) // lateral surface
  {
    G4PolyconeSideRZ p0 = GetCorner(i0);
    G4PolyconeSideRZ p1 = GetCorner(i1);
    if (p1.r < p0.r)
    {
      p0 = GetCorner(i1);
      p1 = GetCorner(i0);
    }
    if (p1.r - p0.r < kCarTolerance) // cylindrical surface
    {
      r = (p1.r - p0.r)*u + p0.r;
      z = (p1.z - p0.z)*u + p0.z;
    }
    else // conical surface
    {
      r = std::sqrt(p1.r*p1.r*u + p0.r*p0.r*(1. - u));
      z = p0.z + (p1.z - p0.z)*(r - p0.r)/(p1.r - p0.r);
    }
    phi = (GetEndPhi() - GetStartPhi())*v + GetStartPhi();
  }
  else // phi cut
  {
    G4int nrz = GetNumRZCorner();
    phi = (i0 < nrz) ? GetStartPhi() : GetEndPhi();
    if (i0 >= nrz) { i0 -= nrz; }
    G4PolyconeSideRZ p0 = GetCorner(i0);
    G4PolyconeSideRZ p1 = GetCorner(i1);
    G4PolyconeSideRZ p2 = GetCorner(i2);
    if (u + v > 1.) { u = 1. - u; v = 1. - v; }
    r = (p1.r - p0.r)*u +  (p2.r - p0.r)*v + p0.r;
    z = (p1.z - p0.z)*u +  (p2.z - p0.z)*v + p0.z;
  }
  return G4ThreeVector(r*std::cos(phi), r*std::sin(phi), z);
}

//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron

G4Polyhedron* G4Polycone::CreatePolyhedron() const
{
  std::vector<G4TwoVector> rz(numCorner);
  for (G4int i = 0; i < numCorner; ++i)
    rz[i].set(corners[i].r, corners[i].z);
  return new G4PolyhedronPcon(startPhi, endPhi - startPhi, rz);
}

// SetOriginalParameters
//
G4bool  G4Polycone::SetOriginalParameters(G4ReduciblePolygon* rz)
{
  G4int numPlanes = numCorner;
  G4bool isConvertible = true;
  G4double Zmax=rz->Bmax();
  rz->StartWithZMin();

  // Prepare vectors for storage
  //
  std::vector<G4double> Z;
  std::vector<G4double> Rmin;
  std::vector<G4double> Rmax;

  G4int countPlanes=1;
  G4int icurr=0;
  G4int icurl=0;

  // first plane Z=Z[0]
  //
  Z.push_back(corners[0].z);
  G4double Zprev=Z[0];
  if (Zprev == corners[1].z)
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back (corners[1].r);icurr=1;
  }
  else if (Zprev == corners[numPlanes-1].z)
  {
    Rmin.push_back(corners[numPlanes-1].r);
    Rmax.push_back (corners[0].r);
    icurl=numPlanes-1;
  }
  else
  {
    Rmin.push_back(corners[0].r);
    Rmax.push_back (corners[0].r);
  }

  // next planes until last
  //
  G4int inextr=0, inextl=0;
  for (G4int i=0; i < numPlanes-2; ++i)
  {
    inextr=1+icurr;
    inextl=(icurl <= 0)? numPlanes-1 : icurl-1;

    if((corners[inextr].z >= Zmax) & (corners[inextl].z >= Zmax))  { break; }

    G4double Zleft = corners[inextl].z;
    G4double Zright = corners[inextr].z;
    if(Zright > Zleft)  // Next plane will be Zleft
    {
      Z.push_back(Zleft);
      countPlanes++;
      G4double difZr=corners[inextr].z - corners[icurr].z;
      G4double difZl=corners[inextl].z - corners[icurl].z;

      if(std::fabs(difZl) < kCarTolerance)
      {
        if(std::fabs(difZr) < kCarTolerance)
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[inextl].r);
          Rmax.push_back(corners[icurr].r + (Zleft-corners[icurr].z)/difZr
                                *(corners[inextr].r - corners[icurr].r));
        }
      }
      else if (difZl >= kCarTolerance)
      {
        if(std::fabs(difZr) < kCarTolerance)
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[icurl].r);
          Rmax.push_back(corners[icurr].r + (Zleft-corners[icurr].z)/difZr
                                *(corners[inextr].r - corners[icurr].r));
        }
      }
      else
      {
        isConvertible=false; break;
      }
      icurl=(icurl == 0)? numPlanes-1 : icurl-1;
    }
    else if(std::fabs(Zright-Zleft)<kCarTolerance)  // Zright=Zleft
    {
      Z.push_back(Zleft);
      ++countPlanes;
      ++icurr;

      icurl=(icurl == 0)? numPlanes-1 : icurl-1;

      Rmin.push_back(corners[inextl].r);
      Rmax.push_back(corners[inextr].r);
    }
    else  // Zright<Zleft
    {
      Z.push_back(Zright);
      ++countPlanes;

      G4double difZr=corners[inextr].z - corners[icurr].z;
      G4double difZl=corners[inextl].z - corners[icurl].z;
      if(std::fabs(difZr) < kCarTolerance)
      {
        if(std::fabs(difZl) < kCarTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back(corners[icurr].r);
        }
        else
        {
          Rmin.push_back(corners[icurl].r + (Zright-corners[icurl].z)/difZl
                                *(corners[inextl].r - corners[icurl].r));
          Rmax.push_back(corners[inextr].r);
        }
        ++icurr;
      }           // plate
      else if (difZr >= kCarTolerance)
      {
        if(std::fabs(difZl) < kCarTolerance)
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back (corners[icurr].r);
        }
        else
        {
          Rmax.push_back(corners[inextr].r);
          Rmin.push_back (corners[icurl].r+(Zright-corners[icurl].z)/difZl
                                  * (corners[inextl].r - corners[icurl].r));
        }
        ++icurr;
      }
      else
      {
        isConvertible=false; break;
      }
    }
  }   // end for loop

  // last plane Z=Zmax
  //
  Z.push_back(Zmax);
  ++countPlanes;
  inextr=1+icurr;
  inextl=(icurl <= 0)? numPlanes-1 : icurl-1;

  Rmax.push_back(corners[inextr].r);
  Rmin.push_back(corners[inextl].r);

  // Set original parameters Rmin,Rmax,Z
  //
  if(isConvertible)
  {
   original_parameters = new G4PolyconeHistorical;
   original_parameters->Z_values = new G4double[countPlanes];
   original_parameters->Rmin = new G4double[countPlanes];
   original_parameters->Rmax = new G4double[countPlanes];

   for(G4int j=0; j < countPlanes; ++j)
   {
     original_parameters->Z_values[j] = Z[j];
     original_parameters->Rmax[j] = Rmax[j];
     original_parameters->Rmin[j] = Rmin[j];
   }
   original_parameters->Start_angle = startPhi;
   original_parameters->Opening_angle = endPhi-startPhi;
   original_parameters->Num_z_planes = countPlanes;

  }
  else  // Set parameters(r,z) with Rmin==0 as convention
  {
#ifdef G4SPECSDEBUG
    std::ostringstream message;
    message << "Polycone " << GetName() << G4endl
            << "cannot be converted to Polycone with (Rmin,Rmaz,Z) parameters!";
    G4Exception("G4Polycone::SetOriginalParameters()", "GeomSolids0002",
                JustWarning, message);
#endif
    original_parameters = new G4PolyconeHistorical;
    original_parameters->Z_values = new G4double[numPlanes];
    original_parameters->Rmin = new G4double[numPlanes];
    original_parameters->Rmax = new G4double[numPlanes];

    for(G4int j=0; j < numPlanes; ++j)
    {
      original_parameters->Z_values[j] = corners[j].z;
      original_parameters->Rmax[j] = corners[j].r;
      original_parameters->Rmin[j] = 0.0;
    }
    original_parameters->Start_angle = startPhi;
    original_parameters->Opening_angle = endPhi-startPhi;
    original_parameters->Num_z_planes = numPlanes;
  }
  return isConvertible;
}

#endif
