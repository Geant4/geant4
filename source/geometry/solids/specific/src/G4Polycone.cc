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
// $Id: G4Polycone.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4Polycone.cc
//
// Implementation of a CSG polycone
//
// --------------------------------------------------------------------

#include "G4Polycone.hh"

#if !defined(G4GEOM_USE_UPOLYCONE)

#include "G4PolyconeSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "Randomize.hh"

#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"
#include "G4VPVParameterisation.hh"

using namespace CLHEP;

//
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

  G4int i;
  for (i=0; i<numZPlanes; i++)
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


//
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
  
  G4ReduciblePolygon *rz = new G4ReduciblePolygon( r, z, numRZ );
  
  Create( phiStart, phiTotal, rz );
  
  // Set original_parameters struct for consistency
  //
  
  G4bool convertible=SetOriginalParameters(rz);

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
                         G4ReduciblePolygon *rz    )
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
    || (!rz->RemoveRedundantVertices( kCarTolerance ))     ) 
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

  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if (phiTotal <= 0 || phiTotal > twopi-1E-10)
  {
    phiIsOpen = false;
    startPhi = 0;
    endPhi = twopi;
  }
  else
  {
    phiIsOpen = true;
    
    //
    // Convert phi into our convention
    //
    startPhi = phiStart;
    while( startPhi < 0 )    // Loop checking, 13.08.2015, G.Cosmo
      startPhi += twopi;
    
    endPhi = phiStart+phiTotal;
    while( endPhi < startPhi )    // Loop checking, 13.08.2015, G.Cosmo
      endPhi += twopi;
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
  G4PolyconeSideRZ *corner = corners,
                   *prev = corners + numCorner-1,
                   *nextNext;
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
  numFace = face-faces;
  
  //
  // Make enclosingCylinder
  //
  enclosingCylinder =
    new G4EnclosingCylinder( rz, phiIsOpen, phiStart, phiTotal );
}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4Polycone::G4Polycone( __void__& a )
  : G4VCSGfaceted(a), startPhi(0.),  endPhi(0.), phiIsOpen(false),
    numCorner(0), corners(0),original_parameters(0), 
    enclosingCylinder(0)
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
}


//
// Copy constructor
//
G4Polycone::G4Polycone( const G4Polycone &source )
  : G4VCSGfaceted( source )
{
  CopyStuff( source );
}


//
// Assignment operator
//
G4Polycone &G4Polycone::operator=( const G4Polycone &source )
{
  if (this == &source) return *this;
  
  G4VCSGfaceted::operator=( source );
  
  delete [] corners;
  if (original_parameters) delete original_parameters;
  
  delete enclosingCylinder;
  
  CopyStuff( source );
  
  return *this;
}


//
// CopyStuff
//
void G4Polycone::CopyStuff( const G4Polycone &source )
{
  //
  // Simple stuff
  //
  startPhi  = source.startPhi;
  endPhi    = source.endPhi;
  phiIsOpen  = source.phiIsOpen;
  numCorner  = source.numCorner;

  //
  // The corner array
  //
  corners = new G4PolyconeSideRZ[numCorner];
  
  G4PolyconeSideRZ  *corn = corners,
        *sourceCorn = source.corners;
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

  fRebuildPolyhedron = false;
  fpPolyhedron = 0;
}


//
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

  return 0;
}


//
// Inside
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
EInside G4Polycone::Inside( const G4ThreeVector &p ) const
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


//
// DistanceToIn
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
G4double G4Polycone::DistanceToIn( const G4ThreeVector &p,
                                   const G4ThreeVector &v ) const
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


//
// DistanceToIn
//
G4double G4Polycone::DistanceToIn( const G4ThreeVector &p ) const
{
  return G4VCSGfaceted::DistanceToIn(p);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

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

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

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

//
// ComputeDimensions
//
void G4Polycone::ComputeDimensions(       G4VPVParameterisation* p,
                                    const G4int n,
                                    const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}

//
// GetEntityType
//
G4GeometryType  G4Polycone::GetEntityType() const
{
  return G4String("G4Polycone");
}

//
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
  G4int oldprc = os.precision(16);
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
    for (i=0; i<numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Z_values[i] << "\n";
    }
    os << "              Tangent distances to inner surface (Rmin): \n";
    for (i=0; i<numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Rmin[i] << "\n";
    }
    os << "              Tangent distances to outer surface (Rmax): \n";
    for (i=0; i<numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << original_parameters->Rmax[i] << "\n";
    }
  
  os << "    number of RZ points: " << numCorner << "\n"
     << "              RZ values (corners): \n";
     for (i=0; i<numCorner; i++)
     {
       os << "                         "
          << corners[i].r << ", " << corners[i].z << "\n";
     }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


//
// GetPointOnCone
//
// Auxiliary method for Get Point On Surface
//
G4ThreeVector G4Polycone::GetPointOnCone(G4double fRmin1, G4double fRmax1,
                                         G4double fRmin2, G4double fRmax2,
                                         G4double zOne,   G4double zTwo,
                                         G4double& totArea) const
{ 
  // declare working variables
  //
  G4double Aone, Atwo, Afive, phi, zRand, fDPhi, cosu, sinu;
  G4double rRand1, rmin, rmax, chose, rone, rtwo, qone, qtwo;
  G4double fDz=(zTwo-zOne)/2., afDz=std::fabs(fDz);
  G4ThreeVector point, offset=G4ThreeVector(0.,0.,0.5*(zTwo+zOne));
  fDPhi = endPhi - startPhi;
  rone = (fRmax1-fRmax2)/(2.*fDz); 
  rtwo = (fRmin1-fRmin2)/(2.*fDz);
  if(fRmax1==fRmax2){qone=0.;}
  else
  {
    qone = fDz*(fRmax1+fRmax2)/(fRmax1-fRmax2);
  }
  if(fRmin1==fRmin2){qtwo=0.;}
  else
  {
    qtwo = fDz*(fRmin1+fRmin2)/(fRmin1-fRmin2);
  }
  Aone   = 0.5*fDPhi*(fRmax2 + fRmax1)*(sqr(fRmin1-fRmin2)+sqr(zTwo-zOne));       
  Atwo   = 0.5*fDPhi*(fRmin2 + fRmin1)*(sqr(fRmax1-fRmax2)+sqr(zTwo-zOne));
  Afive  = fDz*(fRmax1-fRmin1+fRmax2-fRmin2);
  totArea = Aone+Atwo+2.*Afive;
  
  phi  = G4RandFlat::shoot(startPhi,endPhi);
  cosu = std::cos(phi);
  sinu = std::sin(phi);


  if( (startPhi == 0) && (endPhi == twopi) ) { Afive = 0; }
  chose = G4RandFlat::shoot(0.,Aone+Atwo+2.*Afive);
  if( (chose >= 0) && (chose < Aone) )
  {
    if(fRmax1 != fRmax2)
    {
      zRand = G4RandFlat::shoot(-1.*afDz,afDz); 
      point = G4ThreeVector (rone*cosu*(qone-zRand),
                             rone*sinu*(qone-zRand), zRand);
    }
    else
    {
      point = G4ThreeVector(fRmax1*cosu, fRmax1*sinu,
                            G4RandFlat::shoot(-1.*afDz,afDz));
     
    }
  }
  else if(chose >= Aone && chose < Aone + Atwo)
  {
    if(fRmin1 != fRmin2)
    { 
      zRand = G4RandFlat::shoot(-1.*afDz,afDz); 
      point = G4ThreeVector (rtwo*cosu*(qtwo-zRand),
                             rtwo*sinu*(qtwo-zRand), zRand);
      
    }
    else
    {
      point = G4ThreeVector(fRmin1*cosu, fRmin1*sinu,
                            G4RandFlat::shoot(-1.*afDz,afDz));
    }
  }
  else if( (chose >= Aone + Atwo + Afive) && (chose < Aone + Atwo + 2.*Afive) )
  {
    zRand  = G4RandFlat::shoot(-1.*afDz,afDz);
    rmin   = fRmin2-((zRand-fDz)/(2.*fDz))*(fRmin1-fRmin2);
    rmax   = fRmax2-((zRand-fDz)/(2.*fDz))*(fRmax1-fRmax2);
    rRand1 = std::sqrt(G4RandFlat::shoot()*(sqr(rmax)-sqr(rmin))+sqr(rmin));
    point  = G4ThreeVector (rRand1*std::cos(startPhi),
                            rRand1*std::sin(startPhi), zRand);
  }
  else
  { 
    zRand  = G4RandFlat::shoot(-1.*afDz,afDz); 
    rmin   = fRmin2-((zRand-fDz)/(2.*fDz))*(fRmin1-fRmin2);
    rmax   = fRmax2-((zRand-fDz)/(2.*fDz))*(fRmax1-fRmax2);
    rRand1 = std::sqrt(G4RandFlat::shoot()*(sqr(rmax)-sqr(rmin))+sqr(rmin));
    point  = G4ThreeVector (rRand1*std::cos(endPhi),
                            rRand1*std::sin(endPhi), zRand);
   
  }
  
  return point+offset;
}


//
// GetPointOnTubs
//
// Auxiliary method for GetPoint On Surface
//
G4ThreeVector G4Polycone::GetPointOnTubs(G4double fRMin, G4double fRMax,
                                         G4double zOne,  G4double zTwo,
                                         G4double& totArea) const
{ 
  G4double xRand,yRand,zRand,phi,cosphi,sinphi,chose,
           aOne,aTwo,aFou,rRand,fDz,fSPhi,fDPhi;
  fDz = std::fabs(0.5*(zTwo-zOne));
  fSPhi = startPhi;
  fDPhi = endPhi-startPhi;
  
  aOne = 2.*fDz*fDPhi*fRMax;
  aTwo = 2.*fDz*fDPhi*fRMin;
  aFou = 2.*fDz*(fRMax-fRMin);
  totArea = aOne+aTwo+2.*aFou;
  phi    = G4RandFlat::shoot(startPhi,endPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);
  rRand  = fRMin + (fRMax-fRMin)*std::sqrt(G4RandFlat::shoot());
 
  if(startPhi == 0 && endPhi == twopi) 
    aFou = 0;
  
  chose  = G4RandFlat::shoot(0.,aOne+aTwo+2.*aFou);
  if( (chose >= 0) && (chose < aOne) )
  {
    xRand = fRMax*cosphi;
    yRand = fRMax*sinphi;
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector(xRand, yRand, zRand+0.5*(zTwo+zOne));
  }
  else if( (chose >= aOne) && (chose < aOne + aTwo) )
  {
    xRand = fRMin*cosphi;
    yRand = fRMin*sinphi;
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector(xRand, yRand, zRand+0.5*(zTwo+zOne));
  }
  else if( (chose >= aOne+aTwo) && (chose <aOne+aTwo+aFou) )
  {
    xRand = rRand*std::cos(fSPhi+fDPhi);
    yRand = rRand*std::sin(fSPhi+fDPhi);
    zRand = G4RandFlat::shoot(-1.*fDz,fDz);
    return G4ThreeVector(xRand, yRand, zRand+0.5*(zTwo+zOne));
  }

  // else

  xRand = rRand*std::cos(fSPhi+fDPhi);
  yRand = rRand*std::sin(fSPhi+fDPhi);
  zRand = G4RandFlat::shoot(-1.*fDz,fDz);
  return G4ThreeVector(xRand, yRand, zRand+0.5*(zTwo+zOne));
}


//
// GetPointOnRing
//
// Auxiliary method for GetPoint On Surface
//
G4ThreeVector G4Polycone::GetPointOnRing(G4double fRMin1, G4double fRMax1,
                                         G4double fRMin2,G4double fRMax2,
                                         G4double zOne) const
{
  G4double xRand,yRand,phi,cosphi,sinphi,rRand1,rRand2,A1,Atot,rCh;
  phi    = G4RandFlat::shoot(startPhi,endPhi);
  cosphi = std::cos(phi);
  sinphi = std::sin(phi);

  if(fRMin1==fRMin2)
  {
    rRand1 = fRMin1; A1=0.;
  }
  else
  {
    rRand1 = G4RandFlat::shoot(fRMin1,fRMin2);
    A1=std::fabs(fRMin2*fRMin2-fRMin1*fRMin1);
  }
  if(fRMax1==fRMax2)
  {
    rRand2=fRMax1; Atot=A1;
  }
  else
  {
    rRand2 = G4RandFlat::shoot(fRMax1,fRMax2);
    Atot   = A1+std::fabs(fRMax2*fRMax2-fRMax1*fRMax1);
  }
  rCh   = G4RandFlat::shoot(0.,Atot);
 
  if(rCh>A1) { rRand1=rRand2; }
  
  xRand = rRand1*cosphi;
  yRand = rRand1*sinphi;

  return G4ThreeVector(xRand, yRand, zOne);
}


//
// GetPointOnCut
//
// Auxiliary method for Get Point On Surface
//
G4ThreeVector G4Polycone::GetPointOnCut(G4double fRMin1, G4double fRMax1,
                                        G4double fRMin2, G4double fRMax2,
                                        G4double zOne,  G4double zTwo,
                                        G4double& totArea) const
{   if(zOne==zTwo)
    {
      return GetPointOnRing(fRMin1, fRMax1,fRMin2,fRMax2,zOne);
    }
    if( (fRMin1 == fRMin2) && (fRMax1 == fRMax2) )
    {
      return GetPointOnTubs(fRMin1, fRMax1,zOne,zTwo,totArea);
    }
    return GetPointOnCone(fRMin1,fRMax1,fRMin2,fRMax2,zOne,zTwo,totArea);
}


//
// GetPointOnSurface
//
G4ThreeVector G4Polycone::GetPointOnSurface() const
{
    G4double Area=0,totArea=0,Achose1=0,Achose2=0,phi,cosphi,sinphi,rRand;
    G4int i=0;
    G4int numPlanes = original_parameters->Num_z_planes;
  
    phi = G4RandFlat::shoot(startPhi,endPhi);
    cosphi = std::cos(phi);
    sinphi = std::sin(phi);

    rRand = original_parameters->Rmin[0] +
      ( (original_parameters->Rmax[0]-original_parameters->Rmin[0])
        * std::sqrt(G4RandFlat::shoot()) );

    std::vector<G4double> areas;       // (numPlanes+1);
    std::vector<G4ThreeVector> points; // (numPlanes-1);
  
    areas.push_back(pi*(sqr(original_parameters->Rmax[0])
                       -sqr(original_parameters->Rmin[0])));

    for(i=0; i<numPlanes-1; i++)
    {
      Area = (original_parameters->Rmin[i]+original_parameters->Rmin[i+1])
           * std::sqrt(sqr(original_parameters->Rmin[i]
                      -original_parameters->Rmin[i+1])+
                       sqr(original_parameters->Z_values[i+1]
                      -original_parameters->Z_values[i]));

      Area += (original_parameters->Rmax[i]+original_parameters->Rmax[i+1])
            * std::sqrt(sqr(original_parameters->Rmax[i]
                       -original_parameters->Rmax[i+1])+
                        sqr(original_parameters->Z_values[i+1]
                       -original_parameters->Z_values[i]));

      Area *= 0.5*(endPhi-startPhi);
    
      if(startPhi==0.&& endPhi == twopi)
      {
        Area += std::fabs(original_parameters->Z_values[i+1]
                         -original_parameters->Z_values[i])*
                         (original_parameters->Rmax[i]
                         +original_parameters->Rmax[i+1]
                         -original_parameters->Rmin[i]
                         -original_parameters->Rmin[i+1]);
      }
      areas.push_back(Area);
      totArea += Area;
    }
  
    areas.push_back(pi*(sqr(original_parameters->Rmax[numPlanes-1])-
                        sqr(original_parameters->Rmin[numPlanes-1])));
  
    totArea += (areas[0]+areas[numPlanes]);
    G4double chose = G4RandFlat::shoot(0.,totArea);

    if( (chose>=0.) && (chose<areas[0]) )
    {
      return G4ThreeVector(rRand*cosphi, rRand*sinphi,
                           original_parameters->Z_values[0]);
    }
  
    for (i=0; i<numPlanes-1; i++)
    {
      Achose1 += areas[i];
      Achose2 = (Achose1+areas[i+1]);
      if(chose>=Achose1 && chose<Achose2)
      {
        return GetPointOnCut(original_parameters->Rmin[i],
                             original_parameters->Rmax[i],
                             original_parameters->Rmin[i+1],
                             original_parameters->Rmax[i+1],
                             original_parameters->Z_values[i],
                             original_parameters->Z_values[i+1], Area);
      }
    }

    rRand = original_parameters->Rmin[numPlanes-1] +
      ( (original_parameters->Rmax[numPlanes-1]-original_parameters->Rmin[numPlanes-1])
        * std::sqrt(G4RandFlat::shoot()) );

    return G4ThreeVector(rRand*cosphi,rRand*sinphi,
                         original_parameters->Z_values[numPlanes-1]);  
 
}

//
// CreatePolyhedron
//
G4Polyhedron* G4Polycone::CreatePolyhedron() const
{ 
  //
  // This has to be fixed in visualization. Fake it for the moment.
  // 
  
    return new G4PolyhedronPcon( original_parameters->Start_angle,
                                 original_parameters->Opening_angle,
                                 original_parameters->Num_z_planes,
                                 original_parameters->Z_values,
                                 original_parameters->Rmin,
                                 original_parameters->Rmax );
}

G4bool  G4Polycone::SetOriginalParameters(G4ReduciblePolygon *rz)
{
  G4int numPlanes = (G4int)numCorner;  
  G4bool isConvertible=true;
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
  for (G4int i=0; i < numPlanes-2; i++)
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
      countPlanes++;
      icurr++;

      icurl=(icurl == 0)? numPlanes-1 : icurl-1;

      Rmin.push_back(corners[inextl].r);  
      Rmax.push_back(corners[inextr].r);
    }
    else  // Zright<Zleft
    {
      Z.push_back(Zright);  
      countPlanes++;

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
        icurr++;
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
        icurr++;
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
  countPlanes++;
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
  
   for(G4int j=0; j < countPlanes; j++)
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
  
    for(G4int j=0; j < numPlanes; j++)
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
