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
// $Id: G4Polyhedra.cc 104316 2017-05-24 13:04:23Z gcosmo $
//
// 
// --------------------------------------------------------------------
// GEANT 4 class source file
//
//
// G4Polyhedra.cc
//
// Implementation of a CSG polyhedra, as an inherited class of G4VCSGfaceted.
//
// To be done:
//    * Cracks: there are probably small cracks in the seams between the
//      phi face (G4PolyPhiFace) and sides (G4PolyhedraSide) that are not
//      entirely leakproof. Also, I am not sure all vertices are leak proof.
//    * Many optimizations are possible, but not implemented.
//    * Visualization needs to be updated outside of this routine.
//
// Utility classes:
//    * G4EnclosingCylinder: I decided a quick check of geometry would be a
//      good idea (for CPU speed). If the quick check fails, the regular
//      full-blown G4VCSGfaceted version is invoked.
//    * G4ReduciblePolygon: Really meant as a check of input parameters,
//      this utility class also "converts" the GEANT3-like PGON/PCON
//      arguments into the newer ones.
// Both these classes are implemented outside this file because they are
// shared with G4Polycone.
//
// --------------------------------------------------------------------

#include "G4Polyhedra.hh"

#if !defined(G4GEOM_USE_UPOLYHEDRA)

#include "G4PolyhedraSide.hh"
#include "G4PolyPhiFace.hh"

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4BoundingEnvelope.hh"

#include "Randomize.hh"

#include "G4EnclosingCylinder.hh"
#include "G4ReduciblePolygon.hh"
#include "G4VPVParameterisation.hh"

#include <sstream>

using namespace CLHEP;

//
// Constructor (GEANT3 style parameters)
//
// GEANT3 PGON radii are specified in the distance to the norm of each face.
//  
G4Polyhedra::G4Polyhedra( const G4String& name, 
                                G4double phiStart,
                                G4double thePhiTotal,
                                G4int theNumSide,  
                                G4int numZPlanes,
                          const G4double zPlane[],
                          const G4double rInner[],
                          const G4double rOuter[]  )
  : G4VCSGfaceted( name ), genericPgon(false)
{
  if (theNumSide <= 0)
  {
    std::ostringstream message;
    message << "Solid must have at least one side - " << GetName() << G4endl
            << "        No sides specified !";
    G4Exception("G4Polyhedra::G4Polyhedra()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  //
  // Calculate conversion factor from G3 radius to G4 radius
  //
  G4double phiTotal = thePhiTotal;
  if ( (phiTotal <=0) || (phiTotal >= twopi*(1-DBL_EPSILON)) )
    { phiTotal = twopi; }
  G4double convertRad = std::cos(0.5*phiTotal/theNumSide);

  //
  // Some historical stuff
  //
  original_parameters = new G4PolyhedraHistorical;
  
  original_parameters->numSide = theNumSide;
  original_parameters->Start_angle = phiStart;
  original_parameters->Opening_angle = phiTotal;
  original_parameters->Num_z_planes = numZPlanes;
  original_parameters->Z_values = new G4double[numZPlanes];
  original_parameters->Rmin = new G4double[numZPlanes];
  original_parameters->Rmax = new G4double[numZPlanes];

  G4int i;
  for (i=0; i<numZPlanes; i++)
  {
    if (( i < numZPlanes-1) && ( zPlane[i] == zPlane[i+1] ))
    {
      if( (rInner[i]   > rOuter[i+1])
        ||(rInner[i+1] > rOuter[i])   )
      {
        DumpInfo();
        std::ostringstream message;
        message << "Cannot create a Polyhedra with no contiguous segments."
                << G4endl
                << "        Segments are not contiguous !" << G4endl
                << "        rMin[" << i << "] = " << rInner[i]
                << " -- rMax[" << i+1 << "] = " << rOuter[i+1] << G4endl
                << "        rMin[" << i+1 << "] = " << rInner[i+1]
                << " -- rMax[" << i << "] = " << rOuter[i];
        G4Exception("G4Polyhedra::G4Polyhedra()", "GeomSolids0002",
                    FatalErrorInArgument, message);
      }
    }
    original_parameters->Z_values[i] = zPlane[i];
    original_parameters->Rmin[i] = rInner[i]/convertRad;
    original_parameters->Rmax[i] = rOuter[i]/convertRad;
  }
  
  
  //
  // Build RZ polygon using special PCON/PGON GEANT3 constructor
  //
  G4ReduciblePolygon *rz =
    new G4ReduciblePolygon( rInner, rOuter, zPlane, numZPlanes );
  rz->ScaleA( 1/convertRad );
  
  //
  // Do the real work
  //
  Create( phiStart, phiTotal, theNumSide, rz );
  
  delete rz;
}


//
// Constructor (generic parameters)
//
G4Polyhedra::G4Polyhedra( const G4String& name, 
                                G4double phiStart,
                                G4double phiTotal,
                                G4int    theNumSide,  
                                G4int    numRZ,
                          const G4double r[],
                          const G4double z[]   )
  : G4VCSGfaceted( name ), genericPgon(true)
{ 
  if (theNumSide <= 0)
  {
    std::ostringstream message;
    message << "Solid must have at least one side - " << GetName() << G4endl
            << "        No sides specified !";
    G4Exception("G4Polyhedra::G4Polyhedra()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  G4ReduciblePolygon *rz = new G4ReduciblePolygon( r, z, numRZ );
  
  Create( phiStart, phiTotal, theNumSide, rz );
  
  // Set original_parameters struct for consistency
  //
  SetOriginalParameters(rz);
   
  delete rz;
}


//
// Create
//
// Generic create routine, called by each constructor
// after conversion of arguments
//
void G4Polyhedra::Create( G4double phiStart,
                          G4double phiTotal,
                          G4int    theNumSide,  
                          G4ReduciblePolygon *rz  )
{
  //
  // Perform checks of rz values
  //
  if (rz->Amin() < 0.0)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        All R values must be >= 0 !";
    G4Exception("G4Polyhedra::Create()", "GeomSolids0002",
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
    G4Exception("G4Polyhedra::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }
    
  if ( (!rz->RemoveDuplicateVertices( kCarTolerance ))
    || (!rz->RemoveRedundantVertices( kCarTolerance )) ) 
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        Too few unique R/Z values !";
    G4Exception("G4Polyhedra::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  if (rz->CrossesItself( 1/kInfinity )) 
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << G4endl
            << "        R/Z segments cross !";
    G4Exception("G4Polyhedra::Create()", "GeomSolids0002",
                FatalErrorInArgument, message);
  }

  numCorner = rz->NumVertices();


  startPhi = phiStart;
  while( startPhi < 0 )    // Loop checking, 13.08.2015, G.Cosmo
    startPhi += twopi;
  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if ( (phiTotal <= 0) || (phiTotal > twopi*(1-DBL_EPSILON)) )
  {
    phiIsOpen = false;
    endPhi = phiStart+twopi;
  }
  else
  {
    phiIsOpen = true;
    
    //
    // Convert phi into our convention
    //
    endPhi = phiStart+phiTotal;
    while( endPhi < startPhi )    // Loop checking, 13.08.2015, G.Cosmo
      endPhi += twopi;
  }
  
  //
  // Save number sides
  //
  numSide = theNumSide;
  
  //
  // Allocate corner array. 
  //
  corners = new G4PolyhedraSideRZ[numCorner];

  //
  // Copy corners
  //
  G4ReduciblePolygonIterator iterRZ(rz);
  
  G4PolyhedraSideRZ *next = corners;
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
  // Construct side faces
  //
  // To do so properly, we need to keep track of four successive RZ
  // corners.
  //
  // But! Don't construct a face if both points are at zero radius!
  //
  G4PolyhedraSideRZ *corner = corners,
                    *prev = corners + numCorner-1,
                    *nextNext;
  G4VCSGface   **face = faces;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    next = corner+1;
    if (next >= corners+numCorner) next = corners;
    nextNext = next+1;
    if (nextNext >= corners+numCorner) nextNext = corners;
    
    if (corner->r < 1/kInfinity && next->r < 1/kInfinity) continue;
/*
    // We must decide here if we can dare declare one of our faces
    // as having a "valid" normal (i.e. allBehind = true). This
    // is never possible if the face faces "inward" in r *unless*
    // we have only one side
    //
    G4bool allBehind;
    if ((corner->z > next->z) && (numSide > 1))
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
*/
    *face++ = new G4PolyhedraSide( prev, corner, next, nextNext,
                 numSide, startPhi, endPhi-startPhi, phiIsOpen );
  } while( prev=corner, corner=next, corner > corners );
  
  if (phiIsOpen)
  {
    //
    // Construct phi open edges
    //
    *face++ = new G4PolyPhiFace( rz, startPhi, phiTotal/numSide, endPhi );
    *face++ = new G4PolyPhiFace( rz, endPhi,   phiTotal/numSide, startPhi );
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
G4Polyhedra::G4Polyhedra( __void__& a )
  : G4VCSGfaceted(a), numSide(0), startPhi(0.), endPhi(0.),
    phiIsOpen(false), genericPgon(false), numCorner(0), corners(0),
    original_parameters(0), enclosingCylinder(0)
{
}


//
// Destructor
//
G4Polyhedra::~G4Polyhedra()
{
  delete [] corners;
  if (original_parameters) delete original_parameters;
  
  delete enclosingCylinder;
}


//
// Copy constructor
//
G4Polyhedra::G4Polyhedra( const G4Polyhedra &source )
  : G4VCSGfaceted( source )
{
  CopyStuff( source );
}


//
// Assignment operator
//
G4Polyhedra &G4Polyhedra::operator=( const G4Polyhedra &source )
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
void G4Polyhedra::CopyStuff( const G4Polyhedra &source )
{
  //
  // Simple stuff
  //
  numSide    = source.numSide;
  startPhi   = source.startPhi;
  endPhi     = source.endPhi;
  phiIsOpen  = source.phiIsOpen;
  numCorner  = source.numCorner;
  genericPgon= source.genericPgon;

  //
  // The corner array
  //
  corners = new G4PolyhedraSideRZ[numCorner];
  
  G4PolyhedraSideRZ  *corn = corners,
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
      new G4PolyhedraHistorical( *source.original_parameters );
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
// Recalculates and reshapes the solid, given pre-assigned scaled
// original_parameters.
//
G4bool G4Polyhedra::Reset()
{
  if (genericPgon)
  {
    std::ostringstream message;
    message << "Solid " << GetName() << " built using generic construct."
            << G4endl << "Not applicable to the generic construct !";
    G4Exception("G4Polyhedra::Reset()", "GeomSolids1001",
                JustWarning, message, "Parameters NOT resetted.");
    return 1;
  }

  //
  // Clear old setup
  //
  G4VCSGfaceted::DeleteStuff();
  delete [] corners;
  delete enclosingCylinder;

  //
  // Rebuild polyhedra
  //
  G4ReduciblePolygon *rz =
    new G4ReduciblePolygon( original_parameters->Rmin,
                            original_parameters->Rmax,
                            original_parameters->Z_values,
                            original_parameters->Num_z_planes );
  Create( original_parameters->Start_angle,
          original_parameters->Opening_angle,
          original_parameters->numSide, rz );
  delete rz;

  return 0;
}


//
// Inside
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
EInside G4Polyhedra::Inside( const G4ThreeVector &p ) const
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
G4double G4Polyhedra::DistanceToIn( const G4ThreeVector &p,
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
G4double G4Polyhedra::DistanceToIn( const G4ThreeVector &p ) const
{
  return G4VCSGfaceted::DistanceToIn(p);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4Polyhedra::BoundingLimits(G4ThreeVector& pMin,
                                 G4ThreeVector& pMax) const
{
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
    G4Exception("G4Polyhedra::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4Polyhedra::CalculateExtent(const EAxis pAxis,
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
    G4Exception("G4Polyhedra::CalculateExtent()",
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

//
// ComputeDimensions
//
void G4Polyhedra::ComputeDimensions(       G4VPVParameterisation* p,
                                     const G4int n,
                                     const G4VPhysicalVolume* pRep )
{
  p->ComputeDimensions(*this,n,pRep);
}


//
// GetEntityType
//
G4GeometryType G4Polyhedra::GetEntityType() const
{
  return G4String("G4Polyhedra");
}


//
// Make a clone of the object
//
G4VSolid* G4Polyhedra::Clone() const
{
  return new G4Polyhedra(*this);
}


//
// Stream object contents to an output stream
//
std::ostream& G4Polyhedra::StreamInfo( std::ostream& os ) const
{
  G4int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4Polyhedra\n"
     << " Parameters: \n"
     << "    starting phi angle : " << startPhi/degree << " degrees \n"
     << "    ending phi angle   : " << endPhi/degree << " degrees \n"
     << "    number of sides    : " << numSide << " \n";
  G4int i=0;
  if (!genericPgon)
  {
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
// GetPointOnPlane
//
// Auxiliary method for get point on surface
//
G4ThreeVector
G4Polyhedra::GetPointOnPlane(G4ThreeVector p0, G4ThreeVector p1, 
                             G4ThreeVector p2, G4ThreeVector p3) const
{
  G4double lambda1, lambda2, chose,aOne,aTwo;
  G4ThreeVector t, u, v, w, Area, normal;
  aOne = 1.;
  aTwo = 1.;

  t = p1 - p0;
  u = p2 - p1;
  v = p3 - p2;
  w = p0 - p3;

  chose = G4RandFlat::shoot(0.,aOne+aTwo);
  if( (chose>=0.) && (chose < aOne) )
  {
    lambda1 = G4RandFlat::shoot(0.,1.);
    lambda2 = G4RandFlat::shoot(0.,lambda1);
    return (p2+lambda1*v+lambda2*w);    
  }

  lambda1 = G4RandFlat::shoot(0.,1.);
  lambda2 = G4RandFlat::shoot(0.,lambda1);
  return (p0+lambda1*t+lambda2*u);
}


//
// GetPointOnTriangle
//
// Auxiliary method for get point on surface
//
G4ThreeVector G4Polyhedra::GetPointOnTriangle(G4ThreeVector p1,
                                              G4ThreeVector p2,
                                              G4ThreeVector p3) const
{
  G4double lambda1,lambda2;
  G4ThreeVector v=p3-p1, w=p1-p2;

  lambda1 = G4RandFlat::shoot(0.,1.);
  lambda2 = G4RandFlat::shoot(0.,lambda1);

  return (p2 + lambda1*w + lambda2*v);
}


//
// GetPointOnSurface
//
G4ThreeVector G4Polyhedra::GetPointOnSurface() const
{
  if( !genericPgon )  // Polyhedra by faces
  {
    G4int j, numPlanes = original_parameters->Num_z_planes, Flag=0;
    G4double chose, totArea=0., Achose1, Achose2,
             rad1, rad2, sinphi1, sinphi2, cosphi1, cosphi2; 
    G4double a, b, l2, rang, totalPhi, ksi,
             area, aTop=0., aBottom=0., zVal=0.;

    G4ThreeVector p0, p1, p2, p3;
    std::vector<G4double> aVector1;
    std::vector<G4double> aVector2;
    std::vector<G4double> aVector3;

    totalPhi= (phiIsOpen) ? (endPhi-startPhi) : twopi;
    ksi = totalPhi/numSide;
    G4double cosksi = std::cos(ksi/2.);

    // Below we generate the areas relevant to our solid
    //
    for(j=0; j<numPlanes-1; j++)
    {
      a = original_parameters->Rmax[j+1];
      b = original_parameters->Rmax[j];
      l2 = sqr(original_parameters->Z_values[j]
              -original_parameters->Z_values[j+1]) + sqr(b-a);
      area = std::sqrt(l2-sqr((a-b)*cosksi))*(a+b)*cosksi;
      aVector1.push_back(area);
    }

    for(j=0; j<numPlanes-1; j++)
    {
      a = original_parameters->Rmin[j+1];//*cosksi;
      b = original_parameters->Rmin[j];//*cosksi;
      l2 = sqr(original_parameters->Z_values[j]
              -original_parameters->Z_values[j+1]) + sqr(b-a);
      area = std::sqrt(l2-sqr((a-b)*cosksi))*(a+b)*cosksi;
      aVector2.push_back(area);
    }
  
    for(j=0; j<numPlanes-1; j++)
    {
      if(phiIsOpen == true)
      {
        aVector3.push_back(0.5*(original_parameters->Rmax[j]
                               -original_parameters->Rmin[j]
                               +original_parameters->Rmax[j+1]
                               -original_parameters->Rmin[j+1])
        *std::fabs(original_parameters->Z_values[j+1]
                  -original_parameters->Z_values[j]));
      }
      else { aVector3.push_back(0.); } 
    }
  
    for(j=0; j<numPlanes-1; j++)
    {
      totArea += numSide*(aVector1[j]+aVector2[j])+2.*aVector3[j];
    }
  
    // Must include top and bottom areas
    //
    if(original_parameters->Rmax[numPlanes-1] != 0.)
    {
      a = original_parameters->Rmax[numPlanes-1];
      b = original_parameters->Rmin[numPlanes-1];
      l2 = sqr(a-b);
      aTop = std::sqrt(l2-sqr((a-b)*cosksi))*(a+b)*cosksi; 
    }

    if(original_parameters->Rmax[0] != 0.)
    {
      a = original_parameters->Rmax[0];
      b = original_parameters->Rmin[0];
      l2 = sqr(a-b);
      aBottom = std::sqrt(l2-sqr((a-b)*cosksi))*(a+b)*cosksi; 
    }

    Achose1 = 0.;
    Achose2 = numSide*(aVector1[0]+aVector2[0])+2.*aVector3[0];

    chose = G4RandFlat::shoot(0.,totArea+aTop+aBottom);
    if( (chose >= 0.) && (chose < aTop + aBottom) )
    {
      chose = G4RandFlat::shoot(startPhi,startPhi+totalPhi);
      rang = std::floor((chose-startPhi)/ksi-0.01);
      if(rang<0) { rang=0; }
      rang = std::fabs(rang);  
      sinphi1 = std::sin(startPhi+rang*ksi);
      sinphi2 = std::sin(startPhi+(rang+1)*ksi);
      cosphi1 = std::cos(startPhi+rang*ksi);
      cosphi2 = std::cos(startPhi+(rang+1)*ksi);
      chose = G4RandFlat::shoot(0., aTop + aBottom);
      if(chose>=0. && chose<aTop)
      {
        rad1 = original_parameters->Rmin[numPlanes-1];
        rad2 = original_parameters->Rmax[numPlanes-1];
        zVal = original_parameters->Z_values[numPlanes-1]; 
      }
      else 
      {
        rad1 = original_parameters->Rmin[0];
        rad2 = original_parameters->Rmax[0];
        zVal = original_parameters->Z_values[0]; 
      }
      p0 = G4ThreeVector(rad1*cosphi1,rad1*sinphi1,zVal);
      p1 = G4ThreeVector(rad2*cosphi1,rad2*sinphi1,zVal);
      p2 = G4ThreeVector(rad2*cosphi2,rad2*sinphi2,zVal);
      p3 = G4ThreeVector(rad1*cosphi2,rad1*sinphi2,zVal);
      return GetPointOnPlane(p0,p1,p2,p3); 
    }
    else
    {
      for (j=0; j<numPlanes-1; j++)
      {
        if( ((chose >= Achose1) && (chose < Achose2)) || (j == numPlanes-2) )
        { 
          Flag = j; break; 
        }
        Achose1 += numSide*(aVector1[j]+aVector2[j])+2.*aVector3[j];
        Achose2 = Achose1 + numSide*(aVector1[j+1]+aVector2[j+1])
                          + 2.*aVector3[j+1];
      }
    }

    // At this point we have chosen a subsection
    // between to adjacent plane cuts...

    j = Flag; 
    
    totArea = numSide*(aVector1[j]+aVector2[j])+2.*aVector3[j];
    chose = G4RandFlat::shoot(0.,totArea);
  
    if( (chose>=0.) && (chose<numSide*aVector1[j]) )
    {
      chose = G4RandFlat::shoot(startPhi,startPhi+totalPhi);
      rang = std::floor((chose-startPhi)/ksi-0.01);
      if(rang<0) { rang=0; }
      rang = std::fabs(rang);
      rad1 = original_parameters->Rmax[j];
      rad2 = original_parameters->Rmax[j+1];
      sinphi1 = std::sin(startPhi+rang*ksi);
      sinphi2 = std::sin(startPhi+(rang+1)*ksi);
      cosphi1 = std::cos(startPhi+rang*ksi);
      cosphi2 = std::cos(startPhi+(rang+1)*ksi);
      zVal = original_parameters->Z_values[j];
    
      p0 = G4ThreeVector(rad1*cosphi1,rad1*sinphi1,zVal);
      p1 = G4ThreeVector(rad1*cosphi2,rad1*sinphi2,zVal);

      zVal = original_parameters->Z_values[j+1];

      p2 = G4ThreeVector(rad2*cosphi2,rad2*sinphi2,zVal);
      p3 = G4ThreeVector(rad2*cosphi1,rad2*sinphi1,zVal);
      return GetPointOnPlane(p0,p1,p2,p3);
    }
    else if ( (chose >= numSide*aVector1[j])
           && (chose <= numSide*(aVector1[j]+aVector2[j])) )
    {
      chose = G4RandFlat::shoot(startPhi,startPhi+totalPhi);
      rang = std::floor((chose-startPhi)/ksi-0.01);
      if(rang<0) { rang=0; }
      rang = std::fabs(rang);
      rad1 = original_parameters->Rmin[j];
      rad2 = original_parameters->Rmin[j+1];
      sinphi1 = std::sin(startPhi+rang*ksi);
      sinphi2 = std::sin(startPhi+(rang+1)*ksi);
      cosphi1 = std::cos(startPhi+rang*ksi);
      cosphi2 = std::cos(startPhi+(rang+1)*ksi);
      zVal = original_parameters->Z_values[j];

      p0 = G4ThreeVector(rad1*cosphi1,rad1*sinphi1,zVal);
      p1 = G4ThreeVector(rad1*cosphi2,rad1*sinphi2,zVal);

      zVal = original_parameters->Z_values[j+1];
    
      p2 = G4ThreeVector(rad2*cosphi2,rad2*sinphi2,zVal);
      p3 = G4ThreeVector(rad2*cosphi1,rad2*sinphi1,zVal);
      return GetPointOnPlane(p0,p1,p2,p3);
    }

    chose = G4RandFlat::shoot(0.,2.2);
    if( (chose>=0.) && (chose < 1.) )
    {
      rang = startPhi;
    }
    else
    {
      rang = endPhi;
    } 

    cosphi1 = std::cos(rang); rad1 = original_parameters->Rmin[j];
    sinphi1 = std::sin(rang); rad2 = original_parameters->Rmax[j];

    p0 = G4ThreeVector(rad1*cosphi1,rad1*sinphi1,
                       original_parameters->Z_values[j]);
    p1 = G4ThreeVector(rad2*cosphi1,rad2*sinphi1,
                       original_parameters->Z_values[j]);

    rad1 = original_parameters->Rmax[j+1];
    rad2 = original_parameters->Rmin[j+1];

    p2 = G4ThreeVector(rad1*cosphi1,rad1*sinphi1,
                       original_parameters->Z_values[j+1]);
    p3 = G4ThreeVector(rad2*cosphi1,rad2*sinphi1,
                       original_parameters->Z_values[j+1]);
    return GetPointOnPlane(p0,p1,p2,p3);
  }
  else  // Generic polyhedra
  {
    return GetPointOnSurfaceGeneric(); 
  }
}

//
// CreatePolyhedron
//
G4Polyhedron* G4Polyhedra::CreatePolyhedron() const
{ 
  if (!genericPgon)
  {
    return new G4PolyhedronPgon( original_parameters->Start_angle,
                                 original_parameters->Opening_angle,
                                 original_parameters->numSide,
                                 original_parameters->Num_z_planes,
                                 original_parameters->Z_values,
                                 original_parameters->Rmin,
                                 original_parameters->Rmax);
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
    if (phiIsOpen)
    {
      // Triangulate open ends.  Simple ear-chopping algorithm...
      // I'm not sure how robust this algorithm is (J.Allison).
      //
      std::vector<G4bool> chopped(numCorner, false);
      std::vector<G4int*> triQuads;
      G4int remaining = numCorner;
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
            if (++iStepper >= numCorner) iStepper = 0;
          }
          while (chopped[iStepper]);
        }
        while (C < 0 && iStepper != iStarter);

        // Check triangle at B is pointing outward (an "ear").
        // Sign of z cross product determines...

        G4double BAr = corners[A].r - corners[B].r;
        G4double BAz = corners[A].z - corners[B].z;
        G4double BCr = corners[C].r - corners[B].r;
        G4double BCz = corners[C].z - corners[B].z;
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
            if (++iStarter >= numCorner) { iStarter = 0; }
          }
          while (chopped[iStarter]);
        }
      }

      // Transfer to faces...

      nNodes = (numSide + 1) * numCorner;
      nFaces = numSide * numCorner + 2 * triQuads.size();
      faces_vec = new int4[nFaces];
      G4int iface = 0;
      G4int addition = numCorner * numSide;
      G4int d = numCorner - 1;
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
      const G4double dPhi = (endPhi - startPhi) / numSide;
      G4double phi = startPhi;
      G4int ixyz = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < numCorner; ++iCorner)
        {
          xyz[ixyz][0] = corners[iCorner].r * std::cos(phi);
          xyz[ixyz][1] = corners[iCorner].r * std::sin(phi);
          xyz[ixyz][2] = corners[iCorner].z;
          if (iCorner < numCorner - 1)
          {
            faces_vec[iface][0] = ixyz + 1;
            faces_vec[iface][1] = ixyz + numCorner + 1;
            faces_vec[iface][2] = ixyz + numCorner + 2;
            faces_vec[iface][3] = ixyz + 2;
          }
          else
          {
            faces_vec[iface][0] = ixyz + 1;
            faces_vec[iface][1] = ixyz + numCorner + 1;
            faces_vec[iface][2] = ixyz + 2;
            faces_vec[iface][3] = ixyz - numCorner + 2;
          }
          ++iface;
          ++ixyz;
        }
        phi += dPhi;
      }

      // Last corners...

      for (G4int iCorner = 0; iCorner < numCorner; ++iCorner)
      {
        xyz[ixyz][0] = corners[iCorner].r * std::cos(phi);
        xyz[ixyz][1] = corners[iCorner].r * std::sin(phi);
        xyz[ixyz][2] = corners[iCorner].z;
        ++ixyz;
      }
    }
    else  // !phiIsOpen - i.e., a complete 360 degrees.
    {
      nNodes = numSide * numCorner;
      nFaces = numSide * numCorner;;
      xyz = new double3[nNodes];
      faces_vec = new int4[nFaces];
      // const G4double dPhi = (endPhi - startPhi) / numSide;
      const G4double dPhi = twopi / numSide;
      // !phiIsOpen endPhi-startPhi = 360 degrees.
      G4double phi = startPhi;
      G4int ixyz = 0, iface = 0;
      for (G4int iSide = 0; iSide < numSide; ++iSide)
      {
        for (G4int iCorner = 0; iCorner < numCorner; ++iCorner)
        {
          xyz[ixyz][0] = corners[iCorner].r * std::cos(phi);
          xyz[ixyz][1] = corners[iCorner].r * std::sin(phi);
          xyz[ixyz][2] = corners[iCorner].z;
          if (iSide < numSide - 1)
          {
            if (iCorner < numCorner - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + numCorner + 1;
              faces_vec[iface][2] = ixyz + numCorner + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + numCorner + 1;
              faces_vec[iface][2] = ixyz + 2;
              faces_vec[iface][3] = ixyz - numCorner + 2;
            }
          }
          else   // Last side joins ends...
          {
            if (iCorner < numCorner - 1)
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz + numCorner - nFaces + 1;
              faces_vec[iface][2] = ixyz + numCorner - nFaces + 2;
              faces_vec[iface][3] = ixyz + 2;
            }
            else
            {
              faces_vec[iface][0] = ixyz + 1;
              faces_vec[iface][1] = ixyz - nFaces + numCorner + 1;
              faces_vec[iface][2] = ixyz - nFaces + 2;
              faces_vec[iface][3] = ixyz - numCorner + 2;
            }
          }
          ++ixyz;
          ++iface;
        }
        phi += dPhi;
      }
    }
    G4Polyhedron* polyhedron = new G4Polyhedron;
    G4int problem = polyhedron->createPolyhedron(nNodes,nFaces,xyz,faces_vec);
    delete [] faces_vec;
    delete [] xyz;
    if (problem)
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


void  G4Polyhedra::SetOriginalParameters(G4ReduciblePolygon *rz)
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
    if(Zright>Zleft)
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
      Rmax.push_back (corners[inextr].r);
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
                                * (corners[inextl].r - corners[icurl].r));
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
   original_parameters = new G4PolyhedraHistorical;
   original_parameters->numSide = numSide;
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
    message << "Polyhedra " << GetName() << G4endl
      << "cannot be converted to Polyhedra with (Rmin,Rmaz,Z) parameters!";
    G4Exception("G4Polyhedra::SetOriginalParameters()",
                "GeomSolids0002", JustWarning, message);
#endif
    original_parameters = new G4PolyhedraHistorical;
    original_parameters->numSide = numSide;
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
  //return isConvertible;
}

#endif
