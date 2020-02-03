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
// G4VCSGfaceted implementation; a virtual class of a CSG type shape
// that is built entirely out of G4VCSGface faces.
//
// Author: David C. Williams (davidw@scipp.ucsc.edu)
// --------------------------------------------------------------------

#include "G4VCSGfaceted.hh"
#include "G4VCSGface.hh"
#include "G4SolidExtentList.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "Randomize.hh"

#include "G4Polyhedron.hh"   
#include "G4VGraphicsScene.hh"
#include "G4VisExtent.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex polyhedronMutex = G4MUTEX_INITIALIZER;
}

//
// Constructor
//
G4VCSGfaceted::G4VCSGfaceted( const G4String& name )
  : G4VSolid(name),
    fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.)
{
}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4VCSGfaceted::G4VCSGfaceted( __void__& a )
  : G4VSolid(a),
    fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.)
{
}

//
// Destructor
//
G4VCSGfaceted::~G4VCSGfaceted()
{
  DeleteStuff();
  delete fpPolyhedron; fpPolyhedron = nullptr;
}


//
// Copy constructor
//
G4VCSGfaceted::G4VCSGfaceted( const G4VCSGfaceted& source )
  : G4VSolid( source )
{
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  CopyStuff( source );
}


//
// Assignment operator
//
G4VCSGfaceted& G4VCSGfaceted::operator=( const G4VCSGfaceted& source )
{
  if (&source == this) { return *this; }
  
  // Copy base class data
  //
  G4VSolid::operator=(source);

  // Copy data
  //
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  DeleteStuff();
  CopyStuff( source );
  
  return *this;
}


//
// CopyStuff (protected)
//
// Copy the contents of source
//
void G4VCSGfaceted::CopyStuff( const G4VCSGfaceted& source )
{
  numFace = source.numFace;
  if (numFace == 0) { return; }    // odd, but permissable?
  
  faces = new G4VCSGface*[numFace];
  
  G4VCSGface **face = faces,
       **sourceFace = source.faces;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    *face = (*sourceFace)->Clone();
  } while( ++sourceFace, ++face < faces+numFace );
  fCubicVolume = source.fCubicVolume;
  fSurfaceArea = source.fSurfaceArea;
  fRebuildPolyhedron = false;
  fpPolyhedron = nullptr;
}


//
// DeleteStuff (protected)
//
// Delete all allocated objects
//
void G4VCSGfaceted::DeleteStuff()
{
  if (numFace)
  {
    G4VCSGface **face = faces;
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      delete *face;
    } while( ++face < faces + numFace );

    delete [] faces;
  }
  delete fpPolyhedron; fpPolyhedron = nullptr;
}


//
// CalculateExtent
//
G4bool G4VCSGfaceted::CalculateExtent( const EAxis axis,
                                       const G4VoxelLimits& voxelLimit,
                                       const G4AffineTransform& transform,
                                             G4double& min,
                                             G4double& max ) const
{
  G4SolidExtentList  extentList( axis, voxelLimit );

  //
  // Loop over all faces, checking min/max extent as we go.
  //
  G4VCSGface **face = faces;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    (*face)->CalculateExtent( axis, voxelLimit, transform, extentList );
  } while( ++face < faces + numFace );
  
  //
  // Return min/max value
  //
  return extentList.GetExtent( min, max );
}


//
// Inside
//
// It could be a good idea to override this virtual
// member to add first a simple test (such as spherical
// test or whatnot) and to call this version only if
// the simplier test fails.
//
EInside G4VCSGfaceted::Inside( const G4ThreeVector& p ) const
{
  EInside answer=kOutside;
  G4VCSGface **face = faces;
  G4double best = kInfinity;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double distance;
    EInside result = (*face)->Inside( p, kCarTolerance/2, &distance );
    if (result == kSurface) { return kSurface; }
    if (distance < best)
    {
      best = distance;
      answer = result;
    }
  } while( ++face < faces + numFace );

  return answer;
}


//
// SurfaceNormal
//
G4ThreeVector G4VCSGfaceted::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4ThreeVector answer;
  G4VCSGface **face = faces;
  G4double best = kInfinity;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double distance;
    G4ThreeVector normal = (*face)->Normal( p, &distance );
    if (distance < best)
    {
      best = distance;
      answer = normal;
    }
  } while( ++face < faces + numFace );

  return answer;
}


//
// DistanceToIn(p,v)
//
G4double G4VCSGfaceted::DistanceToIn( const G4ThreeVector& p,
                                      const G4ThreeVector& v ) const
{
  G4double distance = kInfinity;
  G4double distFromSurface = kInfinity;
  G4VCSGface **face = faces;
  G4VCSGface *bestFace = *face;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double   faceDistance,
               faceDistFromSurface;
    G4ThreeVector   faceNormal;
    G4bool    faceAllBehind;
    if ((*face)->Intersect( p, v, false, kCarTolerance/2,
                faceDistance, faceDistFromSurface,
                faceNormal, faceAllBehind ) )
    {
      //
      // Intersecting face
      //
      if (faceDistance < distance)
      {
        distance = faceDistance;
        distFromSurface = faceDistFromSurface;
        bestFace = *face;
        if (distFromSurface <= 0) { return 0; }
      }
    }
  } while( ++face < faces + numFace );
  
  if (distance < kInfinity && distFromSurface<kCarTolerance/2)
  {
    if (bestFace->Distance(p,false) < kCarTolerance/2)  { distance = 0; }
  }

  return distance;
}


//
// DistanceToIn(p)
//
G4double G4VCSGfaceted::DistanceToIn( const G4ThreeVector& p ) const
{
  return DistanceTo( p, false );
}


//
// DistanceToOut(p,v)
//
G4double G4VCSGfaceted::DistanceToOut( const G4ThreeVector& p,
                                       const G4ThreeVector& v,
                                       const G4bool calcNorm,
                                             G4bool* validNorm,
                                             G4ThreeVector* n ) const
{
  G4bool allBehind = true;
  G4double distance = kInfinity;
  G4double distFromSurface = kInfinity;
  G4ThreeVector normal;
  
  G4VCSGface **face = faces;
  G4VCSGface *bestFace = *face;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double  faceDistance,
              faceDistFromSurface;
    G4ThreeVector  faceNormal;
    G4bool    faceAllBehind;
    if ((*face)->Intersect( p, v, true, kCarTolerance/2,
                faceDistance, faceDistFromSurface,
                faceNormal, faceAllBehind ) )
    {
      //
      // Intersecting face
      //
      if ( (distance < kInfinity) || (!faceAllBehind) )  { allBehind = false; }
      if (faceDistance < distance)
      {
        distance = faceDistance;
        distFromSurface = faceDistFromSurface;
        normal = faceNormal;
        bestFace = *face;
        if (distFromSurface <= 0.)  { break; }
      }
    }
  } while( ++face < faces + numFace );
  
  if (distance < kInfinity)
  {
    if (distFromSurface <= 0.)
    {
      distance = 0.;
    }
    else if (distFromSurface<kCarTolerance/2)
    {
      if (bestFace->Distance(p,true) < kCarTolerance/2)  { distance = 0.; }
    }

    if (calcNorm)
    {
      *validNorm = allBehind;
      *n = normal;
    }
  }
  else
  { 
    if (Inside(p) == kSurface)  { distance = 0.; }
    if (calcNorm)  { *validNorm = false; }
  }

  return distance;
}


//
// DistanceToOut(p)
//
G4double G4VCSGfaceted::DistanceToOut( const G4ThreeVector& p ) const
{
  return DistanceTo( p, true );
}


//
// DistanceTo
//
// Protected routine called by DistanceToIn and DistanceToOut
//
G4double G4VCSGfaceted::DistanceTo( const G4ThreeVector& p,
                                    const G4bool outgoing ) const
{
  G4VCSGface **face = faces;
  G4double best = kInfinity;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double distance = (*face)->Distance( p, outgoing );
    if (distance < best)  { best = distance; }
  } while( ++face < faces + numFace );

  return (best < 0.5*kCarTolerance) ? 0. : best;
}


//
// DescribeYourselfTo
//
void G4VCSGfaceted::DescribeYourselfTo( G4VGraphicsScene& scene ) const
{
   scene.AddSolid( *this );
}


//
// GetExtent
//
// Define the sides of the box into which our solid instance would fit.
//
G4VisExtent G4VCSGfaceted::GetExtent() const 
{
  static const G4ThreeVector xMax(1,0,0), xMin(-1,0,0),
                             yMax(0,1,0), yMin(0,-1,0),
                             zMax(0,0,1), zMin(0,0,-1);
  static const G4ThreeVector *axes[6] =
     { &xMin, &xMax, &yMin, &yMax, &zMin, &zMax };
  
  G4double answers[6] =
     {-kInfinity, -kInfinity, -kInfinity, -kInfinity, -kInfinity, -kInfinity};

  G4VCSGface **face = faces;
  do    // Loop checking, 13.08.2015, G.Cosmo
  {    
    const G4ThreeVector **axis = axes+5 ;
    G4double* answer = answers+5;
    do    // Loop checking, 13.08.2015, G.Cosmo
    {
      G4double testFace = (*face)->Extent( **axis );
      if (testFace > *answer)  { *answer = testFace; }
    }
    while( --axis, --answer >= answers );
    
  } while( ++face < faces + numFace );
  
    return G4VisExtent( -answers[0], answers[1], 
                        -answers[2], answers[3],
                        -answers[4], answers[5]  );
}


//
// GetEntityType
//
G4GeometryType G4VCSGfaceted::GetEntityType() const
{
  return G4String("G4CSGfaceted");
}


//
// Stream object contents to an output stream
//
std::ostream& G4VCSGfaceted::StreamInfo( std::ostream& os ) const
{
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4VCSGfaceted\n"
     << " Parameters: \n"
     << "    number of faces: " << numFace << "\n"
     << "-----------------------------------------------------------\n";

  return os;
}


//
// GetCubVolStatistics
//
G4int G4VCSGfaceted::GetCubVolStatistics() const
{
  return fStatistics;
}


//
// GetCubVolEpsilon
//
G4double G4VCSGfaceted::GetCubVolEpsilon() const
{
  return fCubVolEpsilon;
}


//
// SetCubVolStatistics
//
void G4VCSGfaceted::SetCubVolStatistics(G4int st)
{
  fCubicVolume=0.;
  fStatistics=st;
}


//
// SetCubVolEpsilon
//
void G4VCSGfaceted::SetCubVolEpsilon(G4double ep)
{
  fCubicVolume=0.;
  fCubVolEpsilon=ep;
}


//
// GetAreaStatistics
//
G4int G4VCSGfaceted::GetAreaStatistics() const
{
  return fStatistics;
}


//
// GetAreaAccuracy
//
G4double G4VCSGfaceted::GetAreaAccuracy() const
{
  return fAreaAccuracy;
}


//
// SetAreaStatistics
//
void G4VCSGfaceted::SetAreaStatistics(G4int st)
{
  fSurfaceArea=0.;
  fStatistics=st;
}


//
// SetAreaAccuracy
//
void G4VCSGfaceted::SetAreaAccuracy(G4double ep)
{
  fSurfaceArea=0.;
  fAreaAccuracy=ep;
}


//
// GetCubicVolume
//
G4double G4VCSGfaceted::GetCubicVolume()
{
  if(fCubicVolume != 0.) {;}
  else   { fCubicVolume = EstimateCubicVolume(fStatistics,fCubVolEpsilon); }
  return fCubicVolume;
}


//
// GetSurfaceArea
//
G4double G4VCSGfaceted::GetSurfaceArea()
{
  if(fSurfaceArea != 0.) {;}
  else   { fSurfaceArea = EstimateSurfaceArea(fStatistics,fAreaAccuracy); }
  return fSurfaceArea;
}


//
// GetPolyhedron
//
G4Polyhedron* G4VCSGfaceted::GetPolyhedron () const
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


//
// GetPointOnSurfaceGeneric proportional to Areas of faces
// in case of GenericPolycone or GenericPolyhedra
//
G4ThreeVector G4VCSGfaceted::GetPointOnSurfaceGeneric( ) const
{
  // Preparing variables
  //
  G4ThreeVector answer=G4ThreeVector(0.,0.,0.);
  G4VCSGface **face = faces;
  G4double area = 0.;
  G4int i;
  std::vector<G4double> areas; 

  // First step: calculate surface areas
  //
  do    // Loop checking, 13.08.2015, G.Cosmo
  {
    G4double result = (*face)->SurfaceArea( );
    areas.push_back(result);
    area=area+result;
  } while( ++face < faces + numFace );

  // Second Step: choose randomly one surface
  //
  G4VCSGface **face1 = faces;
  G4double chose = area*G4UniformRand();
  G4double Achose1, Achose2;
  Achose1=0.; Achose2=0.; 
  i=0;

  do
  {
    Achose2+=areas[i];
    if(chose>=Achose1 && chose<Achose2)
    {
      G4ThreeVector point;
      point= (*face1)->GetPointOnFace();
      return point;
    }
    ++i;
    Achose1=Achose2;
  } while( ++face1 < faces + numFace );

  return answer;
}
