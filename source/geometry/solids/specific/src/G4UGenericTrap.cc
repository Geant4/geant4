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
// Implementation of G4UGenericTrap wrapper class
// --------------------------------------------------------------------

#include "G4GenericTrap.hh"
#include "G4UGenericTrap.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

#include "G4Polyhedron.hh"
#include "G4PolyhedronArbitrary.hh"

using namespace CLHEP;

////////////////////////////////////////////////////////////////////////
//
// Constructor (generic parameters)
//
G4UGenericTrap::G4UGenericTrap(const G4String& name, G4double halfZ,
                               const std::vector<G4TwoVector>& vertices)
  : Base_t(name), fVisSubdivisions(0)
{
  SetZHalfLength(halfZ);
  Initialise(vertices);
}


////////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UGenericTrap::G4UGenericTrap(__void__& a)
  : Base_t(a), fVisSubdivisions(0), fVertices()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Destructor
//
G4UGenericTrap::~G4UGenericTrap()
{
}


//////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UGenericTrap::G4UGenericTrap(const G4UGenericTrap &source)
  : Base_t(source), fVisSubdivisions(source.fVisSubdivisions),
    fVertices(source.fVertices)
    
{
}


//////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UGenericTrap&
G4UGenericTrap::operator=(const G4UGenericTrap &source)
{
  if (this == &source) return *this;
  
  Base_t::operator=( source );
  fVertices = source.fVertices;
  fVisSubdivisions = source.fVisSubdivisions;
  
  return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Accessors & modifiers
//
G4double G4UGenericTrap::GetZHalfLength() const
{
  return GetDZ();
}
G4int G4UGenericTrap::GetNofVertices() const
{
  return fVertices.size();
}
G4TwoVector G4UGenericTrap::GetVertex(G4int index) const
{
  return G4TwoVector(GetVerticesX()[index], GetVerticesY()[index]);
}
const std::vector<G4TwoVector>& G4UGenericTrap::GetVertices() const
{
  return fVertices;
}
G4double G4UGenericTrap::GetTwistAngle(G4int index) const
{
  return GetTwist(index);
}
G4bool G4UGenericTrap::IsTwisted() const
{
  return !IsPlanar();
}
G4int G4UGenericTrap::GetVisSubdivisions() const
{
  return fVisSubdivisions;
}

void G4UGenericTrap::SetVisSubdivisions(G4int subdiv)
{
  fVisSubdivisions = subdiv;
}

void G4UGenericTrap::SetZHalfLength(G4double halfZ)
{
  SetDZ(halfZ);
}

void G4UGenericTrap::Initialise(const std::vector<G4TwoVector>& v)
{
  G4double verticesx[8], verticesy[8];
  for (G4int i=0; i<8; ++i)
  {
    fVertices.push_back(v[i]);
    verticesx[i] = v[i].x();
    verticesy[i] = v[i].y();
  }
  Initialize(verticesx, verticesy, GetZHalfLength());
}

/////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UGenericTrap::BoundingLimits(G4ThreeVector& pMin,
                                    G4ThreeVector& pMax) const
{
  U3Vector vmin, vmax;
  Extent(vmin,vmax);
  pMin.set(vmin.x(),vmin.y(),vmin.z());
  pMax.set(vmax.x(),vmax.y(),vmax.z());

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4UGenericTrap::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UGenericTrap::CalculateExtent(const EAxis pAxis,
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

  // Set bounding envelope (benv) and calculate extent
  //
  // To build the bounding envelope with plane faces each side face of
  // the trapezoid is subdivided in triangles. Subdivision is done by
  // duplication of vertices in the bases in a way that the envelope be
  // a convex polyhedron (some faces of the envelope can be degenerate)
  //
  G4double dz = GetZHalfLength();
  G4ThreeVectorList baseA(8), baseB(8);
  for (G4int i=0; i<4; ++i)
  {
    G4TwoVector va = GetVertex(i);
    G4TwoVector vb = GetVertex(i+4);
    baseA[2*i].set(va.x(),va.y(),-dz);
    baseB[2*i].set(vb.x(),vb.y(), dz);
  }
  for (G4int i=0; i<4; ++i)
  {
    G4int k1=2*i, k2=(2*i+2)%8;
    G4double ax = (baseA[k2].x()-baseA[k1].x());
    G4double ay = (baseA[k2].y()-baseA[k1].y());
    G4double bx = (baseB[k2].x()-baseB[k1].x());
    G4double by = (baseB[k2].y()-baseB[k1].y());
    G4double znorm = ax*by - ay*bx;
    baseA[k1+1] = (znorm < 0.0) ? baseA[k2] : baseA[k1];
    baseB[k1+1] = (znorm < 0.0) ? baseB[k1] : baseB[k2];
  }

  std::vector<const G4ThreeVectorList *> polygons(2);
  polygons[0] = &baseA;
  polygons[1] = &baseB;

  G4BoundingEnvelope benv(bmin,bmax,polygons);
  exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// CreatePolyhedron()
//
G4Polyhedron* G4UGenericTrap::CreatePolyhedron() const
{
  // Approximation of Twisted Side
  // Construct extra Points, if Twisted Side
  //
  G4PolyhedronArbitrary* polyhedron;
  size_t nVertices, nFacets;
  G4double fDz = GetZHalfLength();

  G4int subdivisions=0;
  G4int i;
  if(IsTwisted())
  {
    if ( GetVisSubdivisions()!= 0 )
    {
      subdivisions=GetVisSubdivisions();
    }
    else
    {
      // Estimation of Number of Subdivisions for smooth visualisation
      //
      G4double maxTwist=0.;
      for(i=0; i<4; i++)
      {
        if(GetTwistAngle(i)>maxTwist) { maxTwist=GetTwistAngle(i); }
      }

      // Computes bounding vectors for the shape
      //
      G4double Dx,Dy;
      G4ThreeVector minVec, maxVec;
      BoundingLimits(minVec,maxVec);
      Dx = 0.5*(maxVec.x()- minVec.y());
      Dy = 0.5*(maxVec.y()- minVec.y());
      if (Dy > Dx)  { Dx=Dy; }
    
      subdivisions=8*G4int(maxTwist/(Dx*Dx*Dx)*fDz);
      if (subdivisions<4)  { subdivisions=4; }
      if (subdivisions>30) { subdivisions=30; }
    }
  }
  G4int sub4=4*subdivisions;
  nVertices = 8+subdivisions*4;
  nFacets = 6+subdivisions*4;
  G4double cf=1./(subdivisions+1);
  polyhedron = new G4PolyhedronArbitrary (nVertices, nFacets);

  // Add Vertex
  //
  for (i=0;i<4;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(GetVertex(i).x(),
                                        GetVertex(i).y(),-fDz));
  }
  for( i=0;i<subdivisions;i++)
  {
    for(G4int j=0;j<4;j++)
    {
      G4TwoVector u=GetVertex(j)+cf*(i+1)*( GetVertex(j+4)-GetVertex(j));
      polyhedron->AddVertex(G4ThreeVector(u.x(),u.y(),-fDz+cf*2*fDz*(i+1)));
    }    
  }
  for (i=4;i<8;i++)
  {
    polyhedron->AddVertex(G4ThreeVector(GetVertex(i).x(),
                                        GetVertex(i).y(),fDz));
  }

  // Add Facets
  //
  polyhedron->AddFacet(1,4,3,2);  //Z-plane
  for (i=0;i<subdivisions+1;i++)
  {
    G4int is=i*4;
    polyhedron->AddFacet(5+is,8+is,4+is,1+is);
    polyhedron->AddFacet(8+is,7+is,3+is,4+is);
    polyhedron->AddFacet(7+is,6+is,2+is,3+is);
    polyhedron->AddFacet(6+is,5+is,1+is,2+is); 
  }
  polyhedron->AddFacet(5+sub4,6+sub4,7+sub4,8+sub4);  //Z-plane

  polyhedron->SetReferences();
  polyhedron->InvertFacets();

  return (G4Polyhedron*) polyhedron;
}

#endif  // G4GEOM_USE_USOLIDS
