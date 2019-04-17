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
// 
// Implementation for G4UTet wrapper class
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
                     G4ThreeVector anchor,
                     G4ThreeVector p2,
                     G4ThreeVector p3,
                     G4ThreeVector p4, G4bool* degeneracyFlag)
  : Base_t(pName, U3Vector(anchor.x(),anchor.y(),anchor.z()),
                  U3Vector(p2.x(), p2.y(), p2.z()),
                  U3Vector(p3.x(), p3.y(), p3.z()),
                  U3Vector(p4.x(), p4.y(), p4.z()))
{
  G4double fXMin=std::min(std::min(std::min(anchor.x(), p2.x()),p3.x()),p4.x());
  G4double fXMax=std::max(std::max(std::max(anchor.x(), p2.x()),p3.x()),p4.x());
  G4double fYMin=std::min(std::min(std::min(anchor.y(), p2.y()),p3.y()),p4.y());
  G4double fYMax=std::max(std::max(std::max(anchor.y(), p2.y()),p3.y()),p4.y());
  G4double fZMin=std::min(std::min(std::min(anchor.z(), p2.z()),p3.z()),p4.z());
  G4double fZMax=std::max(std::max(std::max(anchor.z(), p2.z()),p3.z()),p4.z());

  G4ThreeVector fMiddle=G4ThreeVector(fXMax+fXMin,fYMax+fYMin,fZMax+fZMin)*0.5;
  G4double fMaxSize=std::max(std::max(std::max((anchor-fMiddle).mag(),
                                               (p2-fMiddle).mag()),
                                      (p3-fMiddle).mag()),
                             (p4-fMiddle).mag());
  // fV<x><y> is vector from vertex <y> to vertex <x>
  //
  G4ThreeVector fV21=p2-anchor;
  G4ThreeVector fV31=p3-anchor;
  G4ThreeVector fV41=p4-anchor;

  // make sure this is a correctly oriented set of points for the tetrahedron
  //
  G4double signed_vol=fV21.cross(fV31).dot(fV41);
  G4bool degenerate=std::fabs(signed_vol) < 1e-9*fMaxSize*fMaxSize*fMaxSize;

  if(degeneracyFlag) *degeneracyFlag=degenerate;
  else if (degenerate)
  {
    G4Exception("G4UTet::G4UTet()", "GeomSolids0002", FatalException,
                "Degenerate tetrahedron not allowed.");
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
G4UTet::~G4UTet()
{
}

///////////////////////////////////////////////////////////////////////////////
//
// Copy constructor
//
G4UTet::G4UTet(const G4UTet& rhs)
  : Base_t(rhs)
{
}


///////////////////////////////////////////////////////////////////////////////
//
// Assignment operator
//
G4UTet& G4UTet::operator = (const G4UTet& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   return *this;
}

///////////////////////////////////////////////////////////////////////////////
//
// Accessors
//
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

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UTet::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  U3Vector vmin, vmax;
  Base_t::Extent(vmin,vmax);
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
    G4Exception("G4UTet::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }
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
  G4int index = 0;
  G4double array[12];
  Base_t::GetParametersList(index, array);

  G4Polyhedron *ph=new G4Polyhedron;
  G4double xyz[4][3];
  const G4int faces[4][4]={{1,3,2,0},{1,4,3,0},{1,2,4,0},{2,3,4,0}};
  xyz[0][0]=array[0]; xyz[0][1]=array[1]; xyz[0][2]=array[2]; // fAnchor
  xyz[1][0]=array[3]; xyz[1][1]=array[4]; xyz[1][2]=array[5]; // fP2
  xyz[2][0]=array[6]; xyz[2][1]=array[7]; xyz[2][2]=array[8]; // fP3
  xyz[3][0]=array[9]; xyz[3][1]=array[10]; xyz[3][2]=array[11]; // fP4

  ph->createPolyhedron(4,4,xyz,faces);

  return ph;
}

#endif  // G4GEOM_USE_USOLIDS
