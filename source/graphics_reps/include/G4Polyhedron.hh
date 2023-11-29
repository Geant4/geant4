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

#ifndef G4POLYHEDRON_HH
#define G4POLYHEDRON_HH

// Class Description:
// G4Polyhedron is an intermediate class between G4 and visualization
// systems. It is intended to provide some service like:
//   - polygonization of the G4 shapes with triangulization
//     (quadrilaterization) of complex polygons;
//   - calculation of normals for faces and vertices.
//
// Inherits from HepPolyhedron, to which reference should be made for
// functionality.
//
// Public constructors:
//   G4PolyhedronBox(dx,dy,dz)            - create G4Polyhedron for G4 Box;
//   G4PolyhedronTrd1(dx1,dx2,dy,dz)      - create G4Polyhedron for G4 Trd1;
//   G4PolyhedronTrd2(dx1,dx2,dy1,dy2,dz) - create G4Polyhedron for G4 Trd2;
//   G4PolyhedronTrap(dz,theta,phi,
//                    h1,bl1,tl1,alp1,
//                    h2,bl2,tl2,alp2)    - create G4Polyhedron for G4 Trap;
//   G4PolyhedronPara(dx,dy,dz,
//                    alpha,theta,phi)    - create G4Polyhedron for G4 Para;
//
//   G4PolyhedronTube(rmin,rmax,dz)       - create G4Polyhedron for G4 Tube;
//   G4PolyhedronTubs(rmin,rmax,dz,
//                    phi1,dphi)          - create G4Polyhedron for G4 Tubs;
//   G4PolyhedronCone(rmin1,rmax1,
//                    rmin2,rmax2,dz)     - create G4Polyhedron for G4 Cone;
//   G4PolyhedronCons(rmin1,rmax1,
//                    rmin2,rmax2,dz,
//                    phi1,dphi)          - create G4Polyhedron for G4 Cons;
//
//   G4PolyhedronPgon(phi,dphi,npdv,nz,
//                    z(*),rmin(*),rmax(*)) - create G4Polyhedron for G4 Pgon;
//   G4PolyhedronPcon(phi,dphi,nz,
//                    z(*),rmin(*),rmax(*)) - create G4Polyhedron for G4 Pcon;
//
//   G4PolyhedronSphere(rmin,rmax,
//                      phi,dphi,the,dthe)  - create G4Polyhedron for Sphere;
//   G4PolyhedronTorus(rmin,rmax,rtor,
//                     phi,dphi)            - create G4Polyhedron for Torus;
//   G4PolyhedronTet(p0[3],p1[3],p2[3],p3[3]) - create polyhedron for Tet;
//
//   G4PolyhedronEllipsoid(dx,dy,dz,
//                         zcut1,zcut2)     - create G4Polyhedron for Ellipsoid;
//   G4PolyhedronEllipticalCone(dx,dy,z,
//                              zcut1)      - create polyhedron for Elliptical cone;
//   G4PolyhedronParaboloid(r1,r2,dz,
//                          phi,dphi)       - create polyhedron for Paraboloid;
//   G4PolyhedronHype(r1,r2,
//                    tan1,tan2,halfz)      - create polyhedron for Hype;
//   G4PolyhedronHyperbolicMirror(a,h,r)    - create polyhedron for Hyperbolic mirror;
//
//   G4PolyhedronTetMesh(vector<p>)         - create polyhedron for tetrahedron mesh;
//
//   G4PolyhedronBoxMesh(sx,sy,sz,vector<p>) - create polyhedron for box mesh;
//
// Public functions inherited from HepPolyhedron (this list might be
// incomplete):
//   GetNoVertices()  - returns number of vertices
//   GetNoFacets()    - returns number of faces
//   GetNextVertexIndex(index, edgeFlag) - get vertex indeces of the
//                      quadrilaterals in order; returns false when
//                      finished each face;
//   GetVertex(index) - returns vertex by index;
//   GetNextVertex(vertex, edgeFlag) - get vertices with edge visibility
//                      of the quadrilaterals in order;
//                      returns false when finished each face;
//   GetNextVertex(vertex, edgeFlag, normal) - get vertices with edge
//                      visibility and normal of the quadrilaterals
//                      in order; returns false when finished each face;
//   GetNextNormal(normal) - get normals of each face in order;
//                      returns false when finished all faces;
//   GetNextUnitNormal(normal) - get normals of unit length of each face
//                      in order; returns false when finished all faces;
//   GetNextEdgeIndeces(i1, i2, edgeFlag) - get indeces of the next edge;
//                      returns false for the last edge;
//   GetNextEdge(p1, p2, edgeFlag) - get next edge;
//                      returns false for the last edge;
//   SetVertex(index, v) - set vertex;
//   SetFacet(index,iv1,iv2,iv3,iv4) - set facet;
//   SetReferences()     - set references to neighbouring facets;
//   JoinCoplanarFacets(tol) - join couples of triangles into quadrilaterals;
//   InvertFacets()      - invert the order on nodes in facets;
//   SetNumberOfRotationSteps(G4int n) - set number of steps for whole circle;
//
// History:
// 21st February 2000  Evgeni Chernaev, John Allison
// - Re-written to inherit HepPolyhedron.
//
// 11.03.05 J.Allison
// - Added fNumberOfRotationStepsAtTimeOfCreation and access method.
//   (NumberOfRotationSteps is also called number of sides per circle or
//   line segments per circle - see
//   /vis/viewer/set/lineSegmentsPerCircle.)
// 20.06.05 G.Cosmo
// - Added G4PolyhedronEllipsoid.
// 09.03.06 J.Allison
// - Added operator<<.

#include "globals.hh"
#include "HepPolyhedron.h"
#include "G4Visible.hh"

class G4Polyhedron : public HepPolyhedron, public G4Visible {
public:
  // Default constructor
  G4Polyhedron ();

  // Constructors
  G4Polyhedron (G4int Nvert, G4int Nface);
  G4Polyhedron (const HepPolyhedron& from);

  // Copy and move constructors
  G4Polyhedron (const G4Polyhedron& from) = default;
  G4Polyhedron (G4Polyhedron&& from) = default;

  // Assignment and move assignment
  G4Polyhedron & operator=(const G4Polyhedron & from) = default;
  G4Polyhedron & operator=(G4Polyhedron && from) = default;

  // Destructor
  ~G4Polyhedron () override;

  G4int GetNumberOfRotationStepsAtTimeOfCreation() const {
    return fNumberOfRotationStepsAtTimeOfCreation;
  }
private:
  G4int fNumberOfRotationStepsAtTimeOfCreation = fNumberOfRotationSteps;
};

class G4PolyhedronBox: public G4Polyhedron {
public:
  G4PolyhedronBox (G4double dx, G4double dy, G4double dz);
  ~G4PolyhedronBox () override;
};

class G4PolyhedronCone: public G4Polyhedron {
public:
  G4PolyhedronCone (G4double Rmn1, G4double Rmx1,
                    G4double Rmn2, G4double Rmx2, G4double Dz);
  ~G4PolyhedronCone () override;
};

class G4PolyhedronCons: public G4Polyhedron {
public:
  G4PolyhedronCons (G4double Rmn1, G4double Rmx1,
                    G4double Rmn2, G4double Rmx2, G4double Dz,
                    G4double Phi1, G4double Dphi);
  ~G4PolyhedronCons () override;
};

class G4PolyhedronPara: public G4Polyhedron {
public:
  G4PolyhedronPara (G4double Dx, G4double Dy, G4double Dz,
                    G4double Alpha, G4double Theta, G4double Phi);
  ~G4PolyhedronPara () override;
};

class G4PolyhedronPcon: public G4Polyhedron {
public:
  G4PolyhedronPcon (G4double phi, G4double dphi, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
  G4PolyhedronPcon (G4double phi, G4double dphi,
                    const std::vector<G4TwoVector> &rz);
  ~G4PolyhedronPcon () override;
};

class G4PolyhedronPgon: public G4Polyhedron {
public:
  G4PolyhedronPgon (G4double phi, G4double dphi, G4int npdv, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
  G4PolyhedronPgon (G4double phi, G4double dphi, G4int npdv,
                    const std::vector<G4TwoVector> &rz);

  ~G4PolyhedronPgon () override;
};

class G4PolyhedronSphere: public G4Polyhedron {
public:
  G4PolyhedronSphere (G4double rmin, G4double rmax,
                      G4double phi, G4double dphi,
                      G4double the, G4double dthe);
  ~G4PolyhedronSphere () override;
};

class G4PolyhedronTet: public G4Polyhedron {
public:
  G4PolyhedronTet (const G4double p0[3],
                   const G4double p1[3],
                   const G4double p2[3],
                   const G4double p3[3]);
  ~G4PolyhedronTet () override;
};

class G4PolyhedronTorus: public G4Polyhedron {
public:
  G4PolyhedronTorus (G4double rmin, G4double rmax, G4double rtor,
                    G4double phi, G4double dphi);
  ~G4PolyhedronTorus () override;
};

class G4PolyhedronTrap: public G4Polyhedron {
public:
  G4PolyhedronTrap (G4double Dz, G4double Theta, G4double Phi,
                    G4double Dy1,
                    G4double Dx1, G4double Dx2, G4double Alp1,
                    G4double Dy2,
                    G4double Dx3, G4double Dx4, G4double Alp2);
  ~G4PolyhedronTrap () override;
};

class G4PolyhedronTrd1: public G4Polyhedron {
public:
  G4PolyhedronTrd1 (G4double Dx1, G4double Dx2,
                    G4double Dy, G4double Dz);
  ~G4PolyhedronTrd1 () override;
};

class G4PolyhedronTrd2: public G4Polyhedron {
public:
  G4PolyhedronTrd2 (G4double Dx1, G4double Dx2,
                    G4double Dy1, G4double Dy2, G4double Dz);
  ~G4PolyhedronTrd2 () override;
};

class G4PolyhedronTube: public G4Polyhedron {
public:
  G4PolyhedronTube (G4double Rmin, G4double Rmax, G4double Dz);
  ~G4PolyhedronTube () override;
};

class G4PolyhedronTubs: public G4Polyhedron {
public:
  G4PolyhedronTubs (G4double Rmin, G4double Rmax, G4double Dz,
                    G4double Phi1, G4double Dphi);
  ~G4PolyhedronTubs () override;
};

class G4PolyhedronParaboloid: public G4Polyhedron {
 public:
  G4PolyhedronParaboloid(G4double r1, G4double r2, G4double dz,
                         G4double sPhi, G4double dPhi);
  ~G4PolyhedronParaboloid () override;
};

class G4PolyhedronHype: public G4Polyhedron {
 public:
  G4PolyhedronHype(G4double r1, G4double r2, G4double tan1,
                   G4double tan2, G4double halfZ);
  ~G4PolyhedronHype () override;
};

class G4PolyhedronEllipsoid : public G4Polyhedron {
 public:
  G4PolyhedronEllipsoid(G4double dx, G4double dy, G4double dz,
                        G4double zcut1, G4double zcut2);
  ~G4PolyhedronEllipsoid () override;
};

class G4PolyhedronEllipticalCone : public G4Polyhedron {
 public:
  G4PolyhedronEllipticalCone(G4double dx, G4double dy, G4double z,
                             G4double zcut1);
  ~G4PolyhedronEllipticalCone () override;
};

class G4PolyhedronHyperbolicMirror : public G4Polyhedron {
 public:
  G4PolyhedronHyperbolicMirror(G4double a, G4double h, G4double r);
  ~G4PolyhedronHyperbolicMirror () override;
};

class G4PolyhedronTetMesh : public G4Polyhedron {
 public:
  G4PolyhedronTetMesh(const std::vector<G4ThreeVector>& tetrahedra);
  ~G4PolyhedronTetMesh () override;
};

class G4PolyhedronBoxMesh : public G4Polyhedron {
 public:
  G4PolyhedronBoxMesh(G4double sizeX, G4double sizeY, G4double sizeZ,
                      const std::vector<G4ThreeVector>& positions);

  ~G4PolyhedronBoxMesh () override;
};

std::ostream& operator<<(std::ostream& os, const G4Polyhedron&);

#endif /* G4POLYHEDRON_HH */
