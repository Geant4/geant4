// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyhedron.hh,v 1.7 2000-02-23 16:57:36 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

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
//   SetNumberOfRotationSteps(G4int n) - Set number of steps for whole circle;

// History:
// 21st February 2000  Evgeni Chernaev, John Allison
// - Re-written to inherit HepPolyhedron.

#include "globals.hh"
#include "HepPolyhedron.h"
#include "G4VVisPrim.hh"

class G4Polyhedron : public HepPolyhedron, public G4VVisPrim {
public:
  G4Polyhedron ();
  G4Polyhedron (const G4Polyhedron& from);
  G4Polyhedron (const HepPolyhedron& from);
  virtual G4Visible& operator = (const G4Visible& from) {
    return G4Visible::operator = (from);
  }
  virtual G4VVisPrim& operator = (const G4VVisPrim& from) {
    return G4VVisPrim::operator = (from);
  }
  virtual HepPolyhedron& operator = (const HepPolyhedron& from) {
    return HepPolyhedron::operator = (from);
  }
  G4Polyhedron& operator = (const G4Polyhedron& from);
};

class G4PolyhedronBox: public G4Polyhedron {
public:
  G4PolyhedronBox (G4double dx, G4double dy, G4double dz);
};

class G4PolyhedronCone: public G4Polyhedron {
public:
  G4PolyhedronCone (G4double Rmn1, G4double Rmx1, 
                    G4double Rmn2, G4double Rmx2, G4double Dz);
};

class G4PolyhedronCons: public G4Polyhedron {
public:
  G4PolyhedronCons (G4double Rmn1, G4double Rmx1, 
                    G4double Rmn2, G4double Rmx2, G4double Dz,
                    G4double Phi1, G4double Dphi);
};

class G4PolyhedronPara: public G4Polyhedron {
public:
  G4PolyhedronPara (G4double Dx, G4double Dy, G4double Dz,
                    G4double Alpha, G4double Theta, G4double Phi);
};

class G4PolyhedronPcon: public G4Polyhedron {
public:
  G4PolyhedronPcon (G4double phi, G4double dphi, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
};

class G4PolyhedronPgon: public G4Polyhedron {
public:
  G4PolyhedronPgon (G4double phi, G4double dphi, G4int npdv, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
};

class G4PolyhedronSphere: public G4Polyhedron {
public:
  G4PolyhedronSphere (G4double rmin, G4double rmax,
                      G4double phi, G4double dphi,
                      G4double the, G4double dthe);
};

class G4PolyhedronTorus: public G4Polyhedron {
public:
  G4PolyhedronTorus (G4double rmin, G4double rmax, G4double rtor,
                    G4double phi, G4double dphi);
};

class G4PolyhedronTrap: public G4Polyhedron {
public:
  G4PolyhedronTrap (G4double Dz, G4double Theta, G4double Phi,
                    G4double Dy1,
                    G4double Dx1, G4double Dx2, G4double Alp1,
                    G4double Dy2,
                    G4double Dx3, G4double Dx4, G4double Alp2);
};

class G4PolyhedronTrd1: public G4Polyhedron {
public:
  G4PolyhedronTrd1 (G4double Dx1, G4double Dx2,
                    G4double Dy, G4double Dz);
};

class G4PolyhedronTrd2: public G4Polyhedron {
public:
  G4PolyhedronTrd2 (G4double Dx1, G4double Dx2,
                    G4double Dy1, G4double Dy2, G4double Dz);
};

class G4PolyhedronTube: public G4Polyhedron {
public:
  G4PolyhedronTube (G4double Rmin, G4double Rmax, G4double Dz);
};

class G4PolyhedronTubs: public G4Polyhedron {
public:
  G4PolyhedronTubs (G4double Rmin, G4double Rmax, G4double Dz, 
                    G4double Phi1, G4double Dphi);
};

#endif /* G4POLYHEDRON_HH */
