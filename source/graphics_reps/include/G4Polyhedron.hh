// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyhedron.hh,v 1.2 1999-05-19 08:33:41 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// G4Polyhedron is an intermediate class between G4 and visualization
// systems. It is intended to provide some service like:
//   - polygonization of the G4 shapes with triangulization
//     (quadrilaterization) of complex polygons;
//   - calculation of normals for faces and vertices;
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
// Public functions:
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
// 20.06.96 Evgeni Chernyaev <Evgueni.Tcherniaev@cern.ch> - initial version
//
// 23.07.96 John Allison
// - added GetNoVertices, GetNoFacets, GetNextVertex, GetNextNormal
//
// 30.09.96 E.Chernyaev
// - added GetNextVertexIndex, GetVertex by Yasuhide Sawada
// - added GetNextUnitNormal, GetNextEdgeIndeces, GetNextEdge
// - improvements: angles now expected in radians
//                 int -> G4int, double -> G4double  
// - G4ThreeVector replaced by either G4Point3D or G4Normal3D
//
// 15.12.96 E.Chernyaev
// - private functions G4PolyhedronAlloc, G4PolyhedronPrism renamed
//   to AllocateMemory and CreatePrism
// - added private functions GetNumberOfRotationSteps, RotateEdge,
//   RotateAroundZ, SetReferences
// - rewritten G4PolyhedronCons;
// - added G4PolyhedronPara, ...Trap, ...Pgon, ...Pcon, ...Sphere, ...Torus,
//   so full List of implemented shapes now looks like:
//   BOX, TRD1, TRD2, TRAP, TUBE, TUBS, CONE, CONS, PARA, PGON, PCON,
//   SPHERE, TORUS
//
// 01.06.97 E.Chernyaev
// - RotateAroundZ modified and SetSideFacets added to allow Rmin=Rmax
//   in bodies of revolution
//
// 24.06.97 J.Allison
// - added static private member fNumberOfRotationSteps and static public
//   functions void SetNumberOfRotationSteps (G4int n) and
//   void ResetNumberOfRotationSteps ().  Modified
//   GetNumberOfRotationSteps() appropriately.  Made all three functions
//   inline (at end of this .hh file).
//   Usage:
//    G4Polyhedron::SetNumberOfRotationSteps
//     (fpView -> GetViewParameters ().GetNoOfSides ());
//    pPolyhedron = solid.CreatePolyhedron ();
//    G4Polyhedron::ResetNumberOfRotationSteps ();


#ifndef G4POLYHEDRON_HH
#define G4POLYHEDRON_HH

#include "G4ios.hh"
#include "G4VVisPrim.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"

class G4Facet {
  friend class G4Polyhedron;
  friend ostream& operator<<(ostream&, const G4Facet &facet);

 private:
  struct G4Edge { G4int v,f; };
  G4Edge edge[4];

 public:
  G4Facet(G4int v1=0, G4int f1=0, G4int v2=0, G4int f2=0, 
	  G4int v3=0, G4int f3=0, G4int v4=0, G4int f4=0)
  { edge[0].v=v1; edge[0].f=f1; edge[1].v=v2; edge[1].f=f2;
    edge[2].v=v3; edge[2].f=f3; edge[3].v=v4; edge[3].f=f4; }
};

class G4Polyhedron: public G4VVisPrim {
  friend ostream& operator<<(ostream&, const G4Polyhedron &ph);

 private:
  static G4int fNumberOfRotationSteps;

 protected:
  G4int nvert, nface;
  G4Point3D *pV;
  G4Facet   *pF;

  // Allocate memory for G4Polyhedron
  void AllocateMemory(G4int Nvert, G4int Nface);

  // Create G4Polyhedron for prism with quadrilateral base
  void CreatePrism();

  // Get number of steps for whole circle
  G4int GetNumberOfRotationSteps();

  // Generate facets by revolving an edge around Z-axis
  void RotateEdge(G4int k1, G4int k2, G4double r1, G4double r2,
                  G4int v1, G4int v2, G4int vEdge,
                  G4bool ifWholeCircle, G4int ns, G4int &kface);

  // Set side facets for the case of incomplete rotation
  void SetSideFacets(G4int ii[4], G4int vv[4],
                     G4int *kk, G4double *r,
                     G4double dphi, G4int ns, G4int &kface);

  // Create G4Polyhedron for body of revolution around Z-axis
  void RotateAroundZ(G4int nstep, G4double phi, G4double dphi,
                     G4int np1, G4int np2,
                     const G4double *z, G4double *r,
                     G4int nodeVis, G4int edgeVis);

  // For each edge set reference to neighbouring facet
  void SetReferences();

 public:
  // Constructor
  G4Polyhedron(G4int Nvert=0, G4int Nface=0)
    : nvert(Nvert), nface(Nface),
      pV(Nvert ? new G4Point3D[Nvert+1] : 0),
      pF(Nface ? new G4Facet[Nface+1] : 0) {}

  // Copy constructor
  G4Polyhedron(const G4Polyhedron &from);

  // Destructor
  ~G4Polyhedron() { delete [] pV; delete [] pF; }

  // Assignment
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
  G4Polyhedron& operator=(const G4Polyhedron &from);

  // Get number of vertices
  G4int GetNoVertices() const { return nvert; }

  // Get number of facets
  G4int GetNoFacets() const { return nface; }

  // Get next vertex index of the quadrilateral
  G4bool GetNextVertexIndex(G4int &index, G4int &edgeFlag) const;

  // Get vertex by index 
  G4Point3D GetVertex(G4int index) const;

  // Get next vertex + edge visibility of the quadrilateral
  G4bool GetNextVertex(G4Point3D &vertex, G4int &edgeFlag) const;

  // Get next vertex + edge visibility + normal of the quadrilateral
  //G4bool GetNextVertex
  //(G4Point3D &vertex, G4int &edgeFlag, G4Normal3D &normal) const;

  // Get normal of the next face 
  G4bool GetNextNormal(G4Normal3D &normal) const;

  // Get normal of unit length of the next face 
  G4bool GetNextUnitNormal(G4Normal3D &normal) const;

  // Get indeces of the next edge
  G4bool GetNextEdgeIndeces(G4int &i1, G4int &i2, G4int &edgeFlag) const;

  // Get next edge
  G4bool GetNextEdge(G4Point3D &p1, G4Point3D &p2, G4int &edgeFlag) const;

  // Set number of steps for whole circle
  static void SetNumberOfRotationSteps(G4int n);

  // Reset number of steps for whole circle to default value
  static void ResetNumberOfRotationSteps();
};

class G4PolyhedronTrd2 : public G4Polyhedron {
 public:
  G4PolyhedronTrd2(G4double Dx1, G4double Dx2,
		   G4double Dy1, G4double Dy2, G4double Dz);
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronTrd1 : public G4PolyhedronTrd2 {
 public:
  G4PolyhedronTrd1(G4double Dx1, G4double Dx2, G4double Dy, G4double Dz) :
    G4PolyhedronTrd2(Dx1, Dx2, Dy, Dy, Dz) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronBox : public G4PolyhedronTrd2 {
 public:
  G4PolyhedronBox(G4double Dx, G4double Dy, G4double Dz) :
    G4PolyhedronTrd2(Dx, Dx, Dy, Dy, Dz) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronTrap : public G4Polyhedron {
 public:
  G4PolyhedronTrap(G4double Dz, G4double Theta, G4double Phi,
                   G4double Dy1, G4double Dx1, G4double Dx2, G4double Alp1,
                   G4double Dy2, G4double Dx3, G4double Dx4, G4double Alp2);
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronPara : public G4PolyhedronTrap {
 public:
  G4PolyhedronPara(G4double Dx, G4double Dy, G4double Dz,
                   G4double Alpha, G4double Theta, G4double Phi) :
    G4PolyhedronTrap(Dz, Theta, Phi, Dy, Dx, Dx, Alpha, Dy, Dx, Dx, Alpha) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronCons : public G4Polyhedron {
 public:
  G4PolyhedronCons(G4double Rmn1, G4double Rmx1, 
		   G4double Rmn2, G4double Rmx2, G4double Dz,
		   G4double Phi1, G4double Dphi); 
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronCone : public G4PolyhedronCons {
 public:
  G4PolyhedronCone(G4double Rmn1, G4double Rmx1, 
		   G4double Rmn2, G4double Rmx2, G4double Dz)
    : G4PolyhedronCons(Rmn1, Rmx1, Rmn2, Rmx2, Dz, 0*deg, 360*deg) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronTubs : public G4PolyhedronCons {
 public:
  G4PolyhedronTubs(G4double Rmin, G4double Rmax, G4double Dz, 
		   G4double Phi1, G4double Dphi)
    : G4PolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, Phi1, Dphi) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronTube : public G4PolyhedronCons {
 public:
  G4PolyhedronTube (G4double Rmin, G4double Rmax, G4double Dz)
    : G4PolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, 0*deg, 360*deg) {}
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronPgon : public G4Polyhedron {
 public:
  G4PolyhedronPgon(G4double phi, G4double dphi, G4int npdv, G4int nz,
                   const G4double *z,
                   const G4double *rmin,
                   const G4double *rmax);
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronPcon : public G4PolyhedronPgon {
 public:
  G4PolyhedronPcon(G4double phi, G4double dphi, G4int nz,
                   const G4double *z,
                   const G4double *rmin,
                   const G4double *rmax)
    : G4PolyhedronPgon(phi, dphi, 0, nz, z, rmin, rmax) {}	
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronSphere : public G4Polyhedron {
 public:
  G4PolyhedronSphere(G4double rmin, G4double rmax,
                     G4double phi, G4double dphi,
                     G4double the, G4double dthe);
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

class G4PolyhedronTorus : public G4Polyhedron {
 public:
  G4PolyhedronTorus(G4double rmin, G4double rmax, G4double rtor,
                    G4double phi, G4double dphi);
  virtual G4Visible& operator=(const G4Visible &from);
  virtual G4VVisPrim& operator=(const G4VVisPrim &from);
};

inline G4int G4Polyhedron::GetNumberOfRotationSteps()
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNumberOfRotationSteps      Date:    11.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised: 24.06.97 *
 *                                                                     *
 * Function: Get number of steps for whole circle                      *
 *                                                                     *
 ***********************************************************************/
{
  return fNumberOfRotationSteps;
}

inline void G4Polyhedron::ResetNumberOfRotationSteps()
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::ResetNumberOfRotationSteps    Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Reset number of steps for whole circle to default value   *
 *                                                                     *
 ***********************************************************************/
{
  fNumberOfRotationSteps = 24;
}

#endif /* G4POLYHEDRON_HH */
