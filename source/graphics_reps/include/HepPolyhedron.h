// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: HepPolyhedron.h,v 1.1 2000-02-22 15:31:33 johna Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Class Description:
// HepPolyhedron is an intermediate class between description of a shape
// and visualization systems. It is intended to provide some service like:
//   - polygonization of shapes with triangulization (quadrilaterization)
//     of complex polygons;
//   - calculation of normals for faces and vertices;
//   - finding result of boolean operation on polyhedra;
//
// Public constructors:
//
//   HepPolyhedronBox (dx,dy,dz)
//                                        - create polyhedron for Box;
//   HepPolyhedronTrd1 (dx1,dx2,dy,dz)
//                                        - create polyhedron for G3 Trd1;
//   HepPolyhedronTrd2 (dx1,dx2,dy1,dy2,dz)
//                                        - create polyhedron for G3 Trd2;
//   HepPolyhedronTrap (dz,theta,phi, h1,bl1,tl1,alp1, h2,bl2,tl2,alp2)
//                                        - create polyhedron for G3 Trap;
//   HepPolyhedronPara (dx,dy,dz,alpha,theta,phi)
//                                        - create polyhedron for G3 Para;
//   HepPolyhedronTube (rmin,rmax,dz)
//                                        - create polyhedron for G3 Tube;
//   HepPolyhedronTubs (rmin,rmax,dz,phi1,dphi)
//                                        - create polyhedron for G3 Tubs;
//   HepPolyhedronCone (rmin1,rmax1,rmin2,rmax2,dz)
//                                        - create polyhedron for G3 Cone;
//   HepPolyhedronCons (rmin1,rmax1,rmin2,rmax2,dz,phi1,dphi)
//                                        - create polyhedron for G3 Cons;
//   HepPolyhedronPgon (phi,dphi,npdv,nz, z(*),rmin(*),rmax(*))
//                                        - create polyhedron for G3 Pgon;
//   HepPolyhedronPcon (phi,dphi,nz, z(*),rmin(*),rmax(*))
//                                        - create polyhedron for G3 Pcon;
//   HepPolyhedronSphere (rmin,rmax,phi,dphi,the,dthe)
//                                        - create polyhedron for Sphere;
//   HepPolyhedronTorus (rmin,rmax,rtor,phi,dphi)
//                                        - create polyhedron for Torus;
// Public functions:
//
//   GetNoVertices ()       - returns number of vertices;
//   GetNoFacets ()         - returns number of faces;
//   GetNextVertexIndex (index,edgeFlag) - get vertex indeces of the
//                            quadrilaterals in order;
//                            returns false when finished each face;
//   GetVertex (index)      - returns vertex by index;
//   GetNextVertex (vertex,edgeFlag) - get vertices with edge visibility
//                            of the quadrilaterals in order;
//                            returns false when finished each face;
//   GetNextVertex (vertex,edgeFlag,normal) - get vertices with edge
//                            visibility and normal of the quadrilaterals
//                            in order; returns false when finished each face;
//   GetNextEdgeIndeces (i1,i2,edgeFlag) - get indeces of the next edge;
//                            returns false for the last edge;
//   GetNextEdgeIndeces (i1,i2,edgeFlag,iface1,iface2) - get indeces of
//                            the next edge with indeces of the faces
//                            to which the edge belongs;
//                            returns false for the last edge;
//   GetNextEdge (p1,p2,edgeFlag) - get next edge;
//                            returns false for the last edge;
//   GetNextEdge (p1,p2,edgeFlag,iface1,iface2) - get next edge with indeces
//                            of the faces to which the edge belongs;
//                            returns false for the last edge;
//   GetFacet (index,n,nodes,edgeFlags=0,normals=0) - get face by index;
//   GetNextFacet (n,nodes,edgeFlags=0,normals=0) - get next face with normals
//                            at the nodes; returns false for the last face;
//   GetNormal (index)      - get normal of face given by index;
//   GetUnitNormal (index)  - get unit normal of face given by index;
//   GetNextNormal (normal) - get normals of each face in order;
//                            returns false when finished all faces;
//   GetNextUnitNormal (normal) - get normals of unit length of each face
//                            in order; returns false when finished all faces;
//   SetNumberOfRotationSteps (n) - set number of steps for whole circle;
//   ResetNumberOfRotationSteps() - reset number of steps for whole circle
//                            to default value;
// History:
//
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


#ifndef HEP_POLYHEDRON_HH
#define HEP_POLYHEDRON_HH

#include <CLHEP/config/CLHEP.h>
#include <CLHEP/config/iostream.h>
#include <CLHEP/config/TemplateFunctions.h>
#include <CLHEP/Units/SystemOfUnits.h>
#include <CLHEP/Geometry/Point3D.h>
#include <CLHEP/Geometry/Vector3D.h>
#include <CLHEP/Geometry/Normal3D.h>
#include <CLHEP/Geometry/Plane3D.h>
#include <CLHEP/Geometry/Transform3D.h>

#ifndef DEFAULT_NUMBER_OF_STEPS
#define DEFAULT_NUMBER_OF_STEPS 24
#endif

class G4Facet {
  friend class HepPolyhedron;
  friend std::ostream& operator<<(std::ostream&, const G4Facet &facet);

 private:
  struct G4Edge { int v,f; };
  G4Edge edge[4];

 public:
  G4Facet(int v1=0, int f1=0, int v2=0, int f2=0, 
	  int v3=0, int f3=0, int v4=0, int f4=0)
  { edge[0].v=v1; edge[0].f=f1; edge[1].v=v2; edge[1].f=f2;
    edge[2].v=v3; edge[2].f=f3; edge[3].v=v4; edge[3].f=f4; }
};

class HepPolyhedron {
  friend std::ostream& operator<<(std::ostream&, const HepPolyhedron &ph);

 private:
  static int fNumberOfRotationSteps;

 protected:
  int nvert, nface;
  HepPoint3D *pV;
  G4Facet    *pF;

  // Allocate memory for HepPolyhedron
  void AllocateMemory(int Nvert, int Nface);

  // Find neighbouring facet
  int FindNeighbour(int iFace, int iNode, int iOrder) const;

  // Find normal at node
  HepNormal3D FindNodeNormal(int iFace, int iNode) const;

  // Create HepPolyhedron for prism with quadrilateral base
  void CreatePrism();

  // Get number of steps for whole circle
  int GetNumberOfRotationSteps() { return fNumberOfRotationSteps; }

  // Generate facets by revolving an edge around Z-axis
  void RotateEdge(int k1, int k2, HepDouble r1, HepDouble r2,
                  int v1, int v2, int vEdge,
                  HepBoolean ifWholeCircle, int ns, int &kface);

  // Set side facets for the case of incomplete rotation
  void SetSideFacets(int ii[4], int vv[4],
                     int *kk, HepDouble *r,
                     HepDouble dphi, int ns, int &kface);

  // Create HepPolyhedron for body of revolution around Z-axis
  void RotateAroundZ(int nstep, HepDouble phi, HepDouble dphi,
                     int np1, int np2,
                     const HepDouble *z, HepDouble *r,
                     int nodeVis, int edgeVis);

  // For each edge set reference to neighbouring facet
  void SetReferences();

  // Invert the order on nodes in facets
  void InvertFacets();

 public:
  // Constructor
  HepPolyhedron(int Nvert=0, int Nface=0)
    : nvert(Nvert), nface(Nface),
      pV(Nvert ? new HepPoint3D[Nvert+1] : 0),
      pF(Nface ? new G4Facet[Nface+1] : 0) {}

  // Copy constructor
  HepPolyhedron(const HepPolyhedron & from);

  // Destructor
  virtual ~HepPolyhedron() { delete [] pV; delete [] pF; }

  // Assignment
  virtual HepPolyhedron & operator=(const HepPolyhedron & from);

  // Get number of vertices
  int GetNoVertices() const { return nvert; }

  // Get number of facets
  int GetNoFacets() const { return nface; }

  // Transform the polyhedron
  HepPolyhedron & Transform(const HepTransform3D & t);

  // Get next vertex index of the quadrilateral
  HepBoolean GetNextVertexIndex(int & index, int & edgeFlag) const;

  // Get vertex by index 
  HepPoint3D GetVertex(int index) const;

  // Get next vertex + edge visibility of the quadrilateral
  HepBoolean GetNextVertex(HepPoint3D & vertex, int & edgeFlag) const;

  // Get next vertex + edge visibility + normal of the quadrilateral
  HepBoolean GetNextVertex(HepPoint3D & vertex, int & edgeFlag,
			   HepNormal3D & normal) const;

  // Get indeces of the next edge with indeces of the faces
  HepBoolean GetNextEdgeIndeces(int & i1, int & i2, int & edgeFlag,
				int & iface1, int & iface2) const;

  // Get indeces of the next edge
  HepBoolean GetNextEdgeIndeces(int & i1, int & i2, int & edgeFlag) const;

  // Get next edge
  HepBoolean GetNextEdge(HepPoint3D &p1, HepPoint3D &p2, int &edgeFlag) const;

  // Get next edge
  HepBoolean GetNextEdge(HepPoint3D &p1, HepPoint3D &p2, int &edgeFlag,
			 int &iface1, int &iface2) const;

  // Get face by index
  void GetFacet(int iFace, int &n, int *iNodes,
		int *edgeFlags = 0, int *iFaces = 0) const;

  // Get face by index
  void GetFacet(int iFace, int &n, HepPoint3D *nodes,
		int *edgeFlags = 0, HepNormal3D *normals = 0) const;

  // Get next face with normals at the nodes
  HepBoolean GetNextFacet(int &n, HepPoint3D *nodes,
			  int *edgeFlags=0, HepNormal3D *normals=0) const;

  // Get normal of the face given by index
  HepNormal3D GetNormal(int iFace) const;

  // Get unit normal of the face given by index
  HepNormal3D GetUnitNormal(int iFace) const;

  // Get normal of the next face
  HepBoolean GetNextNormal(HepNormal3D &normal) const;

  // Get normal of unit length of the next face 
  HepBoolean GetNextUnitNormal(HepNormal3D &normal) const;

  // Set number of steps for whole circle
  static void SetNumberOfRotationSteps(int n);

  // Reset number of steps for whole circle to default value
  static void ResetNumberOfRotationSteps() {
    fNumberOfRotationSteps = DEFAULT_NUMBER_OF_STEPS;
  }
};

class HepPolyhedronTrd2 : public HepPolyhedron {
 public:
  HepPolyhedronTrd2(HepDouble Dx1, HepDouble Dx2,
		    HepDouble Dy1, HepDouble Dy2, HepDouble Dz);
};

class HepPolyhedronTrd1 : public HepPolyhedronTrd2 {
 public:
  HepPolyhedronTrd1(HepDouble Dx1, HepDouble Dx2,
		    HepDouble Dy, HepDouble Dz) :
    HepPolyhedronTrd2(Dx1, Dx2, Dy, Dy, Dz) {}
};

class HepPolyhedronBox : public HepPolyhedronTrd2 {
 public:
  HepPolyhedronBox(HepDouble Dx, HepDouble Dy, HepDouble Dz) :
    HepPolyhedronTrd2(Dx, Dx, Dy, Dy, Dz) {}
};

class HepPolyhedronTrap : public HepPolyhedron {
public:
  HepPolyhedronTrap(HepDouble Dz, HepDouble Theta, HepDouble Phi,
		    HepDouble Dy1,
		    HepDouble Dx1, HepDouble Dx2, HepDouble Alp1,
		    HepDouble Dy2,
		    HepDouble Dx3, HepDouble Dx4, HepDouble Alp2);
};

class HepPolyhedronPara : public HepPolyhedronTrap {
public:
  HepPolyhedronPara(HepDouble Dx, HepDouble Dy, HepDouble Dz,
		    HepDouble Alpha, HepDouble Theta, HepDouble Phi) :
    HepPolyhedronTrap(Dz, Theta, Phi, Dy, Dx, Dx, Alpha, Dy, Dx, Dx, Alpha) {}
};

class HepPolyhedronCons : public HepPolyhedron {
public:
  HepPolyhedronCons(HepDouble Rmn1, HepDouble Rmx1, 
		    HepDouble Rmn2, HepDouble Rmx2, HepDouble Dz,
		    HepDouble Phi1, HepDouble Dphi); 
};

class HepPolyhedronCone : public HepPolyhedronCons {
public:
  HepPolyhedronCone(HepDouble Rmn1, HepDouble Rmx1, 
		    HepDouble Rmn2, HepDouble Rmx2, HepDouble Dz)
    : HepPolyhedronCons(Rmn1, Rmx1, Rmn2, Rmx2, Dz, 0*deg, 360*deg) {}
};

class HepPolyhedronTubs : public HepPolyhedronCons {
public:
  HepPolyhedronTubs(HepDouble Rmin, HepDouble Rmax, HepDouble Dz, 
		    HepDouble Phi1, HepDouble Dphi)
    : HepPolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, Phi1, Dphi) {}
};

class HepPolyhedronTube : public HepPolyhedronCons {
public:
  HepPolyhedronTube (HepDouble Rmin, HepDouble Rmax, HepDouble Dz)
    : HepPolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, 0*deg, 360*deg) {}
};

class HepPolyhedronPgon : public HepPolyhedron {
public:
  HepPolyhedronPgon(HepDouble phi, HepDouble dphi, int npdv, int nz,
		    const HepDouble *z,
		    const HepDouble *rmin,
		    const HepDouble *rmax);
};

class HepPolyhedronPcon : public HepPolyhedronPgon {
public:
  HepPolyhedronPcon(HepDouble phi, HepDouble dphi, int nz,
		    const HepDouble *z,
		    const HepDouble *rmin,
		    const HepDouble *rmax)
    : HepPolyhedronPgon(phi, dphi, 0, nz, z, rmin, rmax) {}	
};

class HepPolyhedronSphere : public HepPolyhedron {
public:
  HepPolyhedronSphere(HepDouble rmin, HepDouble rmax,
		      HepDouble phi, HepDouble dphi,
		      HepDouble the, HepDouble dthe);
};

class HepPolyhedronTorus : public HepPolyhedron {
public:
  HepPolyhedronTorus(HepDouble rmin, HepDouble rmax, HepDouble rtor,
                    HepDouble phi, HepDouble dphi);
};

#endif /* HEP_POLYHEDRON_HH */
