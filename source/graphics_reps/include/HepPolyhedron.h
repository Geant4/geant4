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
//                                        - create polyhedron for Trd1;
//   HepPolyhedronTrd2 (dx1,dx2,dy1,dy2,dz)
//                                        - create polyhedron for Trd2;
//   HepPolyhedronTrap (dz,theta,phi, h1,bl1,tl1,alp1, h2,bl2,tl2,alp2)
//                                        - create polyhedron for Trap;
//   HepPolyhedronPara (dx,dy,dz,alpha,theta,phi)
//                                        - create polyhedron for Para;
//   HepPolyhedronTube (rmin,rmax,dz)
//                                        - create polyhedron for Tube;
//   HepPolyhedronTubs (rmin,rmax,dz,phi1,dphi)
//                                        - create polyhedron for Tubs;
//   HepPolyhedronCone (rmin1,rmax1,rmin2,rmax2,dz)
//                                        - create polyhedron for Cone;
//   HepPolyhedronCons (rmin1,rmax1,rmin2,rmax2,dz,phi1,dphi)
//                                        - create polyhedron for Cons;
//   HepPolyhedronPgon (phi,dphi,npdv,nz, z(*),rmin(*),rmax(*))
//                                        - create polyhedron for Pgon;
//   HepPolyhedronPcon (phi,dphi,nz, z(*),rmin(*),rmax(*))
//                                        - create polyhedron for Pcon;
//   HepPolyhedronSphere (rmin,rmax,phi,dphi,the,dthe)
//                                        - create polyhedron for Sphere;
//   HepPolyhedronTorus (rmin,rmax,rtor,phi,dphi)
//                                        - create polyhedron for Torus;
//   HepPolyhedronTet (p0[3],p1[3],p2[3],p3[3])
//                                        - create polyhedron for Tet;
//   HepPolyhedronEllipsoid (dx,dy,dz,zcut1,zcut2)
//                                        - create polyhedron for Ellipsoid;
//   HepPolyhedronEllipticalCone(dx,dy,z,zcut1)
//                                        - create polyhedron for Elliptical cone;
//   HepPolyhedronParaboloid (r1,r2,dz,phi,dphi)
//                                        - create polyhedron for Paraboloid;
//   HepPolyhedronHype (r1,r2,tan1,tan2,halfz)
//                                        - create polyhedron for Hype;
//   HepPolyhedronHyperbolicMirror (a,h,r)
//                                        - create polyhedron for Hyperbolic mirror;
//   HepPolyhedronTetMesh (vector<p>)
//                                        - create polyhedron for tetrahedron mesh;
//   HepPolyhedronBoxMesh (sx,sy,sz,vector<p>)
//                                        - create polyhedron for box mesh;
// Public functions:
//
//   GetNoVertices ()       - returns number of vertices;
//   GetNoFacets ()         - returns number of faces;
//   GetNextVertexIndex (index,edgeFlag) - get vertex indices of the
//                            quadrilaterals in order;
//                            returns false when finished each face;
//   GetVertex (index)      - returns vertex by index;
//   GetNextVertex (vertex,edgeFlag) - get vertices with edge visibility
//                            of the quadrilaterals in order;
//                            returns false when finished each face;
//   GetNextVertex (vertex,edgeFlag,normal) - get vertices with edge
//                            visibility and normal of the quadrilaterals
//                            in order; returns false when finished each face;
//   GetNextEdgeIndices (i1,i2,edgeFlag) - get indices of the next edge;
//                            returns false for the last edge;
//   GetNextEdgeIndices (i1,i2,edgeFlag,iface1,iface2) - get indices of
//                            the next edge with indices of the faces
//                            to which the edge belongs;
//                            returns false for the last edge;
//   GetNextEdge (p1,p2,edgeFlag) - get next edge;
//                            returns false for the last edge;
//   GetNextEdge (p1,p2,edgeFlag,iface1,iface2) - get next edge with indices
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
//   GetSurfaceArea()       - get surface area of the polyhedron;
//   GetVolume()            - get volume of the polyhedron;
//   GetNumberOfRotationSteps() - get number of steps for whole circle;
//   SetVertex(index, v)    - set vertex;
//   SetFacet(index,iv1,iv2,iv3,iv4) - set facet;
//   SetReferences()        - set references to neighbouring facets;
//   JoinCoplanarFacets(tolerance) - join coplanar facets where it is possible
//   InvertFacets()         - invert the order on nodes in facets;
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
// - added GetNextUnitNormal, GetNextEdgeIndices, GetNextEdge
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
//
// 19.03.00 E.Chernyaev
// - added boolean operations (add, subtract, intersect) on polyhedra;
//
// 25.05.01 E.Chernyaev
// - added GetSurfaceArea() and GetVolume();
//
// 05.11.02 E.Chernyaev
// - added createTwistedTrap() and createPolyhedron();
//
// 06.03.05 J.Allison
// - added IsErrorBooleanProcess
//
// 20.06.05 G.Cosmo
// - added HepPolyhedronEllipsoid
//
// 18.07.07 T.Nikitina
// - added HepPolyhedronParaboloid;
//
// 21.10.09 J.Allison
// - removed IsErrorBooleanProcess (now error is returned through argument)
//
// 22.02.20 E.Chernyaev
// - added HepPolyhedronTet, HepPolyhedronHyberbolicMirror
//
// 12.05.21 E.Chernyaev
// - added TriangulatePolygon(), RotateContourAroundZ()
// - added HepPolyhedronPgon, HepPolyhedronPcon given by rz-countour
//
// 26.03.22 E.Chernyaev
// - added HepPolyhedronTetMesh
//
// 04.04.22 E.Chernyaev
// - added JoinCoplanarFacets()
//
// 07.04.22 E.Chernyaev
// - added HepPolyhedronBoxMesh

#ifndef HEP_POLYHEDRON_HH
#define HEP_POLYHEDRON_HH

#include <vector>
#include "G4Types.hh"
#include "G4TwoVector.hh"
#include "G4ThreeVector.hh"
#include "G4Point3D.hh"
#include "G4Normal3D.hh"
#include "G4Transform3D.hh"

#ifndef DEFAULT_NUMBER_OF_STEPS
#define DEFAULT_NUMBER_OF_STEPS 24
#endif

class G4Facet {
  friend class HepPolyhedron;
  friend std::ostream& operator<<(std::ostream&, const G4Facet &facet);

 private:
  struct G4Edge { G4int v,f; };
  G4Edge edge[4];

 public:
  G4Facet(G4int v1=0, G4int f1=0, G4int v2=0, G4int f2=0,
          G4int v3=0, G4int f3=0, G4int v4=0, G4int f4=0)
  { edge[0].v=v1; edge[0].f=f1; edge[1].v=v2; edge[1].f=f2;
    edge[2].v=v3; edge[2].f=f3; edge[3].v=v4; edge[3].f=f4; }
};

class HepPolyhedron {
  friend std::ostream& operator<<(std::ostream&, const HepPolyhedron &ph);

 protected:
  static G4ThreadLocal G4int fNumberOfRotationSteps;
  G4int nvert, nface;
  G4Point3D  *pV;
  G4Facet    *pF;

  // Re-allocate memory for HepPolyhedron
  void AllocateMemory(G4int Nvert, G4int Nface);

  // Find neighbouring facet
  G4int FindNeighbour(G4int iFace, G4int iNode, G4int iOrder) const;

  // Find normal at node
  G4Normal3D FindNodeNormal(G4int iFace, G4int iNode) const;

  // Create HepPolyhedron for prism with quadrilateral base
  void CreatePrism();

  // Generate facets by revolving an edge around Z-axis
  void RotateEdge(G4int k1, G4int k2, G4double r1, G4double r2,
                  G4int v1, G4int v2, G4int vEdge,
                  G4bool ifWholeCircle, G4int ns, G4int &kface);

  // Set side facets for the case of incomplete rotation
  void SetSideFacets(G4int ii[4], G4int vv[4],
                     G4int *kk, G4double *r,
                     G4double dphi, G4int ns, G4int &kface);

  // Create HepPolyhedron for body of revolution around Z-axis
  void RotateAroundZ(G4int nstep, G4double phi, G4double dphi,
                     G4int np1, G4int np2,
                     const G4double *z, G4double *r,
                     G4int nodeVis, G4int edgeVis);

  // Create HepPolyhedron for body of revolution around Z-axis
  void RotateContourAroundZ(G4int nstep, G4double phi, G4double dphi,
                            const std::vector<G4TwoVector> &rz,
                            G4int nodeVis, G4int edgeVis);

  // Triangulate closed polygon (contour)
  G4bool TriangulatePolygon(const std::vector<G4TwoVector> &polygon,
                            std::vector<G4int> &result);

  // Helper function for TriangulatePolygon()
  G4bool CheckSnip(const std::vector<G4TwoVector> &contour,
                   G4int a, G4int b, G4int c,
                   G4int n, const G4int* V);

 public:
  // Default constructor
  HepPolyhedron() : nvert(0), nface(0), pV(nullptr), pF(nullptr) {}

  // Constructor with allocation of memory
  HepPolyhedron(G4int Nvert, G4int Nface);

  // Copy constructor
  HepPolyhedron(const HepPolyhedron & from);

  // Move constructor
  HepPolyhedron(HepPolyhedron && from);

  // Destructor
  virtual ~HepPolyhedron() { delete [] pV; delete [] pF; }

  // Assignment
  HepPolyhedron & operator=(const HepPolyhedron & from);

  // Move assignment
  HepPolyhedron & operator=(HepPolyhedron && from);

  // Get number of vertices
  G4int GetNoVertices() const { return nvert; }
  G4int GetNoVerteces() const { return nvert; }  // Old spelling.

  // Get number of facets
  G4int GetNoFacets() const { return nface; }

  // Transform the polyhedron
  HepPolyhedron & Transform(const G4Transform3D & t);

  // Get next vertex index of the quadrilateral
  G4bool GetNextVertexIndex(G4int & index, G4int & edgeFlag) const;

  // Get vertex by index
  G4Point3D GetVertex(G4int index) const;

  // Get next vertex + edge visibility of the quadrilateral
  G4bool GetNextVertex(G4Point3D & vertex, G4int & edgeFlag) const;

  // Get next vertex + edge visibility + normal of the quadrilateral
  G4bool GetNextVertex(G4Point3D & vertex, G4int & edgeFlag,
                       G4Normal3D & normal) const;

  // Get indices of the next edge with indices of the faces
  G4bool GetNextEdgeIndices(G4int & i1, G4int & i2, G4int & edgeFlag,
                            G4int & iface1, G4int & iface2) const;
  G4bool GetNextEdgeIndeces(G4int & i1, G4int & i2, G4int & edgeFlag,
                            G4int & iface1, G4int & iface2) const
  {return GetNextEdgeIndices(i1,i2,edgeFlag,iface1,iface2);}  // Old spelling

  // Get indices of the next edge
  G4bool GetNextEdgeIndices(G4int & i1, G4int & i2, G4int & edgeFlag) const;
  G4bool GetNextEdgeIndeces(G4int & i1, G4int & i2, G4int & edgeFlag) const
  {return GetNextEdgeIndices(i1,i2,edgeFlag);}  // Old spelling.

  // Get next edge
  G4bool GetNextEdge(G4Point3D &p1, G4Point3D &p2, G4int &edgeFlag) const;

  // Get next edge
  G4bool GetNextEdge(G4Point3D &p1, G4Point3D &p2, G4int &edgeFlag,
                     G4int &iface1, G4int &iface2) const;

  // Get face by index
  void GetFacet(G4int iFace, G4int &n, G4int *iNodes,
                G4int *edgeFlags = nullptr, G4int *iFaces = nullptr) const;

  // Get face by index
  void GetFacet(G4int iFace, G4int &n, G4Point3D *nodes,
                G4int *edgeFlags=nullptr, G4Normal3D *normals=nullptr) const;

  // Get next face with normals at the nodes
  G4bool GetNextFacet(G4int &n, G4Point3D *nodes, G4int *edgeFlags=nullptr,
                      G4Normal3D *normals=nullptr) const;

  // Get normal of the face given by index
  G4Normal3D GetNormal(G4int iFace) const;

  // Get unit normal of the face given by index
  G4Normal3D GetUnitNormal(G4int iFace) const;

  // Get normal of the next face
  G4bool GetNextNormal(G4Normal3D &normal) const;

  // Get normal of unit length of the next face
  G4bool GetNextUnitNormal(G4Normal3D &normal) const;

  // Boolean operations
  HepPolyhedron add(const HepPolyhedron &p) const;
  HepPolyhedron subtract(const HepPolyhedron &p) const;
  HepPolyhedron intersect(const HepPolyhedron &p) const;

  // Get area of the surface of the polyhedron
  G4double GetSurfaceArea() const;

  // Get volume of the polyhedron
  G4double GetVolume() const;

  // Get number of steps for whole circle
  static G4int GetNumberOfRotationSteps();

  // Set vertex (1 <= index <= Nvert)
  void SetVertex(G4int index, const G4Point3D& v);

  // Set facet (1 <= index <= Nface)
  void SetFacet(G4int index, G4int iv1, G4int iv2, G4int iv3, G4int iv4 = 0);

  // For each edge set reference to neighbouring facet,
  // call this after all vertices and facets have been set
  void SetReferences();

  // Join couples of triangular facets to quadrangular facets
  // where it is possible
  void JoinCoplanarFacets(G4double tolerance);

  // Invert the order on nodes in facets
  void InvertFacets();

  // Set number of steps for whole circle
  static void SetNumberOfRotationSteps(G4int n);

  // Reset number of steps for whole circle to default value
  static void ResetNumberOfRotationSteps();

  /**
   * Creates polyhedron for twisted trapezoid.
   * The trapezoid is given by two bases perpendicular to the z-axis.
   *
   * @param  Dz  half length in z
   * @param  xy1 1st base (at z = -Dz)
   * @param  xy2 2nd base (at z = +Dz)
   * @return status of the operation - is non-zero in case of problem
   */
  G4int createTwistedTrap(G4double Dz,
                        const G4double xy1[][2], const G4double xy2[][2]);

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
   * @param  faces  faces (quadrilaterals or triangles)
   * @return status of the operation - is non-zero in case of problem
   */
  G4int createPolyhedron(G4int Nnodes, G4int Nfaces,
                         const G4double xyz[][3], const G4int faces[][4]);
};

class HepPolyhedronTrd2 : public HepPolyhedron
{
 public:
  HepPolyhedronTrd2(G4double Dx1, G4double Dx2,
                    G4double Dy1, G4double Dy2, G4double Dz);
  ~HepPolyhedronTrd2() override;
};

class HepPolyhedronTrd1 : public HepPolyhedronTrd2
{
 public:
  HepPolyhedronTrd1(G4double Dx1, G4double Dx2,
                    G4double Dy, G4double Dz);
  ~HepPolyhedronTrd1() override;
};

class HepPolyhedronBox : public HepPolyhedronTrd2
{
 public:
  HepPolyhedronBox(G4double Dx, G4double Dy, G4double Dz);
  ~HepPolyhedronBox() override;
};

class HepPolyhedronTrap : public HepPolyhedron
{
 public:
  HepPolyhedronTrap(G4double Dz, G4double Theta, G4double Phi,
                    G4double Dy1,
                    G4double Dx1, G4double Dx2, G4double Alp1,
                    G4double Dy2,
                    G4double Dx3, G4double Dx4, G4double Alp2);
  ~HepPolyhedronTrap() override;
};

class HepPolyhedronPara : public HepPolyhedronTrap
{
 public:
  HepPolyhedronPara(G4double Dx, G4double Dy, G4double Dz,
                    G4double Alpha, G4double Theta, G4double Phi);
  ~HepPolyhedronPara() override;
};

class HepPolyhedronParaboloid : public HepPolyhedron
{
 public:
  HepPolyhedronParaboloid(G4double r1,
                          G4double r2,
                          G4double dz,
                          G4double Phi1,
                          G4double Dphi);
  ~HepPolyhedronParaboloid() override;
};

class HepPolyhedronHype : public HepPolyhedron
{
 public:
  HepPolyhedronHype(G4double r1,
                    G4double r2,
                    G4double tan1,
                    G4double tan2,
                    G4double halfZ);
  ~HepPolyhedronHype() override;
};

class HepPolyhedronCons : public HepPolyhedron
{
 public:
  HepPolyhedronCons(G4double Rmn1, G4double Rmx1,
                    G4double Rmn2, G4double Rmx2, G4double Dz,
                    G4double Phi1, G4double Dphi);
  ~HepPolyhedronCons() override;
};

class HepPolyhedronCone : public HepPolyhedronCons
{
 public:
  HepPolyhedronCone(G4double Rmn1, G4double Rmx1,
                    G4double Rmn2, G4double Rmx2, G4double Dz);
  ~HepPolyhedronCone() override;
};

class HepPolyhedronTubs : public HepPolyhedronCons
{
 public:
  HepPolyhedronTubs(G4double Rmin, G4double Rmax, G4double Dz,
                    G4double Phi1, G4double Dphi);
  ~HepPolyhedronTubs() override;
};

class HepPolyhedronTube : public HepPolyhedronCons
{
 public:
  HepPolyhedronTube (G4double Rmin, G4double Rmax, G4double Dz);
  ~HepPolyhedronTube() override;
};

class HepPolyhedronPgon : public HepPolyhedron
{
 public:
  HepPolyhedronPgon(G4double phi, G4double dphi, G4int npdv, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
  HepPolyhedronPgon(G4double phi, G4double dphi, G4int npdv,
                    const std::vector<G4TwoVector> &rz);
  ~HepPolyhedronPgon() override;
};

class HepPolyhedronPcon : public HepPolyhedronPgon
{
 public:
  HepPolyhedronPcon(G4double phi, G4double dphi, G4int nz,
                    const G4double *z,
                    const G4double *rmin,
                    const G4double *rmax);
  HepPolyhedronPcon(G4double phi, G4double dphi,
                    const std::vector<G4TwoVector> &rz);
  ~HepPolyhedronPcon() override;
};

class HepPolyhedronSphere : public HepPolyhedron
{
 public:
  HepPolyhedronSphere(G4double rmin, G4double rmax,
                      G4double phi, G4double dphi,
                      G4double the, G4double dthe);
  ~HepPolyhedronSphere() override;
};

class HepPolyhedronTorus : public HepPolyhedron
{
 public:
  HepPolyhedronTorus(G4double rmin, G4double rmax, G4double rtor,
                     G4double phi, G4double dphi);
  ~HepPolyhedronTorus() override;
};

class HepPolyhedronTet : public HepPolyhedron
{
 public:
  HepPolyhedronTet(const G4double p0[3],
                   const G4double p1[3],
                   const G4double p2[3],
                   const G4double p3[3]);
  ~HepPolyhedronTet() override;
};

class HepPolyhedronEllipsoid : public HepPolyhedron
{
 public:
  HepPolyhedronEllipsoid(G4double dx, G4double dy, G4double dz,
                         G4double zcut1, G4double zcut2);
  ~HepPolyhedronEllipsoid() override;
};

class HepPolyhedronEllipticalCone : public HepPolyhedron
{
 public:
  HepPolyhedronEllipticalCone(G4double dx, G4double dy, G4double z,
                              G4double zcut1);
  ~HepPolyhedronEllipticalCone() override;
};

class HepPolyhedronHyperbolicMirror : public HepPolyhedron
{
 public:
  HepPolyhedronHyperbolicMirror(G4double a, G4double h, G4double r);
  ~HepPolyhedronHyperbolicMirror() override;
};

class HepPolyhedronTetMesh : public HepPolyhedron
{
 public:
  HepPolyhedronTetMesh(const std::vector<G4ThreeVector>& tetrahedra);
  ~HepPolyhedronTetMesh() override;
};

class HepPolyhedronBoxMesh : public HepPolyhedron
{
 public:
  HepPolyhedronBoxMesh(G4double sizeX, G4double sizeY, G4double sizeZ,
                       const std::vector<G4ThreeVector>& positions);
  ~HepPolyhedronBoxMesh() override;
};

#endif /* HEP_POLYHEDRON_HH */
