// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Polyhedron.cc,v 1.3 1999-05-19 08:33:49 stesting Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
//
// G4 Polyhedron library
//
// History:
// 23.07.96 E.Chernyaev <Evgueni.Tcherniaev@cern.ch> - initial version
//
// 30.09.96 E.Chernyaev
// - added GetNextVertexIndex, GetVertex by Yasuhide Sawada
// - added GetNextUnitNormal, GetNextEdgeIndeces, GetNextEdge
//
// 15.12.96 E.Chernyaev
// - added GetNumberOfRotationSteps, RotateEdge, RotateAroundZ, SetReferences
// - rewritten G4PolyhedronCons;
// - added G4PolyhedronPara, ...Trap, ...Pgon, ...Pcon, ...Sphere, ...Torus
//
// 01.06.97 E.Chernyaev
// - modified RotateAroundZ, added SetSideFacets
  
#include "G4Polyhedron.hh"

/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron << operator                    Date:    09.05.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Print contents of G4 polyhedron                           *
 *                                                                     *
 ***********************************************************************/
ostream& operator<<(ostream &ostr, const G4Facet &facet) {
  for (G4int k=0; k<4; k++) {
    ostr << " " << facet.edge[k].v << "/" << facet.edge[k].f;
  }
  return ostr;
}

ostream& operator<<(ostream &ostr, const G4Polyhedron &ph) {
  G4int i;

  ostr << endl;
  ostr << "Nverteces=" << ph.nvert << ", Nfacets=" << ph.nface << endl;
  for (i=1; i<=ph.nvert; i++) {
     ostr << "xyz(" << i << ")=" <<
       ph.pV[i].x() << ' ' << ph.pV[i].y() << ' ' << ph.pV[i].z() << endl;
  }
  for (i=1; i<=ph.nface; i++) {
    ostr << "face(" << i << ")=" << ph.pF[i] << endl;
  }
  return ostr;
}

G4Polyhedron::G4Polyhedron(const G4Polyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron copy constructor              Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 ***********************************************************************/
{
  if (from.nvert > 0 && from.nface > 0) {
    G4int i;
    nvert = from.nvert;
    nface = from.nface;
    pV = new G4Point3D[nvert + 1];
    pF = new G4Facet[nface + 1];
    for (i=1; i<=nvert; i++) pV[i] = from.pV[i];
    for (i=1; i<=nface; i++) pF[i] = from.pF[i];
  }else{
    nvert = 0; nface = 0; pV = 0; pF = 0;
  }
}

G4Visible & G4Polyhedron::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4Polyhedron::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTrd2::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTrd2::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTrd1::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTrd1::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronBox::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronBox::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTrap::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTrap::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronPara::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronPara::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronCons::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronCons::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronCone::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronCone::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTubs::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTubs::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTube::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTube::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronPgon::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronPgon::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronPcon::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronPcon::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronSphere::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronSphere::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Visible & G4PolyhedronTorus::operator=(const G4Visible &from) {
  return G4Visible::operator = (from);
}

G4VVisPrim & G4PolyhedronTorus::operator=(const G4VVisPrim &from) {
  return G4VVisPrim::operator = (from);
}

G4Polyhedron & G4Polyhedron::operator=(const G4Polyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron = operator                    Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 * Function: Copy contents of one GEANT4 polyhedron to another         *
 *                                                                     *
 ***********************************************************************/
{
  if (this == &from) return *this;
  G4VVisPrim::operator=(from);
  delete [] pV;
  delete [] pF;
  if (from.nvert > 0  && from.nface > 0) {
    G4int i;
    nvert = from.nvert;
    nface = from.nface;
    pV = new G4Point3D[nvert + 1];
    pF = new G4Facet[nface + 1];
    for (i=1; i<=nvert; i++) pV[i] = from.pV[i];
    for (i=1; i<=nface; i++) pF[i] = from.pF[i];
  }else{
    nvert = 0; nface = 0; pV = 0; pF = 0;
  }
  return *this;
}

void G4Polyhedron::SetNumberOfRotationSteps(G4int n)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::SetNumberOfRotationSteps      Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Set number of steps for whole circle                      *
 *                                                                     *
 ***********************************************************************/
{
  const G4int  nMin = 3;
  if (n < nMin) {
    n = nMin;
    G4cout << "G4Polyhedron::SetNumberOfRotationSteps: attempt to set the"
      "\nnumber of steps per circle < " << nMin << "; forced to " << n << endl;
  }
  fNumberOfRotationSteps = n;
}

void G4Polyhedron::AllocateMemory(G4int Nvert, G4int Nface)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::AllocateMemory                Date:    19.06.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Allocate memory for GEANT4 polyhedron                     *
 *                                                                     *
 * Input: Nvert - number of nodes                                      *
 *        Nface - number of faces                                      *
 *                                                                     *
 ***********************************************************************/
{
  nvert = Nvert;
  nface = Nface;
  pV    = new G4Point3D[nvert+1];
  pF    = new G4Facet[nface+1];
}

void G4Polyhedron::CreatePrism()
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::CreatePrism                   Date:    15.07.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Set facets for a prism                                    *
 *                                                                     *
 ***********************************************************************/
{
  enum {DUMMY, BOTTOM, LEFT, BACK, RIGHT, FRONT, TOP};

  pF[1] = G4Facet(1,LEFT,  4,BACK,  3,RIGHT,  2,FRONT);
  pF[2] = G4Facet(5,TOP,   8,BACK,  4,BOTTOM, 1,FRONT);
  pF[3] = G4Facet(8,TOP,   7,RIGHT, 3,BOTTOM, 4,LEFT);
  pF[4] = G4Facet(7,TOP,   6,FRONT, 2,BOTTOM, 3,BACK);
  pF[5] = G4Facet(6,TOP,   5,LEFT,  1,BOTTOM, 2,RIGHT);
  pF[6] = G4Facet(5,FRONT, 6,RIGHT, 7,BACK,   8,LEFT);
}

void G4Polyhedron::RotateEdge(G4int k1, G4int k2, G4double r1, G4double r2,
			      G4int v1, G4int v2, G4int vEdge,
                              G4bool ifWholeCircle, G4int ns, G4int &kface)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::RotateEdge                    Date:    05.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Create set of facets by rotation of an edge around Z-axis *
 *                                                                     *
 * Input: k1, k2 - end vertices of the edge                            *
 *        r1, r2 - radiuses of the end vertices                        *
 *        v1, v2 - visibility of edges produced by rotation of the end *
 *                 vertices                                            *
 *        vEdge  - visibility of the edge                              *
 *        ifWholeCircle - is true in case of whole circle rotation     *
 *        ns     - number of discrete steps                            *
 *        r[]    - r-coordinates                                       *
 *        kface  - current free cell in the pF array                   *
 *                                                                     *
 ***********************************************************************/
{
  if (r1 == 0. && r2 == 0) return;

  G4int i;
  G4int i1  = k1;
  G4int i2  = k2;
  G4int ii1 = ifWholeCircle ? i1 : i1+ns;
  G4int ii2 = ifWholeCircle ? i2 : i2+ns;
  G4int vv  = ifWholeCircle ? vEdge : 1;

  if (ns == 1) {
    if (r1 == 0.) {
      pF[kface++]   = G4Facet(i1,0,    v2*i2,0, (i2+1),0);
    }else if (r2 == 0.) {
      pF[kface++]   = G4Facet(i1,0,    i2,0,    v1*(i1+1),0);
    }else{
      pF[kface++]   = G4Facet(i1,0,    v2*i2,0, (i2+1),0, v1*(i1+1),0);
    }
  }else{
    if (r1 == 0.) {
      pF[kface++]   = G4Facet(vv*i1,0,    v2*i2,0, vEdge*(i2+1),0);
      for (i2++,i=1; i<ns-1; i2++,i++) {
	pF[kface++] = G4Facet(vEdge*i1,0, v2*i2,0, vEdge*(i2+1),0);
      }
      pF[kface++]   = G4Facet(vEdge*i1,0, v2*i2,0, vv*ii2,0);
    }else if (r2 == 0.) {
      pF[kface++]   = G4Facet(vv*i1,0,    vEdge*i2,0, v1*(i1+1),0);
      for (i1++,i=1; i<ns-1; i1++,i++) {
	pF[kface++] = G4Facet(vEdge*i1,0, vEdge*i2,0, v1*(i1+1),0);
      }
      pF[kface++]   = G4Facet(vEdge*i1,0, vv*i2,0,    v1*ii1,0);
    }else{
      pF[kface++]   = G4Facet(vv*i1,0,    v2*i2,0, vEdge*(i2+1),0,v1*(i1+1),0);
      for (i1++,i2++,i=1; i<ns-1; i1++,i2++,i++) {
	pF[kface++] = G4Facet(vEdge*i1,0, v2*i2,0, vEdge*(i2+1),0,v1*(i1+1),0);
      }  
      pF[kface++]   = G4Facet(vEdge*i1,0, v2*i2,0, vv*ii2,0,      v1*ii1,0);
    }
  }
}

void G4Polyhedron::SetSideFacets(G4int ii[4], G4int vv[4],
				 G4int *kk, G4double *r,
                                 G4double dphi, G4int ns, G4int &kface)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::SetSideFacets                 Date:    20.05.97 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Set side facets for the case of incomplete rotation       *
 *                                                                     *
 * Input: ii[4] - indeces of original verteces                         *
 *        vv[4] - visibility of edges                                  *
 *        kk[]  - indeces of nodes                                     *
 *        r[]   - radiuses                                             *
 *        dphi  - delta phi                                            *
 *        ns     - number of discrete steps                            *
 *        kface  - current free cell in the pF array                   *
 *                                                                     *
 ***********************************************************************/
{
  G4int k1, k2, k3, k4;
  
  if (abs((G4double)(dphi-M_PI)) < perMillion) {          // half a circle
    for (G4int i=0; i<4; i++) {
      k1 = ii[i];
      k2 = (i == 3) ? ii[0] : ii[i+1];
      if (r[k1] == 0. && r[k2] == 0.) vv[i] = -1;      
    }
  }

  if (ii[1] == ii[2]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[2]];
    k3 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[2]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += ns;
    if (r[ii[2]] != 0.) k2 += ns;
    if (r[ii[3]] != 0.) k3 += ns;
    pF[kface++] = G4Facet(vv[2]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }else if (kk[ii[0]] == kk[ii[1]]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[2]];
    k3 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[1]*k1,0, vv[2]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += ns;
    if (r[ii[2]] != 0.) k2 += ns;
    if (r[ii[3]] != 0.) k3 += ns;
    pF[kface++] = G4Facet(vv[2]*k3,0, vv[1]*k2,0, vv[3]*k1,0);
  }else if (kk[ii[2]] == kk[ii[3]]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[1]];
    k3 = kk[ii[2]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[1]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += ns;
    if (r[ii[1]] != 0.) k2 += ns;
    if (r[ii[2]] != 0.) k3 += ns;
    pF[kface++] = G4Facet(vv[1]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }else{
    k1 = kk[ii[0]];
    k2 = kk[ii[1]];
    k3 = kk[ii[2]];
    k4 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[1]*k2,0, vv[2]*k3,0, vv[3]*k4,0);
    if (r[ii[0]] != 0.) k1 += ns;
    if (r[ii[1]] != 0.) k2 += ns;
    if (r[ii[2]] != 0.) k3 += ns;
    if (r[ii[3]] != 0.) k4 += ns;
    pF[kface++] = G4Facet(vv[2]*k4,0, vv[1]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }
}

void G4Polyhedron::RotateAroundZ(G4int nstep, G4double phi, G4double dphi,
                                 G4int np1, G4int np2,
				 const G4double *z, G4double *r,
				 G4int nodeVis, G4int edgeVis)
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::RotateAroundZ                 Date:    27.11.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Create G4Polyhedron for a solid produced by rotation of   *
 *           two polylines around Z-axis                               *
 *                                                                     *
 * Input: nstep - number of discrete steps, if 0 then default          *
 *        phi   - starting phi angle                                   *
 *        dphi  - delta phi                                            *
 *        np1   - number of points in external polyline                *
 *                (must be negative in case of closed polyline)        *
 *        np2   - number of points in internal polyline (may be 1)     *
 *        z[]   - z-coordinates (+z >>> -z for both polylines)         *
 *        r[]   - r-coordinates                                        *
 *        nodeVis - how to Draw edges joing consecutive positions of   *
 *                  node during rotation                               *
 *        edgeVis - how to Draw edges                                  *
 *                                                                     *
 ***********************************************************************/
{
  static G4double wholeCircle   = 2*M_PI;
    
  //   S E T   R O T A T I O N   P A R A M E T E R S

  G4bool   ifWholeCircle = abs(dphi-wholeCircle) < perMillion ? true : false;
  G4double delPhi        = ifWholeCircle ? wholeCircle : dphi;  
  G4int    nSphi         = nstep > 0 ? nstep :
    G4int(delPhi*GetNumberOfRotationSteps()/wholeCircle+.5);
  if (nSphi == 0)  nSphi = 1;
  G4int    nVphi         = ifWholeCircle ? nSphi : nSphi+1;
  G4bool   ifClosed      = np1 > 0 ? false : true;
  
  //   C O U N T   V E R T E C E S

  G4int absNp1 = abs(np1);
  G4int absNp2 = abs(np2);
  G4int i1beg = 0;
  G4int i1end = absNp1-1;
  G4int i2beg = absNp1;
  G4int i2end = absNp1+absNp2-1; 
  G4int i, j, k;

  for(i=i1beg; i<=i2end; i++) {
    if (abs(r[i]) < perMillion) r[i] = 0.;
  }

  j = 0;                                                // external nodes
  for (i=i1beg; i<=i1end; i++) {
    j += (r[i] == 0.) ? 1 : nVphi;
  }

  G4bool ifSide1 = false;                               // internal nodes
  G4bool ifSide2 = false;

  if (r[i2beg] != r[i1beg] || z[i2beg] != z[i1beg]) {
    j += (r[i2beg] == 0.) ? 1 : nVphi;
    ifSide1 = true;
  }

  for(i=i2beg+1; i<i2end; i++) {
    j += (r[i] == 0.) ? 1 : nVphi;
  }
  
  if (r[i2end] != r[i1end] || z[i2end] != z[i1end]) {
    if (absNp2 > 1) j += (r[i2end] == 0.) ? 1 : nVphi;
    ifSide2 = true;
  }

  //   C O U N T   F A C E S

  k = ifClosed ? absNp1*nSphi : (absNp1-1)*nSphi;       // external faces

  if (absNp2 > 1) {                                     // internal faces
    for(i=i2beg; i<i2end; i++) {
      if (r[i] > 0. || r[i+1] > 0.)       k += nSphi;
    }

    if (ifClosed) {
      if (r[i2end] > 0. || r[i2beg] > 0.) k += nSphi;
    }
  }

  if (!ifClosed) {                                      // side faces
    if (ifSide1 && (r[i1beg] > 0. || r[i2beg] > 0.)) k += nSphi;
    if (ifSide2 && (r[i1end] > 0. || r[i2end] > 0.)) k += nSphi;
  }

  if (!ifWholeCircle) {                                 // phi_side faces
    k += ifClosed ? 2*absNp1 : 2*(absNp1-1);
  }

  //   A L L O C A T E   M E M O R Y

  AllocateMemory(j, k);

  //   G E N E R A T E   V E R T E C E S

  G4int *kk;
  kk = new G4int[absNp1+absNp2];

  k = 1;
  for(i=i1beg; i<=i1end; i++) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = G4Point3D(0., 0., z[i]); } else { k += nVphi; }
  }

  i = i2beg;
  if (ifSide1) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = G4Point3D(0., 0., z[i]); } else { k += nVphi; }
  }else{
    kk[i] = kk[i1beg];
  }

  for(i=i2beg+1; i<i2end; i++) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = G4Point3D(0., 0., z[i]); } else { k += nVphi; }
  }

  if (absNp2 > 1) {
    i = i2end;
    if (ifSide2) {
      kk[i] = k;
      if (r[i] == 0.) pV[k] = G4Point3D(0., 0., z[i]);
    }else{
      kk[i] = kk[i1end];
    }
  }

  G4double cosPhi, sinPhi;

  for(j=0; j<nVphi; j++) {
    cosPhi = cos(phi+j*delPhi/nSphi);
    sinPhi = sin(phi+j*delPhi/nSphi);
    for(i=i1beg; i<=i2end; i++) {
      if (r[i] != 0.) pV[kk[i]+j] = G4Point3D(r[i]*cosPhi, r[i]*sinPhi, z[i]);
    }
  }

  //   G E N E R A T E   E X T E R N A L   F A C E S

  G4int v1,v2;

  k = 1;
  v2 = ifClosed ? nodeVis : 1;
  for(i=i1beg; i<i1end; i++) {
    v1 = v2;
    if (!ifClosed && i == i1end-1) {
      v2 = 1;
    }else{
      v2 = (r[i] == r[i+1] && r[i+1] == r[i+2]) ? -1 : nodeVis;
    }
    RotateEdge(kk[i], kk[i+1], r[i], r[i+1], v1, v2,
	       edgeVis, ifWholeCircle, nSphi, k);
  }
  if (ifClosed) {
    RotateEdge(kk[i1end], kk[i1beg], r[i1end],r[i1beg], nodeVis, nodeVis,
	       edgeVis, ifWholeCircle, nSphi, k);
  }

  //   G E N E R A T E   I N T E R N A L   F A C E S

  if (absNp2 > 1) {
    v2 = ifClosed ? nodeVis : 1;
    for(i=i2beg; i<i2end; i++) {
      v1 = v2;
      if (!ifClosed && i==i2end-1) {
	v2 = 1;
      }else{
	v2 = (r[i] == r[i+1] && r[i+1] == r[i+2]) ? -1 :  nodeVis;
      }
      RotateEdge(kk[i+1], kk[i], r[i+1], r[i], v2, v1,
		 edgeVis, ifWholeCircle, nSphi, k);
    }
    if (ifClosed) {
      RotateEdge(kk[i2beg], kk[i2end], r[i2beg], r[i2end], nodeVis, nodeVis,
		 edgeVis, ifWholeCircle, nSphi, k);
    }
  }

  //   G E N E R A T E   S I D E   F A C E S

  if (!ifClosed) {
    if (ifSide1) {
      RotateEdge(kk[i2beg], kk[i1beg], r[i2beg], r[i1beg], 1, 1,
		 -1, ifWholeCircle, nSphi, k);
    }
    if (ifSide2) {
      RotateEdge(kk[i1end], kk[i2end], r[i1end], r[i2end], 1, 1,
		 -1, ifWholeCircle, nSphi, k);
    }
  }

  //   G E N E R A T E   S I D E   F A C E S  for the case of incomplete circle

  if (!ifWholeCircle) {

    G4int  ii[4], vv[4];

    if (ifClosed) {
      for (i=i1beg; i<=i1end; i++) {
	ii[0] = i;
	ii[3] = (i == i1end) ? i1beg : i+1;
	ii[1] = (absNp2 == 1) ? i2beg : ii[0]+absNp1;
	ii[2] = (absNp2 == 1) ? i2beg : ii[3]+absNp1;
	vv[0] = -1;
	vv[1] = 1;
	vv[2] = -1;
	vv[3] = 1;
	SetSideFacets(ii, vv, kk, r, dphi, nSphi, k);
      }
    }else{
      for (i=i1beg; i<i1end; i++) {
	ii[0] = i;
	ii[3] = i+1;
	ii[1] = (absNp2 == 1) ? i2beg : ii[0]+absNp1;
	ii[2] = (absNp2 == 1) ? i2beg : ii[3]+absNp1;
	vv[0] = (i == i1beg)   ? 1 : -1;
	vv[1] = 1;
	vv[2] = (i == i1end-1) ? 1 : -1;
	vv[3] = 1;
	SetSideFacets(ii, vv, kk, r, dphi, nSphi, k);
      }
    }      
  }

  delete [] kk;

  if (k-1 != nface) {
    G4cerr << "Polyhedron::RotateAroundZ: number of generated faces ("
         << k-1 << ") is not equal to the number of allocated faces ("
         << nface << ")" << endl;
  }
}

void G4Polyhedron::SetReferences()
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::SetReferences                 Date:    04.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: For each edge set reference to neighbouring facet         *
 *                                                                     *
 ***********************************************************************/
{
  if (nface <= 0) return;

  struct edgeListMember {
    edgeListMember *next;
    G4int v2;
    G4int iface;
    G4int iedge;
  } *edgeList, *freeList, **headList;

  
  //   A L L O C A T E   A N D   I N I T I A T E   L I S T S

  edgeList = new edgeListMember[2*nface];
  headList = new edgeListMember*[nvert];
  
  G4int i;
  for (i=0; i<nvert; i++) {
    headList[i] = 0;
  }
  freeList = edgeList;
  for (i=0; i<2*nface-1; i++) {
    edgeList[i].next = &edgeList[i+1];
  }
  edgeList[2*nface-1].next = 0;

  //   L O O P   A L O N G   E D G E S

  G4int iface, iedge, nedge, i1, i2, k1, k2;
  edgeListMember *prev, *cur;
  
  for(iface=1; iface<=nface; iface++) {
    nedge = (pF[iface].edge[3].v == 0) ? 3 : 4;
    for (iedge=0; iedge<nedge; iedge++) {
      i1 = iedge;
      i2 = (iedge < nedge-1) ? iedge+1 : 0;
      i1 = abs(pF[iface].edge[i1].v);
      i2 = abs(pF[iface].edge[i2].v);
      k1 = min(i1,i2);
      k2 = max(i1,i2);
      
      // check head of the List corresponding to k1
      cur = headList[k1];
      if (cur == 0) {
	headList[k1] = freeList;
	freeList = freeList->next;
        cur = headList[k1];
	cur->next = 0;
	cur->v2 = k2;
	cur->iface = iface;
	cur->iedge = iedge;
        continue;
      }	

      if (cur->v2 == k2) {
        headList[k1] = cur->next;
	cur->next = freeList;
        freeList = cur;      
        pF[iface].edge[iedge].f = cur->iface;
	pF[cur->iface].edge[cur->iedge].f = iface;
        i1 = (pF[iface].edge[iedge].v < 0) ? -1 : 1;
	i2 = (pF[cur->iface].edge[cur->iedge].v < 0) ? -1 : 1;
        if (i1 != i2) {
	  G4cerr << "Polyhedron::SetReferences: different edge visibility "
	       << iface << "/" << iedge << "/"
	       << pF[iface].edge[iedge].v << " and "
	       << cur->iface << "/" << cur->iedge << "/"
	       << pF[cur->iface].edge[cur->iedge].v << endl;
	}
	continue;
      }

      // check List itself
      for (;;) {
	prev = cur;
	cur = prev->next;
	if (cur == 0) {
	  prev->next = freeList;
	  freeList = freeList->next;
	  cur = prev->next;
	  cur->next = 0;
	  cur->v2 = k2;
	  cur->iface = iface;
	  cur->iedge = iedge;
	  break;
	}

        if (cur->v2 == k2) {
	  prev->next = cur->next;
	  cur->next = freeList;
	  freeList = cur;      
	  pF[iface].edge[iedge].f = cur->iface;
	  pF[cur->iface].edge[cur->iedge].f = iface;
	  i1 = (pF[iface].edge[iedge].v < 0) ? -1 : 1;
	  i2 = (pF[cur->iface].edge[cur->iedge].v < 0) ? -1 : 1;
	    if (i1 != i2) {
	      G4cerr << "Polyhedron::SetReferences: different edge visibility "
		   << iface << "/" << iedge << "/"
		   << pF[iface].edge[iedge].v << " and "
		   << cur->iface << "/" << cur->iedge << "/"
		   << pF[cur->iface].edge[cur->iedge].v << endl;
	    }
	  break;
	}
      }
    }
  }

  //  C H E C K   T H A T   A L L   L I S T S   A R E   E M P T Y

  for (i=0; i<nvert; i++) {
    if (headList[i] != 0) {
      G4cerr << "Polyhedron::SetReferences: List " << i << " is not empty"
	   << endl;
    }
  }

  //   F R E E   M E M O R Y

  delete [] edgeList;
  delete [] headList;
}

G4bool G4Polyhedron::GetNextVertexIndex(G4int &index, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextVertexIndex           Date:    03.09.96  *
 * Author: Yasuhide Sawada                          Revised:           *
 *                                                                     *
 * Function:                                                           *
 *                                                                     *
 ***********************************************************************/
{
  static G4int iFace = 1;
  static G4int iQVertex = 0;
  G4int vIndex = pF[iFace].edge[iQVertex].v;

//edgeFlag = vIndex > 0 ? true : false;
  edgeFlag = vIndex > 0 ? 1 : 0;
  index = abs(vIndex);

  if (iQVertex >= 3 || pF[iFace].edge[iQVertex+1].v == 0) {
    iQVertex = 0;
    if (++iFace > nface) iFace = 1;
    return false;  // Last Edge
  }else{
    ++iQVertex;
    return true;  // not Last Edge
  }
}

G4Point3D G4Polyhedron::GetVertex(G4int index) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetVertex                    Date:    03.09.96  *
 * Author: Yasuhide Sawada                          Revised:           *
 *                                                                     *
 * Function: Get vertex of the index.                                  *
 *                                                                     *
 ***********************************************************************/
{
  if (index <= 0 || index > nvert){
    G4cerr << "Error: irrelevant vertex label " << index << endl;
      return G4Point3D();
  }
  return pV[index];
}

G4bool G4Polyhedron::GetNextVertex(G4Point3D &vertex, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextVertex                Date:    22.07.96  *
 * Author: John Allison                             Revised:           *
 *                                                                     *
 * Function: Get vertices of the quadrilaterals in order for each      *
 *           face in face order.  Returns false when finished each     *
 *           face.                                                     *
 *                                                                     *
 ***********************************************************************/
{
  G4int index;
  G4bool rep = GetNextVertexIndex(index, edgeFlag);
  vertex = pV[index];
  return rep;
}

G4bool G4Polyhedron::GetNextNormal(G4Normal3D &normal) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextNormal                Date:    22.07.96  *
 * Author: John Allison                             Revised:           *
 *                                                                     *
 * Function: Get normals of each face in face order.  Returns false    *
 *           when finished all faces.                                  *
 *                                                                     *
 ***********************************************************************/
{
  static G4int iFace = 1;
  G4int vIndex;

  vIndex       = pF[iFace].edge[0].v;
  G4Point3D v0 = pV[abs(vIndex)];
  vIndex       = pF[iFace].edge[1].v;
  G4Point3D v1 = pV[abs(vIndex)];
  vIndex       = pF[iFace].edge[2].v;
  G4Point3D v2 = pV[abs(vIndex)];
  vIndex       = pF[iFace].edge[3].v;
  G4Point3D v3 = vIndex > 0 ? pV[abs(vIndex)] : v0;

  // Use diagonals - safer.
  normal = (v2 - v0).cross(v3 - v1);

  if (++iFace <= nface) {
    return true;
  }
  else {
    iFace = 1;
    return false;
  }
}

G4bool G4Polyhedron::GetNextUnitNormal(G4Normal3D &normal) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextUnitNormal            Date:    16.09.96  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get normals of unit length of each face in face order.    *
 *           Returns false when finished all faces.                    *
 *                                                                     *
 ***********************************************************************/
{
  G4bool rep = GetNextNormal(normal);
  normal = normal.unit();
  return rep;
}

G4bool
G4Polyhedron::GetNextEdgeIndeces(G4int &i1, G4int &i2, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextEdgeIndeces           Date:    30.09.96  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get indeces of the next edge.                             *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  static G4int iFace    = 1;
  static G4int iQVertex = 0;
  static G4int iOrder   = 1;
  G4int  k1, k2, kflag;

  if (iFace == 1 && iQVertex == 0) {
    k2 = pF[nface].edge[0].v;
    k1 = pF[nface].edge[3].v;
    if (k1 == 0) k1 = pF[nface].edge[2].v;
    if (abs(k1) > abs(k2)) iOrder = -1;
  }

  do {
    k1    = pF[iFace].edge[iQVertex].v;
    kflag = k1;
    k1    = abs(k1);
    if (iQVertex >= 3 || pF[iFace].edge[iQVertex+1].v == 0) {
      iQVertex = 0;
      k2 = abs(pF[iFace].edge[iQVertex].v);
      iFace++;
    }else{
      iQVertex++;
      k2 = abs(pF[iFace].edge[iQVertex].v);
    }
  } while (iOrder*k1 > iOrder*k2);

  i1 = k1; i2 = k2; edgeFlag = kflag > 0 ? 1 : -1;

  if (iFace > nface) {
    iFace  = 1; iOrder = 1;
    return false;
  }else{
    return true;
  }
}

G4bool
G4Polyhedron::GetNextEdge(G4Point3D &p1, G4Point3D &p2, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::GetNextEdge                  Date:    30.09.96  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get next edge.                                            *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  G4int i1,i2;
  G4bool rep = GetNextEdgeIndeces(i1,i2,edgeFlag);
  p1 = pV[i1];
  p2 = pV[i2];
  return rep;
}

G4PolyhedronTrd2::G4PolyhedronTrd2(G4double Dx1, G4double Dx2,
                                   G4double Dy1, G4double Dy2, G4double Dz)
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronTrd2                            Date:    22.07.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Create GEANT4 TRD2-trapezoid                              *
 *                                                                     *
 * Input: Dx1 - half-length along X at -Dz           8----7            *
 *        Dx2 - half-length along X ay +Dz        5----6  !            *
 *        Dy1 - half-length along Y ay -Dz        !  4-!--3            *
 *        Dy2 - half-length along Y ay +Dz        1----2               *
 *        Dz  - half-length along Z                                    *
 *                                                                     *
 ***********************************************************************/
{
  AllocateMemory(8,6);

  pV[1] = G4Point3D(-Dx1,-Dy1,-Dz);
  pV[2] = G4Point3D( Dx1,-Dy1,-Dz);
  pV[3] = G4Point3D( Dx1, Dy1,-Dz);
  pV[4] = G4Point3D(-Dx1, Dy1,-Dz);
  pV[5] = G4Point3D(-Dx2,-Dy2, Dz);
  pV[6] = G4Point3D( Dx2,-Dy2, Dz);
  pV[7] = G4Point3D( Dx2, Dy2, Dz);
  pV[8] = G4Point3D(-Dx2, Dy2, Dz);

  CreatePrism();
}

G4PolyhedronTrap::G4PolyhedronTrap(G4double Dz, G4double Theta, G4double Phi,
                     G4double Dy1, G4double Dx1, G4double Dx2, G4double Alp1,
                     G4double Dy2, G4double Dx3, G4double Dx4, G4double Alp2)
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronTrap                            Date:    20.11.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Create GEANT4 TRAP-trapezoid                              *
 *                                                                     *
 * Input: DZ   - half-length in Z                                      *
 *        Theta,Phi - polar angles of the line joining centres of the  *
 *                    faces at Z=-Dz and Z=+Dz                         *
 *        Dy1  - half-length in Y of the face at Z=-Dz                 *
 *        Dx1  - half-length in X of low edge of the face at Z=-Dz     *
 *        Dx2  - half-length in X of top edge of the face at Z=-Dz     *
 *        Alp1 - angle between Y-axis and the median joining top and   *
 *               low edges of the face at Z=-Dz                        *
 *        Dy2  - half-length in Y of the face at Z=+Dz                 *
 *        Dx3  - half-length in X of low edge of the face at Z=+Dz     *
 *        Dx4  - half-length in X of top edge of the face at Z=+Dz     *
 *        Alp2 - angle between Y-axis and the median joining top and   *
 *               low edges of the face at Z=+Dz                        *
 *                                                                     *
 ***********************************************************************/
{
  G4double DzTthetaCphi = Dz*tan(Theta)*cos(Phi);
  G4double DzTthetaSphi = Dz*tan(Theta)*sin(Phi);
  G4double Dy1Talp1 = Dy1*tan(Alp1);
  G4double Dy2Talp2 = Dy2*tan(Alp2);
  
  AllocateMemory(8,6);

  pV[1] = G4Point3D(-DzTthetaCphi-Dy1Talp1-Dx1,-DzTthetaSphi-Dy1,-Dz);
  pV[2] = G4Point3D(-DzTthetaCphi-Dy1Talp1+Dx1,-DzTthetaSphi-Dy1,-Dz);
  pV[3] = G4Point3D(-DzTthetaCphi+Dy1Talp1+Dx2,-DzTthetaSphi+Dy1,-Dz);
  pV[4] = G4Point3D(-DzTthetaCphi+Dy1Talp1-Dx2,-DzTthetaSphi+Dy1,-Dz);
  pV[5] = G4Point3D( DzTthetaCphi-Dy2Talp2-Dx3, DzTthetaSphi-Dy2, Dz);
  pV[6] = G4Point3D( DzTthetaCphi-Dy2Talp2+Dx3, DzTthetaSphi-Dy2, Dz);
  pV[7] = G4Point3D( DzTthetaCphi+Dy2Talp2+Dx4, DzTthetaSphi+Dy2, Dz);
  pV[8] = G4Point3D( DzTthetaCphi+Dy2Talp2-Dx4, DzTthetaSphi+Dy2, Dz);

  CreatePrism();
}

G4PolyhedronCons::G4PolyhedronCons(G4double Rmn1, G4double Rmx1,
                                   G4double Rmn2, G4double Rmx2, 
                                   G4double Dz, G4double Phi1, G4double Dphi) 
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronCons::G4PolyhedronCons          Date:    15.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised: 15.12.96 *
 *                                                                     *
 * Function: Constructor for CONS, TUBS, CONE, TUBE                    *
 *                                                                     *
 * Input: Rmn1, Rmx1 - inside and outside radiuses at -Dz              *
 *        Rmn2, Rmx2 - inside and outside radiuses at +Dz              *
 *        Dz         - half length in Z                                *
 *        Phi1       - starting angle of the segment                   *
 *        Dphi       - segment range                                   *
 *                                                                     *
 ***********************************************************************/
{
  static G4double wholeCircle=2*M_PI;

  //   C H E C K   I N P U T   P A R A M E T E R S

  G4int k = 0;
  if (Rmn1 < 0. || Rmx1 < 0. || Rmn2 < 0. || Rmx2 < 0.)        k = 1;
  if (Rmn1 > Rmx1 || Rmn2 > Rmx2)                              k = 1;
  if (Rmn1 == Rmx1 && Rmn2 == Rmx2)                            k = 1;

  if (Dz <= 0.) k += 2;
 
  G4double phi1, phi2, dphi;
  if (Dphi < 0.) {
    phi2 = Phi1; phi1 = phi2 - Dphi;
  }else if (Dphi == 0.) {
    phi1 = Phi1; phi2 = phi1 + wholeCircle;
  }else{
    phi1 = Phi1; phi2 = phi1 + Dphi;
  }
  dphi  = phi2 - phi1;
  if (abs(dphi-wholeCircle) < perMillion) dphi = wholeCircle;
  if (dphi > wholeCircle) k += 4; 

  if (k != 0) {
    G4cerr << "G4PolyhedronCone(s)/Tube(s): error in input parameters";
    if ((k & 1) != 0) G4cerr << " (radiuses)";
    if ((k & 2) != 0) G4cerr << " (half-length)";
    if ((k & 4) != 0) G4cerr << " (angles)";
    G4cerr << endl;
    G4cerr << " Rmn1=" << Rmn1 << " Rmx1=" << Rmx1;
    G4cerr << " Rmn2=" << Rmn2 << " Rmx2=" << Rmx2;
    G4cerr << " Dz=" << Dz << " Phi1=" << Phi1 << " Dphi=" << Dphi << endl;
    return;
  }
  
  //   P R E P A R E   T W O   P O L Y L I N E S

  G4double zz[4], rr[4];
  zz[0] =  Dz; 
  zz[1] = -Dz; 
  zz[2] =  Dz; 
  zz[3] = -Dz; 
  rr[0] =  Rmx2;
  rr[1] =  Rmx1;
  rr[2] =  Rmn2;
  rr[3] =  Rmn1;

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, phi1, dphi, 2, 2, zz, rr, 1, -1); 
  SetReferences();
}

G4PolyhedronPgon::G4PolyhedronPgon(G4double phi,
				   G4double dphi,
				   G4int    npdv,
				   G4int    nz,
				   const G4double *z,
				   const G4double *rmin,
				   const G4double *rmax)
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronPgon                            Date:    09.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Constructor of polyhedron for PGON, PCON                  *
 *                                                                     *
 * Input: phi  - initial phi                                           *
 *        dphi - delta phi                                             *
 *        npdv - number of steps along phi                             *
 *        nz   - number of z-planes (at least two)                     *
 *        z[]  - z coordinates of the slices                           *
 *        rmin[] - smaller r at the slices                             *
 *        rmax[] - bigger  r at the slices                             *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  if (dphi <= 0. || dphi > 2*M_PI) {
    G4cerr << "G4PolyhedronPgon/Pcon: wrong delta phi = " << dphi << endl;
    return;
  }    
    
  if (nz < 2) {
    G4cerr << "G4PolyhedronPgon/Pcon: number of z-planes less than two = " << nz
         << endl;
    return;
  }

  if (npdv < 0) {
    G4cerr << "G4PolyhedronPgon/Pcon: error in number of phi-steps =" << npdv
         << endl;
    return;
  }

  G4int i;
  for (i=0; i<nz; i++) {
    if (rmin[i] < 0. || rmax[i] < 0. || rmin[i] > rmax[i]) {
      G4cerr << "G4PolyhedronPgon: error in radiuses rmin[" << i << "]="
	   << rmin[i] << " rmax[" << i << "]=" << rmax[i] << endl;
      return;
    }
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4double *zz, *rr;
  zz = new G4double[2*nz];
  rr = new G4double[2*nz];

  if (z[0] > z[nz-1]) {
    for (i=0; i<nz; i++) {
      zz[i]    = z[i];
      rr[i]    = rmax[i];
      zz[i+nz] = z[i];
      rr[i+nz] = rmin[i];
    }
  }else{
    for (i=0; i<nz; i++) {
      zz[i]    = z[nz-i-1];
      rr[i]    = rmax[nz-i-1];
      zz[i+nz] = z[nz-i-1];
      rr[i+nz] = rmin[nz-i-1];
    }
  }

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(npdv, phi, dphi, nz, nz, zz, rr, 1, (npdv == 0) ? -1 : 1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

G4PolyhedronSphere::G4PolyhedronSphere(G4double rmin, G4double rmax,
				       G4double phi, G4double dphi,
				       G4double the, G4double dthe)
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronSphere                          Date:    11.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Constructor of polyhedron for SPHERE                      *
 *                                                                     *
 * Input: rmin - internal radius                                       *
 *        rmax - external radius                                       *
 *        phi  - initial phi                                           *
 *        dphi - delta phi                                             *
 *        the  - initial theta                                         *
 *        dthe - delta theta                                           *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  if (dphi <= 0. || dphi > 2*M_PI) {
    G4cerr << "G4PolyhedronSphere: wrong delta phi = " << dphi << endl;
    return;
  }    

  if (the < 0. || the > M_PI) {
    G4cerr << "G4PolyhedronSphere: wrong theta = " << the << endl;
    return;
  }    
  
  if (dthe <= 0. || dthe > M_PI) {
    G4cerr << "G4PolyhedronSphere: wrong delta theta = " << dthe << endl;
    return;
  }    

  if (the+dthe > M_PI) {
    G4cerr << "G4PolyhedronSphere: wrong theta + delta theta = "
         << the << " " << dthe << endl;
    return;
  }    
  
  if (rmin < 0. || rmin >= rmax) {
    G4cerr << "G4PolyhedronSphere: error in radiuses"
         << " rmin=" << rmin << " rmax=" << rmax << endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4int ns = (GetNumberOfRotationSteps() + 1) / 2;
  G4int np1 = G4int(dthe*ns/M_PI+.5) + 1;
  G4int np2 = rmin < perMillion ? 1 : np1;

  G4double *zz, *rr;
  zz = new G4double[np1+np2];
  rr = new G4double[np1+np2];

  G4double a = dthe/(np1-1);
  G4double cosa, sina;
  for (G4int i=0; i<np1; i++) {
    cosa  = cos(the+i*a);
    sina  = sin(the+i*a);
    zz[i] = rmax*cosa;
    rr[i] = rmax*sina;
    if (np2 > 1) {
      zz[i+np1] = rmin*cosa;
      rr[i+np1] = rmin*sina;
    }
  }
  if (np2 == 1) {
    zz[np1] = 0.;
    rr[np1] = 0.;
  }

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, phi, dphi, np1, np2, zz, rr, -1, 1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

G4PolyhedronTorus::G4PolyhedronTorus(G4double rmin, G4double rmax,
				     G4double rtor,
                                     G4double phi, G4double dphi)
/***********************************************************************
 *                                                                     *
 * Name: G4PolyhedronTorus                           Date:    11.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Constructor of polyhedron for TORUS                       *
 *                                                                     *
 * Input: rmin - internal radius                                       *
 *        rmax - external radius                                       *
 *        rtor - radius of torus                                       *
 *        phi  - initial phi                                           *
 *        dphi - delta phi                                             *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  if (dphi <= 0. || dphi > 2*M_PI) {
    G4cerr << "G4PolyhedronTorus: wrong delta phi = " << dphi << endl;
    return;
  }

  if (rmin < 0. || rmin >= rmax || rmax >= rtor) {
    G4cerr << "G4PolyhedronTorus: error in radiuses"
         << " rmin=" << rmin << " rmax=" << rmax
         << " rtorus=" << rtor << endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4int np1 = GetNumberOfRotationSteps();
  G4int np2 = rmin < perMillion ? 1 : np1;

  G4double *zz, *rr;
  zz = new G4double[np1+np2];
  rr = new G4double[np1+np2];

  G4double a = 2*M_PI/np1;
  G4double cosa, sina;
  for (G4int i=0; i<np1; i++) {
    cosa  = cos(i*a);
    sina  = sin(i*a);
    zz[i] = rmax*cosa;
    rr[i] = rtor+rmax*sina;
    if (np2 > 1) {
      zz[i+np1] = rmin*cosa;
      rr[i+np1] = rtor+rmin*sina;
    }
  }
  if (np2 == 1) {
    zz[np1] = 0.;
    rr[np1] = rtor;
    np2 = -1;
  }

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, phi, dphi, -np1, -np2, zz, rr, 1,-1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

G4int G4Polyhedron::fNumberOfRotationSteps = 24;
/***********************************************************************
 *                                                                     *
 * Name: G4Polyhedron::fNumberOfRotationSteps        Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Number of steps for whole circle                          *
 *                                                                     *
 ***********************************************************************/
