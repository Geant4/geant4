// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: HepPolyhedron.cc,v 1.3 2000-04-04 13:35:30 evc Exp $
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
//
// 19.03.00 E.Chernyaev
// - implemented boolean operations (add, subtract, intersect) on polyhedra;
//
  
#include "HepPolyhedron.h"

/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron operator <<                   Date:    09.05.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Print contents of G4 polyhedron                           *
 *                                                                     *
 ***********************************************************************/
std::ostream & operator<<(std::ostream & ostr, const G4Facet & facet) {
  for (int k=0; k<4; k++) {
    ostr << " " << facet.edge[k].v << "/" << facet.edge[k].f;
  }
  return ostr;
}

std::ostream & operator<<(std::ostream & ostr, const HepPolyhedron & ph) {
  ostr << std::endl;
  ostr << "Nverteces=" << ph.nvert << ", Nfacets=" << ph.nface << std::endl;
  int i;
  for (i=1; i<=ph.nvert; i++) {
     ostr << "xyz(" << i << ")="
	  << ph.pV[i].x() << ' ' << ph.pV[i].y() << ' ' << ph.pV[i].z()
	  << std::endl;
  }
  for (i=1; i<=ph.nface; i++) {
    ostr << "face(" << i << ")=" << ph.pF[i] << std::endl;
  }
  return ostr;
}

HepPolyhedron::HepPolyhedron(const HepPolyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron copy constructor             Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 ***********************************************************************/
{
  if (from.nvert > 0 && from.nface > 0) {
    nvert = from.nvert;
    nface = from.nface;
    pV = new HepPoint3D[nvert + 1];
    pF = new G4Facet[nface + 1];
    int i;
    for (i=1; i<=nvert; i++) pV[i] = from.pV[i];
    for (i=1; i<=nface; i++) pF[i] = from.pF[i];
  }else{
    nvert = 0; nface = 0; pV = 0; pF = 0;
  }
}

HepPolyhedron & HepPolyhedron::operator=(const HepPolyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron operator =                   Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 * Function: Copy contents of one GEANT4 polyhedron to another         *
 *                                                                     *
 ***********************************************************************/
{
  if (this == &from) return *this;
  delete [] pV;
  delete [] pF;
  if (from.nvert > 0  && from.nface > 0) {
    nvert = from.nvert;
    nface = from.nface;
    pV = new HepPoint3D[nvert + 1];
    pF = new G4Facet[nface + 1];
    int i;
    for (i=1; i<=nvert; i++) pV[i] = from.pV[i];
    for (i=1; i<=nface; i++) pF[i] = from.pF[i];
  }else{
    nvert = 0; nface = 0; pV = 0; pF = 0;
  }
  return *this;
}

int
HepPolyhedron::FindNeighbour(int iFace, int iNode, int iOrder) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::FindNeighbour                Date:    22.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Find neighbouring face                                    *
 *                                                                     *
 ***********************************************************************/
{
  int i;
  for (i=0; i<4; i++) {
    if (iNode == abs(pF[iFace].edge[i].v)) break;
  }
  if (i == 4) {
    std::cerr
      << "HepPolyhedron::FindNeighbour: face " << iFace
      << " has no node " << iNode
      << std::endl; 
    return 0;
  }
  if (iOrder < 0) {
    if ( --i < 0) i = 3;
    if (pF[iFace].edge[i].v == 0) i = 2;
  }
  return (pF[iFace].edge[i].v > 0) ? 0 : pF[iFace].edge[i].f;
}

HepNormal3D HepPolyhedron::FindNodeNormal(int iFace, int iNode) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::FindNodeNormal               Date:    22.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Find normal at given node                                 *
 *                                                                     *
 ***********************************************************************/
{
  HepNormal3D  normal = GetUnitNormal(iFace);
  int          k = iFace, iOrder = 1, n = 1;

  for(;;) {
    k = FindNeighbour(k, iNode, iOrder);
    if (k == iFace) break; 
    if (k > 0) {
      n++;
      normal += GetUnitNormal(k);
    }else{
      if (iOrder < 0) break;
      k = iFace;
      iOrder = -iOrder;
    }
  }
  return normal.unit();
}

void HepPolyhedron::SetNumberOfRotationSteps(int n)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetNumberOfRotationSteps     Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Set number of steps for whole circle                      *
 *                                                                     *
 ***********************************************************************/
{
  const int nMin = 3;
  if (n < nMin) {
    std::cerr 
      << "HepPolyhedron::SetNumberOfRotationSteps: attempt to set the\n"
      << "number of steps per circle < " << nMin << "; forced to " << nMin
      << std::endl;
    fNumberOfRotationSteps = nMin;
  }else{
    fNumberOfRotationSteps = n;
  }    
}

void HepPolyhedron::AllocateMemory(int Nvert, int Nface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::AllocateMemory               Date:    19.06.96 *
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
  pV    = new HepPoint3D[nvert+1];
  pF    = new G4Facet[nface+1];
}

void HepPolyhedron::CreatePrism()
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::CreatePrism                  Date:    15.07.96 *
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

void HepPolyhedron::RotateEdge(int k1, int k2, HepDouble r1, HepDouble r2,
			      int v1, int v2, int vEdge,
                              HepBoolean ifWholeCircle, int ns, int &kface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::RotateEdge                   Date:    05.12.96 *
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

  int i;
  int i1  = k1;
  int i2  = k2;
  int ii1 = ifWholeCircle ? i1 : i1+ns;
  int ii2 = ifWholeCircle ? i2 : i2+ns;
  int vv  = ifWholeCircle ? vEdge : 1;

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

void HepPolyhedron::SetSideFacets(int ii[4], int vv[4],
				 int *kk, HepDouble *r,
                                 HepDouble dphi, int ns, int &kface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetSideFacets                Date:    20.05.97 *
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
  int k1, k2, k3, k4;
  
  if (abs((HepDouble)(dphi-M_PI)) < perMillion) {          // half a circle
    for (int i=0; i<4; i++) {
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

void HepPolyhedron::RotateAroundZ(int nstep, HepDouble phi, HepDouble dphi,
                                 int np1, int np2,
				 const HepDouble *z, HepDouble *r,
				 int nodeVis, int edgeVis)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::RotateAroundZ                Date:    27.11.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Create HepPolyhedron for a solid produced by rotation of  *
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
  static HepDouble wholeCircle   = 2*M_PI;
    
  //   S E T   R O T A T I O N   P A R A M E T E R S

  HepBoolean ifWholeCircle = (abs(dphi-wholeCircle) < perMillion) ?
    true : false;
  HepDouble   delPhi  = ifWholeCircle ? wholeCircle : dphi;  
  int        nSphi    = (nstep > 0) ?
    nstep : int(delPhi*GetNumberOfRotationSteps()/wholeCircle+.5);
  if (nSphi == 0) nSphi = 1;
  int        nVphi    = ifWholeCircle ? nSphi : nSphi+1;
  HepBoolean ifClosed = np1 > 0 ? false : true;
  
  //   C O U N T   V E R T E C E S

  int absNp1 = abs(np1);
  int absNp2 = abs(np2);
  int i1beg = 0;
  int i1end = absNp1-1;
  int i2beg = absNp1;
  int i2end = absNp1+absNp2-1; 
  int i, j, k;

  for(i=i1beg; i<=i2end; i++) {
    if (abs(r[i]) < perMillion) r[i] = 0.;
  }

  j = 0;                                                // external nodes
  for (i=i1beg; i<=i1end; i++) {
    j += (r[i] == 0.) ? 1 : nVphi;
  }

  HepBoolean ifSide1 = false;                           // internal nodes
  HepBoolean ifSide2 = false;

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

  int *kk;
  kk = new int[absNp1+absNp2];

  k = 1;
  for(i=i1beg; i<=i1end; i++) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = HepPoint3D(0, 0, z[i]); } else { k += nVphi; }
  }

  i = i2beg;
  if (ifSide1) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = HepPoint3D(0, 0, z[i]); } else { k += nVphi; }
  }else{
    kk[i] = kk[i1beg];
  }

  for(i=i2beg+1; i<i2end; i++) {
    kk[i] = k;
    if (r[i] == 0.) { pV[k++] = HepPoint3D(0, 0, z[i]); } else { k += nVphi; }
  }

  if (absNp2 > 1) {
    i = i2end;
    if (ifSide2) {
      kk[i] = k;
      if (r[i] == 0.) pV[k] = HepPoint3D(0, 0, z[i]);
    }else{
      kk[i] = kk[i1end];
    }
  }

  HepDouble cosPhi, sinPhi;

  for(j=0; j<nVphi; j++) {
    cosPhi = cos(phi+j*delPhi/nSphi);
    sinPhi = sin(phi+j*delPhi/nSphi);
    for(i=i1beg; i<=i2end; i++) {
      if (r[i] != 0.) pV[kk[i]+j] = HepPoint3D(r[i]*cosPhi,r[i]*sinPhi,z[i]);
    }
  }

  //   G E N E R A T E   E X T E R N A L   F A C E S

  int v1,v2;

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

    int  ii[4], vv[4];

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
    std::cerr
      << "Polyhedron::RotateAroundZ: number of generated faces ("
      << k-1 << ") is not equal to the number of allocated faces ("
      << nface << ")"
      << std::endl;
  }
}

void HepPolyhedron::SetReferences()
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetReferences                Date:    04.12.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: For each edge set reference to neighbouring facet         *
 *                                                                     *
 ***********************************************************************/
{
  if (nface <= 0) return;

  struct edgeListMember {
    edgeListMember *next;
    int v2;
    int iface;
    int iedge;
  } *edgeList, *freeList, **headList;

  
  //   A L L O C A T E   A N D   I N I T I A T E   L I S T S

  edgeList = new edgeListMember[2*nface];
  headList = new edgeListMember*[nvert];
  
  int i;
  for (i=0; i<nvert; i++) {
    headList[i] = 0;
  }
  freeList = edgeList;
  for (i=0; i<2*nface-1; i++) {
    edgeList[i].next = &edgeList[i+1];
  }
  edgeList[2*nface-1].next = 0;

  //   L O O P   A L O N G   E D G E S

  int iface, iedge, nedge, i1, i2, k1, k2;
  edgeListMember *prev, *cur;
  
  for(iface=1; iface<=nface; iface++) {
    nedge = (pF[iface].edge[3].v == 0) ? 3 : 4;
    for (iedge=0; iedge<nedge; iedge++) {
      i1 = iedge;
      i2 = (iedge < nedge-1) ? iedge+1 : 0;
      i1 = abs(pF[iface].edge[i1].v);
      i2 = abs(pF[iface].edge[i2].v);
      k1 = (i1 < i2) ? i1 : i2;          // k1 = ::min(i1,i2);
      k2 = (i1 > i2) ? i1 : i2;          // k2 = ::max(i1,i2);
      
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
	  std::cerr
	    << "Polyhedron::SetReferences: different edge visibility "
	    << iface << "/" << iedge << "/"
	    << pF[iface].edge[iedge].v << " and "
	    << cur->iface << "/" << cur->iedge << "/"
	    << pF[cur->iface].edge[cur->iedge].v
	    << std::endl;
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
	      std::cerr
		<< "Polyhedron::SetReferences: different edge visibility "
		<< iface << "/" << iedge << "/"
		<< pF[iface].edge[iedge].v << " and "
		<< cur->iface << "/" << cur->iedge << "/"
		<< pF[cur->iface].edge[cur->iedge].v
		<< std::endl;
	    }
	  break;
	}
      }
    }
  }

  //  C H E C K   T H A T   A L L   L I S T S   A R E   E M P T Y

  for (i=0; i<nvert; i++) {
    if (headList[i] != 0) {
      std::cerr
	<< "Polyhedron::SetReferences: List " << i << " is not empty"
	<< std::endl;
    }
  }

  //   F R E E   M E M O R Y

  delete [] edgeList;
  delete [] headList;
}

void HepPolyhedron::InvertFacets()
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::InvertFacets                Date:    01.12.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Invert the order of the nodes in the facets               *
 *                                                                     *
 ***********************************************************************/
{
  if (nface <= 0) return;
  int i, k, nnode, v[4],f[4];
  for (i=1; i<=nface; i++) {
    nnode =  (pF[i].edge[3].v == 0) ? 3 : 4;
    for (k=0; k<nnode; k++) {
      v[k] = (k+1 == nnode) ? pF[i].edge[0].v : pF[i].edge[k+1].v;
      if (v[k] * pF[i].edge[k].v < 0) v[k] = -v[k];
      f[k] = pF[i].edge[k].f;
    }
    for (k=0; k<nnode; k++) {
      pF[i].edge[nnode-1-k].v = v[k];
      pF[i].edge[nnode-1-k].f = f[k];
    }
  }
}

HepPolyhedron & HepPolyhedron::Transform(const HepTransform3D &t)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::Transform                    Date:    01.12.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Make transformation of the polyhedron                     *
 *                                                                     *
 ***********************************************************************/
{
  if (nvert > 0) {
    for (int i=1; i<=nvert; i++) { pV[i] = t * pV[i]; }

    //  C H E C K   D E T E R M I N A N T   A N D
    //  I N V E R T   F A C E T S   I F   I T   I S   N E G A T I V E

    HepVector3D d = t * HepVector3D(0,0,0);
    HepVector3D x = t * HepVector3D(1,0,0) - d;
    HepVector3D y = t * HepVector3D(0,1,0) - d;
    HepVector3D z = t * HepVector3D(0,0,1) - d;
    if ((x.cross(y))*z < 0) InvertFacets();
  }
  return *this;
}

HepBoolean HepPolyhedron::GetNextVertexIndex(int &index, int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextVertexIndex          Date:    03.09.96  *
 * Author: Yasuhide Sawada                          Revised:           *
 *                                                                     *
 * Function:                                                           *
 *                                                                     *
 ***********************************************************************/
{
  static int iFace = 1;
  static int iQVertex = 0;
  int vIndex = pF[iFace].edge[iQVertex].v;

  edgeFlag = (vIndex > 0) ? 1 : 0;
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

HepPoint3D HepPolyhedron::GetVertex(int index) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetVertex                   Date:    03.09.96  *
 * Author: Yasuhide Sawada                          Revised: 17.11.99  *
 *                                                                     *
 * Function: Get vertex of the index.                                  *
 *                                                                     *
 ***********************************************************************/
{
  if (index <= 0 || index > nvert) {
    std::cerr
      << "HepPolyhedron::GetVertex: irrelevant index " << index
      << std::endl;
    return HepPoint3D();
  }
  return pV[index];
}

HepBoolean
HepPolyhedron::GetNextVertex(HepPoint3D &vertex, int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextVertex               Date:    22.07.96  *
 * Author: John Allison                             Revised:           *
 *                                                                     *
 * Function: Get vertices of the quadrilaterals in order for each      *
 *           face in face order.  Returns false when finished each     *
 *           face.                                                     *
 *                                                                     *
 ***********************************************************************/
{
  int index;
  HepBoolean rep = GetNextVertexIndex(index, edgeFlag);
  vertex = pV[index];
  return rep;
}

HepBoolean HepPolyhedron::GetNextVertex(HepPoint3D &vertex, int &edgeFlag,
				       HepNormal3D &normal) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextVertex               Date:    26.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get vertices with normals of the quadrilaterals in order  *
 *           for each face in face order.                              *
 *           Returns false when finished each face.                    *
 *                                                                     *
 ***********************************************************************/
{
  static int iFace = 1;
  static int iNode = 0;

  if (nface == 0) return false;  // empty polyhedron

  int k = pF[iFace].edge[iNode].v;
  if (k > 0) { edgeFlag = 1; } else { edgeFlag = -1; k = -k; }
  vertex = pV[k];
  normal = FindNodeNormal(iFace,k);
  if (iNode >= 3 || pF[iFace].edge[iNode+1].v == 0) {
    iNode = 0;
    if (++iFace > nface) iFace = 1;
    return false;                // last node
  }else{
    ++iNode;
    return true;                 // not last node
  }
}

HepBoolean HepPolyhedron::GetNextEdgeIndeces(int &i1, int &i2, int &edgeFlag,
					    int &iface1, int &iface2) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdgeIndeces          Date:    30.09.96  *
 * Author: E.Chernyaev                              Revised: 17.11.99  *
 *                                                                     *
 * Function: Get indeces of the next edge together with indeces of     *
 *           of the faces which share the edge.                        *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  static int iFace    = 1;
  static int iQVertex = 0;
  static int iOrder   = 1;
  int  k1, k2, kflag, kface1, kface2;

  if (iFace == 1 && iQVertex == 0) {
    k2 = pF[nface].edge[0].v;
    k1 = pF[nface].edge[3].v;
    if (k1 == 0) k1 = pF[nface].edge[2].v;
    if (abs(k1) > abs(k2)) iOrder = -1;
  }

  do {
    k1     = pF[iFace].edge[iQVertex].v;
    kflag  = k1;
    k1     = abs(k1);
    kface1 = iFace; 
    kface2 = pF[iFace].edge[iQVertex].f;
    if (iQVertex >= 3 || pF[iFace].edge[iQVertex+1].v == 0) {
      iQVertex = 0;
      k2 = abs(pF[iFace].edge[iQVertex].v);
      iFace++;
    }else{
      iQVertex++;
      k2 = abs(pF[iFace].edge[iQVertex].v);
    }
  } while (iOrder*k1 > iOrder*k2);

  i1 = k1; i2 = k2; edgeFlag = (kflag > 0) ? 1 : 0;
  iface1 = kface1; iface2 = kface2; 

  if (iFace > nface) {
    iFace  = 1; iOrder = 1;
    return false;
  }else{
    return true;
  }
}

HepBoolean
HepPolyhedron::GetNextEdgeIndeces(int &i1, int &i2, int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdgeIndeces          Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get indeces of the next edge.                             *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  int kface1, kface2;
  return GetNextEdgeIndeces(i1, i2, edgeFlag, kface1, kface2);
}

HepBoolean
HepPolyhedron::GetNextEdge(HepPoint3D &p1,
			   HepPoint3D &p2,
			   int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdge                 Date:    30.09.96  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get next edge.                                            *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  int i1,i2;
  HepBoolean rep = GetNextEdgeIndeces(i1,i2,edgeFlag);
  p1 = pV[i1];
  p2 = pV[i2];
  return rep;
}

HepBoolean
HepPolyhedron::GetNextEdge(HepPoint3D &p1, HepPoint3D &p2,
			  int &edgeFlag, int &iface1, int &iface2) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdge                 Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get next edge with indeces of the faces which share       *
 *           the edge.                                                 *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  int i1,i2;
  HepBoolean rep = GetNextEdgeIndeces(i1,i2,edgeFlag,iface1,iface2);
  p1 = pV[i1];
  p2 = pV[i2];
  return rep;
}

void HepPolyhedron::GetFacet(int iFace, int &n, int *iNodes,
			    int *edgeFlags, int *iFaces) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetFacet                    Date:    15.12.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get face by index                                         *
 *                                                                     *
 ***********************************************************************/
{
  if (iFace < 1 || iFace > nface) {
    std::cerr 
      << "HepPolyhedron::GetFacet: irrelevant index " << iFace
      << std::endl;
    n = 0;
  }else{
    int i, k;
    for (i=0; i<4; i++) { 
      k = pF[iFace].edge[i].v;
      if (k == 0) break;
      if (iFaces != 0) iFaces[i] = pF[iFace].edge[i].f;
      if (k > 0) { 
	iNodes[i] = k;
	if (edgeFlags != 0) edgeFlags[i] = 1;
      }else{
	iNodes[i] = -k;
	if (edgeFlags != 0) edgeFlags[i] = -1;
      }
    }
    n = i;
  }
}

void HepPolyhedron::GetFacet(int index, int &n, HepPoint3D *nodes,
			    int *edgeFlags, HepNormal3D *normals) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetFacet                    Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get face by index                                         *
 *                                                                     *
 ***********************************************************************/
{
  int iNodes[4];
  GetFacet(index, n, iNodes, edgeFlags);
  if (n != 0) {
    for (int i=0; i<4; i++) {
      nodes[i] = pV[iNodes[i]];
      if (normals != 0) normals[i] = FindNodeNormal(index,iNodes[i]);
    }
  }
}

HepBoolean
HepPolyhedron::GetNextFacet(int &n, HepPoint3D *nodes,
			   int *edgeFlags, HepNormal3D *normals) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextFacet                Date:    19.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get next face with normals of unit length at the nodes.   *
 *           Returns false when finished all faces.                    *
 *                                                                     *
 ***********************************************************************/
{
  static int iFace = 1;

  if (edgeFlags == 0) {
    GetFacet(iFace, n, nodes);
  }else if (normals == 0) {
    GetFacet(iFace, n, nodes, edgeFlags);
  }else{
    GetFacet(iFace, n, nodes, edgeFlags, normals);
  }

  if (++iFace > nface) {
    iFace  = 1;
    return false;
  }else{
    return true;
  }
}

HepNormal3D HepPolyhedron::GetNormal(int iFace) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNormal                    Date:    19.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Get normal of the face given by index                     *
 *                                                                     *
 ***********************************************************************/
{
  if (iFace < 1 || iFace > nface) {
    std::cerr 
      << "HepPolyhedron::GetNormal: irrelevant index " << iFace 
      << std::endl;
    return HepNormal3D();
  }

  int i0  = abs(pF[iFace].edge[0].v);
  int i1  = abs(pF[iFace].edge[1].v);
  int i2  = abs(pF[iFace].edge[2].v);
  int i3  = abs(pF[iFace].edge[3].v);
  if (i3 == 0) i3 = i0;
  return (pV[i2] - pV[i0]).cross(pV[i3] - pV[i1]);
}

HepNormal3D HepPolyhedron::GetUnitNormal(int iFace) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNormal                    Date:    19.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Get unit normal of the face given by index                *
 *                                                                     *
 ***********************************************************************/
{
  if (iFace < 1 || iFace > nface) {
    std::cerr 
      << "HepPolyhedron::GetUnitNormal: irrelevant index " << iFace
      << std::endl;
    return HepNormal3D();
  }

  int i0  = abs(pF[iFace].edge[0].v);
  int i1  = abs(pF[iFace].edge[1].v);
  int i2  = abs(pF[iFace].edge[2].v);
  int i3  = abs(pF[iFace].edge[3].v);
  if (i3 == 0) i3 = i0;
  return ((pV[i2] - pV[i0]).cross(pV[i3] - pV[i1])).unit();
}

HepBoolean HepPolyhedron::GetNextNormal(HepNormal3D &normal) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextNormal               Date:    22.07.96  *
 * Author: John Allison                             Revised: 19.11.99  *
 *                                                                     *
 * Function: Get normals of each face in face order.  Returns false    *
 *           when finished all faces.                                  *
 *                                                                     *
 ***********************************************************************/
{
  static int iFace = 1;
  normal = GetNormal(iFace);
  if (++iFace > nface) {
    iFace = 1;
    return false;
  }else{
    return true;
  }
}

HepBoolean HepPolyhedron::GetNextUnitNormal(HepNormal3D &normal) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextUnitNormal           Date:    16.09.96  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get normals of unit length of each face in face order.    *
 *           Returns false when finished all faces.                    *
 *                                                                     *
 ***********************************************************************/
{
  HepBoolean rep = GetNextNormal(normal);
  normal = normal.unit();
  return rep;
}

HepPolyhedronTrd2::HepPolyhedronTrd2(HepDouble Dx1, HepDouble Dx2,
				     HepDouble Dy1, HepDouble Dy2,
				     HepDouble Dz)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronTrd2                           Date:    22.07.96 *
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

  pV[1] = HepPoint3D(-Dx1,-Dy1,-Dz);
  pV[2] = HepPoint3D( Dx1,-Dy1,-Dz);
  pV[3] = HepPoint3D( Dx1, Dy1,-Dz);
  pV[4] = HepPoint3D(-Dx1, Dy1,-Dz);
  pV[5] = HepPoint3D(-Dx2,-Dy2, Dz);
  pV[6] = HepPoint3D( Dx2,-Dy2, Dz);
  pV[7] = HepPoint3D( Dx2, Dy2, Dz);
  pV[8] = HepPoint3D(-Dx2, Dy2, Dz);

  CreatePrism();
}

HepPolyhedronTrap::HepPolyhedronTrap(HepDouble Dz,
				     HepDouble Theta,
				     HepDouble Phi,
				     HepDouble Dy1,
				     HepDouble Dx1,
				     HepDouble Dx2,
				     HepDouble Alp1,
				     HepDouble Dy2,
				     HepDouble Dx3,
				     HepDouble Dx4,
				     HepDouble Alp2)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronTrap                           Date:    20.11.96 *
 * Author: E.Chernyaev                               Revised:          *
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
  HepDouble DzTthetaCphi = Dz*tan(Theta)*cos(Phi);
  HepDouble DzTthetaSphi = Dz*tan(Theta)*sin(Phi);
  HepDouble Dy1Talp1 = Dy1*tan(Alp1);
  HepDouble Dy2Talp2 = Dy2*tan(Alp2);
  
  AllocateMemory(8,6);

  pV[1] = HepPoint3D(-DzTthetaCphi-Dy1Talp1-Dx1,-DzTthetaSphi-Dy1,-Dz);
  pV[2] = HepPoint3D(-DzTthetaCphi-Dy1Talp1+Dx1,-DzTthetaSphi-Dy1,-Dz);
  pV[3] = HepPoint3D(-DzTthetaCphi+Dy1Talp1+Dx2,-DzTthetaSphi+Dy1,-Dz);
  pV[4] = HepPoint3D(-DzTthetaCphi+Dy1Talp1-Dx2,-DzTthetaSphi+Dy1,-Dz);
  pV[5] = HepPoint3D( DzTthetaCphi-Dy2Talp2-Dx3, DzTthetaSphi-Dy2, Dz);
  pV[6] = HepPoint3D( DzTthetaCphi-Dy2Talp2+Dx3, DzTthetaSphi-Dy2, Dz);
  pV[7] = HepPoint3D( DzTthetaCphi+Dy2Talp2+Dx4, DzTthetaSphi+Dy2, Dz);
  pV[8] = HepPoint3D( DzTthetaCphi+Dy2Talp2-Dx4, DzTthetaSphi+Dy2, Dz);

  CreatePrism();
}

HepPolyhedronCons::HepPolyhedronCons(HepDouble Rmn1,
				     HepDouble Rmx1,
				     HepDouble Rmn2,
				     HepDouble Rmx2, 
				     HepDouble Dz,
				     HepDouble Phi1,
				     HepDouble Dphi) 
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronCons::HepPolyhedronCons        Date:    15.12.96 *
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
  static HepDouble wholeCircle=2*M_PI;

  //   C H E C K   I N P U T   P A R A M E T E R S

  int k = 0;
  if (Rmn1 < 0. || Rmx1 < 0. || Rmn2 < 0. || Rmx2 < 0.)        k = 1;
  if (Rmn1 > Rmx1 || Rmn2 > Rmx2)                              k = 1;
  if (Rmn1 == Rmx1 && Rmn2 == Rmx2)                            k = 1;

  if (Dz <= 0.) k += 2;
 
  HepDouble phi1, phi2, dphi;
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
    std::cerr << "HepPolyhedronCone(s)/Tube(s): error in input parameters";
    if ((k & 1) != 0) std::cerr << " (radiuses)";
    if ((k & 2) != 0) std::cerr << " (half-length)";
    if ((k & 4) != 0) std::cerr << " (angles)";
    std::cerr << std::endl;
    std::cerr << " Rmn1=" << Rmn1 << " Rmx1=" << Rmx1;
    std::cerr << " Rmn2=" << Rmn2 << " Rmx2=" << Rmx2;
    std::cerr << " Dz=" << Dz << " Phi1=" << Phi1 << " Dphi=" << Dphi
	      << std::endl;
    return;
  }
  
  //   P R E P A R E   T W O   P O L Y L I N E S

  HepDouble zz[4], rr[4];
  zz[0] =  Dz; 
  zz[1] = -Dz; 
  zz[2] =  Dz; 
  zz[3] = -Dz; 
  rr[0] =  Rmx2;
  rr[1] =  Rmx1;
  rr[2] =  Rmn2;
  rr[3] =  Rmn1;

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, phi1, dphi, 2, 2, zz, rr, -1, -1); 
  SetReferences();
}

HepPolyhedronPgon::HepPolyhedronPgon(HepDouble phi,
				     HepDouble dphi,
				     int    npdv,
				     int    nz,
				     const HepDouble *z,
				     const HepDouble *rmin,
				     const HepDouble *rmax)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronPgon                           Date:    09.12.96 *
 * Author: E.Chernyaev                               Revised:          *
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
    std::cerr
      << "HepPolyhedronPgon/Pcon: wrong delta phi = " << dphi
      << std::endl;
    return;
  }    
    
  if (nz < 2) {
    std::cerr
      << "HepPolyhedronPgon/Pcon: number of z-planes less than two = " << nz
      << std::endl;
    return;
  }

  if (npdv < 0) {
    std::cerr
      << "HepPolyhedronPgon/Pcon: error in number of phi-steps =" << npdv
      << std::endl;
    return;
  }

  int i;
  for (i=0; i<nz; i++) {
    if (rmin[i] < 0. || rmax[i] < 0. || rmin[i] > rmax[i]) {
      std::cerr
	<< "HepPolyhedronPgon: error in radiuses rmin[" << i << "]="
	<< rmin[i] << " rmax[" << i << "]=" << rmax[i]
	<< std::endl;
      return;
    }
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  HepDouble *zz, *rr;
  zz = new HepDouble[2*nz];
  rr = new HepDouble[2*nz];

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

  RotateAroundZ(npdv, phi, dphi, nz, nz, zz, rr, -1, (npdv == 0) ? -1 : 1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

HepPolyhedronSphere::HepPolyhedronSphere(HepDouble rmin, HepDouble rmax,
				       HepDouble phi, HepDouble dphi,
				       HepDouble the, HepDouble dthe)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronSphere                         Date:    11.12.96 *
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
    std::cerr
      << "HepPolyhedronSphere: wrong delta phi = " << dphi
      << std::endl;
    return;
  }    

  if (the < 0. || the > M_PI) {
    std::cerr
      << "HepPolyhedronSphere: wrong theta = " << the
      << std::endl;
    return;
  }    
  
  if (dthe <= 0. || dthe > M_PI) {
    std::cerr
      << "HepPolyhedronSphere: wrong delta theta = " << dthe
      << std::endl;
    return;
  }    

  if (the+dthe > M_PI) {
    std::cerr
      << "HepPolyhedronSphere: wrong theta + delta theta = "
      << the << " " << dthe
      << std::endl;
    return;
  }    
  
  if (rmin < 0. || rmin >= rmax) {
    std::cerr
      << "HepPolyhedronSphere: error in radiuses"
      << " rmin=" << rmin << " rmax=" << rmax
      << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  int ns = (GetNumberOfRotationSteps() + 1) / 2;
  int np1 = int(dthe*ns/M_PI+.5) + 1;
  if (np1 <= 1) np1 = 2;
  int np2 = rmin < perMillion ? 1 : np1;

  HepDouble *zz, *rr;
  zz = new HepDouble[np1+np2];
  rr = new HepDouble[np1+np2];

  HepDouble a = dthe/(np1-1);
  HepDouble cosa, sina;
  for (int i=0; i<np1; i++) {
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

  RotateAroundZ(0, phi, dphi, np1, np2, zz, rr, -1, -1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

HepPolyhedronTorus::HepPolyhedronTorus(HepDouble rmin,
				       HepDouble rmax,
				       HepDouble rtor,
				       HepDouble phi,
				       HepDouble dphi)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronTorus                          Date:    11.12.96 *
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
    std::cerr
      << "HepPolyhedronTorus: wrong delta phi = " << dphi
      << std::endl;
    return;
  }

  if (rmin < 0. || rmin >= rmax || rmax >= rtor) {
    std::cerr
      << "HepPolyhedronTorus: error in radiuses"
      << " rmin=" << rmin << " rmax=" << rmax << " rtorus=" << rtor
      << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  int np1 = GetNumberOfRotationSteps();
  int np2 = rmin < perMillion ? 1 : np1;

  HepDouble *zz, *rr;
  zz = new HepDouble[np1+np2];
  rr = new HepDouble[np1+np2];

  HepDouble a = 2*M_PI/np1;
  HepDouble cosa, sina;
  for (int i=0; i<np1; i++) {
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

  RotateAroundZ(0, phi, dphi, -np1, -np2, zz, rr, -1,-1); 
  SetReferences();
  
  delete [] zz;
  delete [] rr;
}

int HepPolyhedron::fNumberOfRotationSteps = DEFAULT_NUMBER_OF_STEPS;
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::fNumberOfRotationSteps       Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Number of steps for whole circle                          *
 *                                                                     *
 ***********************************************************************/

#include "BooleanProcessor.src"
static BooleanProcessor processor;

HepPolyhedron HepPolyhedron::add(const HepPolyhedron & p) const 
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::add                          Date:    19.03.00 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Boolean "union" of two polyhedra                          *
 *                                                                     *
 ***********************************************************************/
{
  return processor.execute(OP_UNION, *this, p);
}

HepPolyhedron HepPolyhedron::intersect(const HepPolyhedron & p) const 
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::intersect                    Date:    19.03.00 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Boolean "intersection" of two polyhedra                   *
 *                                                                     *
 ***********************************************************************/
{
  return processor.execute(OP_INTERSECTION, *this, p);
}

HepPolyhedron HepPolyhedron::subtract(const HepPolyhedron & p) const 
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::add                          Date:    19.03.00 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Boolean "subtraction" of "p" from "this"                  *
 *                                                                     *
 ***********************************************************************/
{
  return processor.execute(OP_SUBTRACTION, *this, p);
}
