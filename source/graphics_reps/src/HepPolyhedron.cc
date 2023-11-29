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
// G4 Polyhedron library
//
// History:
// 23.07.96 E.Chernyaev <Evgueni.Tcherniaev@cern.ch> - initial version
//
// 30.09.96 E.Chernyaev
// - added GetNextVertexIndex, GetVertex by Yasuhide Sawada
// - added GetNextUnitNormal, GetNextEdgeIndices, GetNextEdge
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
// 25.05.01 E.Chernyaev
// - added GetSurfaceArea() and GetVolume()
//
// 05.11.02 E.Chernyaev
// - added createTwistedTrap() and createPolyhedron()
//
// 20.06.05 G.Cosmo
// - added HepPolyhedronEllipsoid
//
// 18.07.07 T.Nikitina
// - added HepPolyhedronParaboloid
//
// 22.02.20 E.Chernyaev
// - added HepPolyhedronTet, HepPolyhedronHyberbolicMirror
//
// 12.05.21 E.Chernyaev
// - added TriangulatePolygon(), RotateContourAroundZ()
// - added HepPolyhedronPgon, HepPolyhedronPcon given by rz-contour
//
// 26.03.22 E.Chernyaev
// - added SetVertex(), SetFacet()
// - added HepPolyhedronTetMesh
//
// 04.04.22 E.Chernyaev
// - added JoinCoplanarFacets()
//
// 07.04.22 E.Chernyaev
// - added HepPolyhedronBoxMesh

#include "HepPolyhedron.h"
#include "G4PhysicalConstants.hh"
#include "G4Vector3D.hh"

#include <cstdlib>  // Required on some compilers for std::abs(int) ...
#include <cmath>
#include <algorithm>

using CLHEP::perMillion;
using CLHEP::deg;
using CLHEP::pi;
using CLHEP::twopi;
using CLHEP::nm;
const G4double spatialTolerance = 0.01*nm;

/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron operator <<                   Date:    09.05.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Print contents of G4 polyhedron                           *
 *                                                                     *
 ***********************************************************************/
std::ostream & operator<<(std::ostream & ostr, const G4Facet & facet) {
  for (const auto& edge : facet.edge) {
    ostr << " " << edge.v << "/" << edge.f;
  }
  return ostr;
}

std::ostream & operator<<(std::ostream & ostr, const HepPolyhedron & ph) {
  ostr << std::endl;
  ostr << "Nvertices=" << ph.nvert << ", Nfacets=" << ph.nface << std::endl;
  G4int i;
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

HepPolyhedron::HepPolyhedron(G4int Nvert, G4int Nface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron constructor with           Date:    26.03.2022  *
 *       allocation of memory                     Revised:             *
 * Author: E.Tcherniaev (E.Chernyaev)                                  *
 *                                                                     *
 ***********************************************************************/
: nvert(0), nface(0), pV(nullptr), pF(nullptr)
{
  AllocateMemory(Nvert, Nface);
}

HepPolyhedron::HepPolyhedron(const HepPolyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron copy constructor             Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 ***********************************************************************/
: nvert(0), nface(0), pV(nullptr), pF(nullptr)
{
  AllocateMemory(from.nvert, from.nface);
  for (G4int i=1; i<=nvert; i++) pV[i] = from.pV[i];
  for (G4int k=1; k<=nface; k++) pF[k] = from.pF[k];
}

HepPolyhedron::HepPolyhedron(HepPolyhedron&& from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron move constructor           Date:    04.11.2019  *
 * Author: E.Tcherniaev (E.Chernyaev)             Revised:             *
 *                                                                     *
 ***********************************************************************/
: nvert(0), nface(0), pV(nullptr), pF(nullptr)
{
  nvert = from.nvert;
  nface = from.nface;
  pV = from.pV;
  pF = from.pF;

  // Release the data from the source object
  from.nvert = 0;
  from.nface = 0;
  from.pV = nullptr;
  from.pF = nullptr;
}

HepPolyhedron & HepPolyhedron::operator=(const HepPolyhedron &from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron operator =                   Date:    23.07.96  *
 * Author: E.Chernyaev (IHEP/Protvino)              Revised:           *
 *                                                                     *
 * Function: Copy contents of one polyhedron to another                *
 *                                                                     *
 ***********************************************************************/
{
  if (this != &from) {
    AllocateMemory(from.nvert, from.nface);
    for (G4int i=1; i<=nvert; i++) pV[i] = from.pV[i];
    for (G4int k=1; k<=nface; k++) pF[k] = from.pF[k];
  }
  return *this;
}

HepPolyhedron & HepPolyhedron::operator=(HepPolyhedron&& from)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron move operator =              Date:   04.11.2019 *
 * Author: E.Tcherniaev (E.Chernyaev)               Revised:           *
 *                                                                     *
 * Function: Move contents of one polyhedron to another                *
 *                                                                     *
 ***********************************************************************/
{
  if (this != &from) {
    delete [] pV;
    delete [] pF;
    nvert = from.nvert;
    nface = from.nface;
    pV = from.pV;
    pF = from.pF;

    // Release the data from the source object
    from.nvert = 0;
    from.nface = 0;
    from.pV = nullptr;
    from.pF = nullptr;
  }
  return *this;
}

G4int
HepPolyhedron::FindNeighbour(G4int iFace, G4int iNode, G4int iOrder) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::FindNeighbour                Date:    22.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Find neighbouring face                                    *
 *                                                                     *
 ***********************************************************************/
{
  G4int i;
  for (i=0; i<4; i++) {
    if (iNode == std::abs(pF[iFace].edge[i].v)) break;
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

G4Normal3D HepPolyhedron::FindNodeNormal(G4int iFace, G4int iNode) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::FindNodeNormal               Date:    22.11.99 *
 * Author: E.Chernyaev                               Revised:          *
 *                                                                     *
 * Function: Find normal at given node                                 *
 *                                                                     *
 ***********************************************************************/
{
  G4Normal3D normal = GetUnitNormal(iFace);
  G4int      k = iFace, iOrder = 1;

  for(;;) {
    k = FindNeighbour(k, iNode, iOrder);
    if (k == iFace) break;
    if (k > 0) {
      normal += GetUnitNormal(k);
    }else{
      if (iOrder < 0) break;
      k = iFace;
      iOrder = -iOrder;
    }
  }
  return normal.unit();
}

G4int HepPolyhedron::GetNumberOfRotationSteps()
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNumberOfRotationSteps     Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Get number of steps for whole circle                      *
 *                                                                     *
 ***********************************************************************/
{
  return fNumberOfRotationSteps;
}

void HepPolyhedron::SetVertex(G4int index, const G4Point3D& v)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetVertex                    Date:    26.03.22 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Set vertex                                                *
 *                                                                     *
 ***********************************************************************/
{
  if (index < 1 || index > nvert)
  {
    std::cerr
      << "HepPolyhedron::SetVertex: vertex index = " << index
      << " is out of range\n"
      << "   N. of vertices = " << nvert << "\n"
      << "   N. of facets = " << nface << std::endl;
    return;
  }
  pV[index] = v;
}

void
HepPolyhedron::SetFacet(G4int index, G4int iv1, G4int iv2, G4int iv3, G4int iv4)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetFacet                     Date:    26.03.22 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Set facet                                                *
 *                                                                     *
 ***********************************************************************/
{
  if (index < 1 || index > nface)
  {
    std::cerr
      << "HepPolyhedron::SetFacet: facet index = " << index
      << " is out of range\n"
      << "   N. of vertices = " << nvert << "\n"
      << "   N. of facets = " << nface << std::endl;
    return;
  }
  if (iv1 < 1 || iv1 > nvert ||
      iv2 < 1 || iv2 > nvert ||
      iv3 < 1 || iv3 > nvert ||
      iv4 < 0 || iv4 > nvert)
  {
    std::cerr
      << "HepPolyhedron::SetFacet: incorrectly specified facet"
      << " (" << iv1 << ", " << iv2 << ", " << iv3 << ", " << iv4 << ")\n"
      << "   N. of vertices = " << nvert << "\n"
      << "   N. of facets = " << nface << std::endl;
    return;
  }
  pF[index] = G4Facet(iv1, 0, iv2, 0, iv3, 0, iv4, 0);
}

void HepPolyhedron::SetNumberOfRotationSteps(G4int n)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetNumberOfRotationSteps     Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Set number of steps for whole circle                      *
 *                                                                     *
 ***********************************************************************/
{
  const G4int nMin = 3;
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

void HepPolyhedron::ResetNumberOfRotationSteps()
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNumberOfRotationSteps     Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Reset number of steps for whole circle to default value   *
 *                                                                     *
 ***********************************************************************/
{
  fNumberOfRotationSteps = DEFAULT_NUMBER_OF_STEPS;
}

void HepPolyhedron::AllocateMemory(G4int Nvert, G4int Nface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::AllocateMemory               Date:    19.06.96 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised: 05.11.02 *
 *                                                                     *
 * Function: Allocate memory for GEANT4 polyhedron                     *
 *                                                                     *
 * Input: Nvert - number of nodes                                      *
 *        Nface - number of faces                                      *
 *                                                                     *
 ***********************************************************************/
{
  if (nvert == Nvert && nface == Nface) return;
  delete [] pV;
  delete [] pF;
  if (Nvert > 0 && Nface > 0) {
    nvert = Nvert;
    nface = Nface;
    pV    = new G4Point3D[nvert+1];
    pF    = new G4Facet[nface+1];
  }else{
    nvert = 0; nface = 0; pV = nullptr; pF = nullptr;
  }
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

void HepPolyhedron::RotateEdge(G4int k1, G4int k2, G4double r1, G4double r2,
                              G4int v1, G4int v2, G4int vEdge,
                              G4bool ifWholeCircle, G4int nds, G4int &kface)
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
 *        nds    - number of discrete steps                            *
 *        r[]    - r-coordinates                                       *
 *        kface  - current free cell in the pF array                   *
 *                                                                     *
 ***********************************************************************/
{
  if (r1 == 0. && r2 == 0.) return;

  G4int i;
  G4int i1  = k1;
  G4int i2  = k2;
  G4int ii1 = ifWholeCircle ? i1 : i1+nds;
  G4int ii2 = ifWholeCircle ? i2 : i2+nds;
  G4int vv  = ifWholeCircle ? vEdge : 1;

  if (nds == 1) {
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
      for (i2++,i=1; i<nds-1; i2++,i++) {
        pF[kface++] = G4Facet(vEdge*i1,0, v2*i2,0, vEdge*(i2+1),0);
      }
      pF[kface++]   = G4Facet(vEdge*i1,0, v2*i2,0, vv*ii2,0);
    }else if (r2 == 0.) {
      pF[kface++]   = G4Facet(vv*i1,0,    vEdge*i2,0, v1*(i1+1),0);
      for (i1++,i=1; i<nds-1; i1++,i++) {
        pF[kface++] = G4Facet(vEdge*i1,0, vEdge*i2,0, v1*(i1+1),0);
      }
      pF[kface++]   = G4Facet(vEdge*i1,0, vv*i2,0,    v1*ii1,0);
    }else{
      pF[kface++]   = G4Facet(vv*i1,0,    v2*i2,0, vEdge*(i2+1),0,v1*(i1+1),0);
      for (i1++,i2++,i=1; i<nds-1; i1++,i2++,i++) {
        pF[kface++] = G4Facet(vEdge*i1,0, v2*i2,0, vEdge*(i2+1),0,v1*(i1+1),0);
      }
      pF[kface++]   = G4Facet(vEdge*i1,0, v2*i2,0, vv*ii2,0,      v1*ii1,0);
    }
  }
}

void HepPolyhedron::SetSideFacets(G4int ii[4], G4int vv[4],
                                 G4int *kk, G4double *r,
                                 G4double dphi, G4int nds, G4int &kface)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::SetSideFacets                Date:    20.05.97 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Set side facets for the case of incomplete rotation       *
 *                                                                     *
 * Input: ii[4] - indices of original vertices                         *
 *        vv[4] - visibility of edges                                  *
 *        kk[]  - indices of nodes                                     *
 *        r[]   - radiuses                                             *
 *        dphi  - delta phi                                            *
 *        nds    - number of discrete steps                            *
 *        kface  - current free cell in the pF array                   *
 *                                                                     *
 ***********************************************************************/
{
  G4int k1, k2, k3, k4;

  if (std::abs(dphi-pi) < perMillion) { // half a circle
    for (G4int i=0; i<4; i++) {
      k1 = ii[i];
      k2 = ii[(i+1)%4];
      if (r[k1] == 0. && r[k2] == 0.) vv[i] = -1;
    }
  }

  if (ii[1] == ii[2]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[2]];
    k3 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[2]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += nds;
    if (r[ii[2]] != 0.) k2 += nds;
    if (r[ii[3]] != 0.) k3 += nds;
    pF[kface++] = G4Facet(vv[2]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }else if (kk[ii[0]] == kk[ii[1]]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[2]];
    k3 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[1]*k1,0, vv[2]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += nds;
    if (r[ii[2]] != 0.) k2 += nds;
    if (r[ii[3]] != 0.) k3 += nds;
    pF[kface++] = G4Facet(vv[2]*k3,0, vv[1]*k2,0, vv[3]*k1,0);
  }else if (kk[ii[2]] == kk[ii[3]]) {
    k1 = kk[ii[0]];
    k2 = kk[ii[1]];
    k3 = kk[ii[2]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[1]*k2,0, vv[3]*k3,0);
    if (r[ii[0]] != 0.) k1 += nds;
    if (r[ii[1]] != 0.) k2 += nds;
    if (r[ii[2]] != 0.) k3 += nds;
    pF[kface++] = G4Facet(vv[1]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }else{
    k1 = kk[ii[0]];
    k2 = kk[ii[1]];
    k3 = kk[ii[2]];
    k4 = kk[ii[3]];
    pF[kface++] = G4Facet(vv[0]*k1,0, vv[1]*k2,0, vv[2]*k3,0, vv[3]*k4,0);
    if (r[ii[0]] != 0.) k1 += nds;
    if (r[ii[1]] != 0.) k2 += nds;
    if (r[ii[2]] != 0.) k3 += nds;
    if (r[ii[3]] != 0.) k4 += nds;
    pF[kface++] = G4Facet(vv[2]*k4,0, vv[1]*k3,0, vv[0]*k2,0, vv[3]*k1,0);
  }
}

void HepPolyhedron::RotateAroundZ(G4int nstep, G4double phi, G4double dphi,
                                 G4int np1, G4int np2,
                                 const G4double *z, G4double *r,
                                 G4int nodeVis, G4int edgeVis)
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
  static const G4double wholeCircle   = twopi;

  //   S E T   R O T A T I O N   P A R A M E T E R S

  G4bool ifWholeCircle = std::abs(dphi-wholeCircle) < perMillion;
  G4double delPhi = ifWholeCircle ? wholeCircle : dphi;
  G4int nSphi = nstep;
  if (nSphi <= 0) nSphi = GetNumberOfRotationSteps()*delPhi/wholeCircle + 0.5;
  if (nSphi == 0) nSphi = 1;
  G4int nVphi = ifWholeCircle ? nSphi : nSphi + 1;
  G4bool ifClosed = np1 <= 0; // true if external polyline is closed

  //   C O U N T   V E R T I C E S

  G4int absNp1 = std::abs(np1);
  G4int absNp2 = std::abs(np2);
  G4int i1beg = 0;
  G4int i1end = absNp1-1;
  G4int i2beg = absNp1;
  G4int i2end = absNp1+absNp2-1;
  G4int i, j, k;

  for(i=i1beg; i<=i2end; i++) {
    if (std::abs(r[i]) < spatialTolerance) r[i] = 0.;
  }

  // external polyline - check position of nodes relative to Z
  //
  G4int Nverts = 0;
  for (i=i1beg; i<=i1end; i++) {
    Nverts += (r[i] == 0.) ? 1 : nVphi;
  }

  // internal polyline
  //
  G4bool ifSide1 = false; // whether to create bottom faces
  G4bool ifSide2 = false; // whether to create top faces

  if (r[i2beg] != r[i1beg] || z[i2beg] != z[i1beg]) { // first node
    Nverts += (r[i2beg] == 0.) ? 1 : nVphi;
    ifSide1 = true;
  }

  for(i=i2beg+1; i<i2end; i++) { // intermediate nodes
    Nverts += (r[i] == 0.) ? 1 : nVphi;
  }

  if (r[i2end] != r[i1end] || z[i2end] != z[i1end]) { // last node
    if (absNp2 > 1) Nverts += (r[i2end] == 0.) ? 1 : nVphi;
    ifSide2 = true;
  }

  //   C O U N T   F A C E S

  // external lateral faces
  //
  G4int Nfaces = ifClosed ? absNp1*nSphi : (absNp1-1)*nSphi;

  // internal lateral faces
  //
  if (absNp2 > 1) {
    for(i=i2beg; i<i2end; i++) {
      if (r[i] > 0. || r[i+1] > 0.) Nfaces += nSphi;
    }

    if (ifClosed) {
      if (r[i2end] > 0. || r[i2beg] > 0.) Nfaces += nSphi;
    }
  }

  // bottom and top faces
  //
  if (!ifClosed) {
    if (ifSide1 && (r[i1beg] > 0. || r[i2beg] > 0.)) Nfaces += nSphi;
    if (ifSide2 && (r[i1end] > 0. || r[i2end] > 0.)) Nfaces += nSphi;
  }

  // phi_wedge faces
  //
  if (!ifWholeCircle) {
    Nfaces += ifClosed ? 2*absNp1 : 2*(absNp1-1);
  }

  //   A L L O C A T E   M E M O R Y

  AllocateMemory(Nverts, Nfaces);
  if (pV == nullptr || pF == nullptr) return;

  //   G E N E R A T E   V E R T I C E S

  G4int *kk; // array of start indices along polylines
  kk = new G4int[absNp1+absNp2];

  // external polyline
  //
  k = 1; // free position in array of vertices pV
  for(i=i1beg; i<=i1end; i++) {
    kk[i] = k;
    if (r[i] == 0.)
    { pV[k++] = G4Point3D(0, 0, z[i]); } else { k += nVphi; }
  }

  // first point of internal polyline
  //
  i = i2beg;
  if (ifSide1) {
    kk[i] = k;
    if (r[i] == 0.)
    { pV[k++] = G4Point3D(0, 0, z[i]); } else { k += nVphi; }
  }else{
    kk[i] = kk[i1beg];
  }

  // intermediate points of internal polyline
  //
  for(i=i2beg+1; i<i2end; i++) {
    kk[i] = k;
    if (r[i] == 0.)
    { pV[k++] = G4Point3D(0, 0, z[i]); } else { k += nVphi; }
  }

  // last point of internal polyline
  //
  if (absNp2 > 1) {
    i = i2end;
    if (ifSide2) {
      kk[i] = k;
      if (r[i] == 0.) pV[k] = G4Point3D(0, 0, z[i]);
    }else{
      kk[i] = kk[i1end];
    }
  }

  // set vertices
  //
  G4double cosPhi, sinPhi;

  for(j=0; j<nVphi; j++) {
    cosPhi = std::cos(phi+j*delPhi/nSphi);
    sinPhi = std::sin(phi+j*delPhi/nSphi);
    for(i=i1beg; i<=i2end; i++) {
      if (r[i] != 0.)
        pV[kk[i]+j] = G4Point3D(r[i]*cosPhi,r[i]*sinPhi,z[i]);
    }
  }

  //   G E N E R A T E   F A C E S

  //  external faces
  //
  G4int v1,v2;

  k = 1; // free position in array of faces pF
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

  // internal faces
  //
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

  // bottom and top faces
  //
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

  // phi_wedge faces in case of incomplete circle
  //
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
        SetSideFacets(ii, vv, kk, r, delPhi, nSphi, k);
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
        SetSideFacets(ii, vv, kk, r, delPhi, nSphi, k);
      }
    }
  }

  delete [] kk; // free memory

  // final check
  //
  if (k-1 != nface) {
    std::cerr
      << "HepPolyhedron::RotateAroundZ: number of generated faces ("
      << k-1 << ") is not equal to the number of allocated faces ("
      << nface << ")"
      << std::endl;
  }
}

void
HepPolyhedron::RotateContourAroundZ(G4int nstep,
                                    G4double phi,
                                    G4double dphi,
                                    const std::vector<G4TwoVector> &rz,
                                    G4int nodeVis,
                                    G4int edgeVis)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::RotateContourAroundZ         Date:    12.05.21 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Create HepPolyhedron for a solid produced by rotation of  *
 *           a closed polyline (rz-contour) around Z-axis              *
 *                                                                     *
 * Input: nstep - number of discrete steps, if 0 then default          *
 *        phi   - starting phi angle                                   *
 *        dphi  - delta phi                                            *
 *        rz    - rz-contour                                           *
 *        nodeVis - how to Draw edges joing consecutive positions of   *
 *                  node during rotation                               *
 *        edgeVis - how to Draw edges                                  *
 *                                                                     *
 ***********************************************************************/
{
  //   S E T   R O T A T I O N   P A R A M E T E R S

  G4bool ifWholeCircle = std::abs(dphi - twopi) < perMillion;
  G4double delPhi = (ifWholeCircle) ? twopi : dphi;
  G4int nSphi = nstep;
  if (nSphi <= 0) nSphi = GetNumberOfRotationSteps()*delPhi/twopi + 0.5;
  if (nSphi == 0) nSphi = 1;
  G4int nVphi = (ifWholeCircle) ? nSphi : nSphi + 1;

  //   C A L C U L A T E   A R E A

  G4int Nrz = (G4int)rz.size();
  G4double area = 0;
  for (G4int i = 0; i < Nrz; ++i)
  {
    G4int k = (i == 0) ? Nrz - 1 : i - 1;
    area += rz[k].x()*rz[i].y() - rz[i].x()*rz[k].y();
  }

  //   P R E P A R E   P O L Y L I N E

  auto r = new G4double[Nrz];
  auto z = new G4double[Nrz];
  for (G4int i = 0; i < Nrz; ++i)
  {
    r[i] = rz[i].x();
    z[i] = rz[i].y();
    if (std::abs(r[i]) < spatialTolerance) r[i] = 0.;
  }

  //   C O U N T   V E R T I C E S   A N D   F A C E S

  G4int Nverts = 0;
  for(G4int i = 0; i < Nrz; ++i) Nverts += (r[i] == 0.) ? 1 : nVphi;

  G4int Nedges = Nrz;
  for (G4int i = 0; i < Nrz; ++i)
  {
    G4int k = (i == 0) ? Nrz - 1 : i - 1;
    Nedges -= static_cast<int>(r[k] == 0 && r[i] == 0);
  }

  G4int Nfaces = Nedges*nSphi;               // lateral faces
  if (!ifWholeCircle) Nfaces += 2*(Nrz - 2); // phi_wedge faces

  //   A L L O C A T E   M E M O R Y

  AllocateMemory(Nverts, Nfaces);
  if (pV == nullptr || pF == nullptr)
  {
    delete [] r;
    delete [] z;
    return;
  }

  //   S E T   V E R T I C E S

  auto kk = new G4int[Nrz]; // start indices along contour
  G4int kfree = 1; // current free position in array of vertices pV

  // set start indices, set vertices for nodes with r == 0
  for(G4int i = 0; i < Nrz; ++i)
  {
    kk[i] = kfree;
    if (r[i] == 0.) pV[kfree++] = G4Point3D(0, 0, z[i]);
    if (r[i] != 0.) kfree += nVphi;
  }

  // set vertices by rotating r
  for(G4int j = 0; j < nVphi; ++j)
  {
    G4double cosPhi = std::cos(phi + j*delPhi/nSphi);
    G4double sinPhi = std::sin(phi + j*delPhi/nSphi);
    for(G4int i = 0; i < Nrz; ++i)
    {
      if (r[i] != 0.)
        pV[kk[i] + j] = G4Point3D(r[i]*cosPhi, r[i]*sinPhi, z[i]);
    }
  }

  //   S E T   F A C E S

  kfree = 1; // current free position in array of faces pF
  for(G4int i = 0; i < Nrz; ++i)
  {
    G4int i1 = (i < Nrz - 1) ? i + 1 : 0; // inverse order if area > 0
    G4int i2 = i;
    if (area < 0.) std::swap(i1, i2);
    RotateEdge(kk[i1], kk[i2], r[i1], r[i2], nodeVis, nodeVis,
               edgeVis, ifWholeCircle, nSphi, kfree);
  }

  //    S E T   P H I _ W E D G E   F A C E S

  if (!ifWholeCircle)
  {
    std::vector<G4int> triangles;
    TriangulatePolygon(rz, triangles);

    G4int ii[4], vv[4];
    G4int ntria = G4int(triangles.size()/3);
    for (G4int i = 0; i < ntria; ++i)
    {
      G4int i1 = triangles[0 + i*3];
      G4int i2 = triangles[1 + i*3];
      G4int i3 = triangles[2 + i*3];
      if (area < 0.) std::swap(i1, i3);
      G4int v1 = (std::abs(i2-i1) == 1 || std::abs(i2-i1) == Nrz-1) ? 1 : -1;
      G4int v2 = (std::abs(i3-i2) == 1 || std::abs(i3-i2) == Nrz-1) ? 1 : -1;
      G4int v3 = (std::abs(i1-i3) == 1 || std::abs(i1-i3) == Nrz-1) ? 1 : -1;
      ii[0] = i1; ii[1] = i2; ii[2] = i2; ii[3] = i3;
      vv[0] = v1; vv[1] = -1; vv[2] = v2; vv[3] = v3;
      SetSideFacets(ii, vv, kk, r, delPhi, nSphi, kfree);
    }
  }

  // free memory
  delete [] r;
  delete [] z;
  delete [] kk;

  // final check
  if (kfree - 1 != nface)
  {
    std::cerr
      << "HepPolyhedron::RotateContourAroundZ: number of generated faces ("
      << kfree-1 << ") is not equal to the number of allocated faces ("
      << nface << ")"
      << std::endl;
  }
}

G4bool
HepPolyhedron::TriangulatePolygon(const std::vector<G4TwoVector> &polygon,
                                  std::vector<G4int> &result)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::TriangulatePolygon           Date:    12.05.21 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Simple implementation of "ear clipping" algorithm for     *
 *           triangulation of a simple contour/polygon, it places      *
 *           the result in a std::vector as triplets of vertex indices *
 *                                                                     *
 *           If triangulation is sucsessfull then the function         *
 *           returns true, otherwise false                             *
 *                                                                     *
 * Remark:   It's a copy of G4GeomTools::TriangulatePolygon()          *
 *                                                                     *
 ***********************************************************************/
{
  result.resize(0);
  G4int n = (G4int)polygon.size();
  if (n < 3) return false;

  // calculate area
  //
  G4double area = 0.;
  for(G4int i = 0; i < n; ++i)
  {
    G4int k = (i == 0) ? n - 1 : i - 1;
    area += polygon[k].x()*polygon[i].y() - polygon[i].x()*polygon[k].y();
  }

  // allocate and initialize list of Vertices
  // we want a counter-clockwise polygon in V
  //
  auto  V = new G4int[n];
  if (area > 0.)
    for (G4int i = 0; i < n; ++i) V[i] = i;
  else
    for (G4int i = 0; i < n; ++i) V[i] = (n - 1) - i;

  //  Triangulation: remove nv-2 Vertices, creating 1 triangle every time
  //
  G4int nv = n;
  G4int count = 2*nv; // error detection counter
  for(G4int b = nv - 1; nv > 2; )
  {
    // ERROR: if we loop, it is probably a non-simple polygon
    if ((count--) <= 0)
    {
      delete [] V;
      if (area < 0.) std::reverse(result.begin(),result.end());
      return false;
    }

    // three consecutive vertices in current polygon, <a,b,c>
    G4int a = (b   < nv) ? b   : 0; // previous
          b = (a+1 < nv) ? a+1 : 0; // current
    G4int c = (b+1 < nv) ? b+1 : 0; // next

    if (CheckSnip(polygon, a,b,c, nv,V))
    {
      // output Triangle
      result.push_back(V[a]);
      result.push_back(V[b]);
      result.push_back(V[c]);

      // remove vertex b from remaining polygon
      nv--;
      for(G4int i = b; i < nv; ++i) V[i] = V[i+1];

      count = 2*nv; // resest error detection counter
    }
  }
  delete [] V;
  if (area < 0.) std::reverse(result.begin(),result.end());
  return true;
}

G4bool HepPolyhedron::CheckSnip(const std::vector<G4TwoVector> &contour,
                                G4int a, G4int b, G4int c,
                                G4int n, const G4int* V)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::CheckSnip                    Date:    12.05.21 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Check for a valid snip,                                   *
 *           it is a helper functionfor TriangulatePolygon()           *
 *                                                                     *
 ***********************************************************************/
{
  static const G4double kCarTolerance = 1.e-9;

  // check orientation of Triangle
  G4double Ax = contour[V[a]].x(), Ay = contour[V[a]].y();
  G4double Bx = contour[V[b]].x(), By = contour[V[b]].y();
  G4double Cx = contour[V[c]].x(), Cy = contour[V[c]].y();
  if ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax) < kCarTolerance) return false;

  // check that there is no point inside Triangle
  G4double xmin = std::min(std::min(Ax,Bx),Cx);
  G4double xmax = std::max(std::max(Ax,Bx),Cx);
  G4double ymin = std::min(std::min(Ay,By),Cy);
  G4double ymax = std::max(std::max(Ay,By),Cy);

  for (G4int i=0; i<n; ++i)
  {
    if((i == a) || (i == b) || (i == c)) continue;
    G4double Px = contour[V[i]].x();
    if (Px < xmin || Px > xmax) continue;
    G4double Py = contour[V[i]].y();
    if (Py < ymin || Py > ymax) continue;
    // if (PointInTriangle(Ax,Ay,Bx,By,Cx,Cy,Px,Py)) return false;
    if ((Bx-Ax)*(Cy-Ay) - (By-Ay)*(Cx-Ax) > 0.)
    {
      if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) < 0.) continue;
      if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) < 0.) continue;
      if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) < 0.) continue;
    }
    else
    {
      if ((Ax-Cx)*(Py-Cy) - (Ay-Cy)*(Px-Cx) > 0.) continue;
      if ((Bx-Ax)*(Py-Ay) - (By-Ay)*(Px-Ax) > 0.) continue;
      if ((Cx-Bx)*(Py-By) - (Cy-By)*(Px-Bx) > 0.) continue;
    }
    return false;
  }
  return true;
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
    G4int v2;
    G4int iface;
    G4int iedge;
  } *edgeList, *freeList, **headList;


  //   A L L O C A T E   A N D   I N I T I A T E   L I S T S

  edgeList = new edgeListMember[2*nface];
  headList = new edgeListMember*[nvert];

  G4int i;
  for (i=0; i<nvert; i++) {
    headList[i] = nullptr;
  }
  freeList = edgeList;
  for (i=0; i<2*nface-1; i++) {
    edgeList[i].next = &edgeList[i+1];
  }
  edgeList[2*nface-1].next = nullptr;

  //   L O O P   A L O N G   E D G E S

  G4int iface, iedge, nedge, i1, i2, k1, k2;
  edgeListMember *prev, *cur;

  for(iface=1; iface<=nface; iface++) {
    nedge = (pF[iface].edge[3].v == 0) ? 3 : 4;
    for (iedge=0; iedge<nedge; iedge++) {
      i1 = iedge;
      i2 = (iedge < nedge-1) ? iedge+1 : 0;
      i1 = std::abs(pF[iface].edge[i1].v);
      i2 = std::abs(pF[iface].edge[i2].v);
      k1 = (i1 < i2) ? i1 : i2;          // k1 = ::min(i1,i2);
      k2 = (i1 > i2) ? i1 : i2;          // k2 = ::max(i1,i2);

      // check head of the List corresponding to k1
      cur = headList[k1];
      if (cur == nullptr) {
        headList[k1] = freeList;
        if (freeList == nullptr) {
          std::cerr
          << "Polyhedron::SetReferences: bad link "
          << std::endl;
          break;
        }
        freeList = freeList->next;
        cur = headList[k1];
        cur->next = nullptr;
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
        if (cur == nullptr) {
          prev->next = freeList;
          if (freeList == nullptr) {
            std::cerr
            << "Polyhedron::SetReferences: bad link "
            << std::endl;
            break;
          }
          freeList = freeList->next;
          cur = prev->next;
          cur->next = nullptr;
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
    if (headList[i] != nullptr) {
      std::cerr
        << "Polyhedron::SetReferences: List " << i << " is not empty"
        << std::endl;
    }
  }

  //   F R E E   M E M O R Y

  delete [] edgeList;
  delete [] headList;
}

void HepPolyhedron::JoinCoplanarFacets(G4double tolerance)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::JoinCoplanarFacets          Date:    04.04.22  *
 * Author: E.Tcherniaev (E.Chernyaev)               Revised:           *
 *                                                                     *
 * Function: Join couples of triangular facets to quadrangular facets  *
 *           where it is possible                                      *
 *                                                                     *
 ***********************************************************************/
{
  G4int njoin = 0;
  for (G4int icur = 1; icur <= nface; ++icur)
  {
    // skip if already joined or quadrangle
    if (pF[icur].edge[0].v == 0) continue;
    if (pF[icur].edge[3].v != 0) continue;
    // skip if all references point to already checked facets
    if (pF[icur].edge[0].f < icur &&
        pF[icur].edge[1].f < icur &&
        pF[icur].edge[2].f < icur) continue;
    // compute plane equation
    G4Normal3D norm = GetUnitNormal(icur);
    G4double dd = norm.dot(pV[pF[icur].edge[0].v]);
    G4int vcur0 = std::abs(pF[icur].edge[0].v);
    G4int vcur1 = std::abs(pF[icur].edge[1].v);
    G4int vcur2 = std::abs(pF[icur].edge[2].v);
    // select neighbouring facet
    G4int kcheck = 0, icheck = 0, vcheck = 0;
    G4double dist = DBL_MAX;
    for (G4int k = 0; k < 3; ++k)
    {
      G4int itmp = pF[icur].edge[k].f;
      // skip if already checked, joined or quadrangle
      if (itmp < icur) continue;
      if (pF[itmp].edge[0].v == 0 ||
          pF[itmp].edge[3].v != 0) continue;
      // get candidate vertex
      G4int vtmp = 0;
      for (G4int j = 0; j < 3; ++j)
      {
        vtmp = std::abs(pF[itmp].edge[j].v);
	if (vtmp != vcur0 && vtmp != vcur1 && vtmp != vcur2) break;
      }
      // check distance to the plane
      G4double dtmp = std::abs(norm.dot(pV[vtmp]) - dd);
      if (dtmp > tolerance || dtmp >= dist) continue;
      dist = dtmp;
      kcheck = k;
      icheck = itmp;
      vcheck = vtmp;
    }
    if (icheck == 0) continue; // no facet selected
    // join facets
    njoin++;
    pF[icheck].edge[0].v = 0; // mark facet as joined
    if (kcheck == 0)
    {
      pF[icur].edge[3].v = pF[icur].edge[2].v;
      pF[icur].edge[2].v = pF[icur].edge[1].v;
      pF[icur].edge[1].v = vcheck;
    }
    else if (kcheck == 1)
    {
      pF[icur].edge[3].v = pF[icur].edge[2].v;
      pF[icur].edge[2].v = vcheck;
    }
    else
    {
      pF[icur].edge[3].v = vcheck;
    }
  }
  if (njoin == 0) return; // no joined facets

  // restructure facets
  G4int nnew = 0;
  for (G4int icur = 1; icur <= nface; ++icur)
  {
    if (pF[icur].edge[0].v == 0) continue;
    nnew++;
    pF[nnew].edge[0].v = pF[icur].edge[0].v;
    pF[nnew].edge[1].v = pF[icur].edge[1].v;
    pF[nnew].edge[2].v = pF[icur].edge[2].v;
    pF[nnew].edge[3].v = pF[icur].edge[3].v;
  }
  nface = nnew;
  SetReferences();
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
  G4int i, k, nnode, v[4],f[4];
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

HepPolyhedron & HepPolyhedron::Transform(const G4Transform3D &t)
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
    for (G4int i=1; i<=nvert; i++) { pV[i] = t * pV[i]; }

    //  C H E C K   D E T E R M I N A N T   A N D
    //  I N V E R T   F A C E T S   I F   I T   I S   N E G A T I V E

    G4Vector3D d = t * G4Vector3D(0,0,0);
    G4Vector3D x = t * G4Vector3D(1,0,0) - d;
    G4Vector3D y = t * G4Vector3D(0,1,0) - d;
    G4Vector3D z = t * G4Vector3D(0,0,1) - d;
    if ((x.cross(y))*z < 0) InvertFacets();
  }
  return *this;
}

G4bool HepPolyhedron::GetNextVertexIndex(G4int &index, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextVertexIndex          Date:    03.09.96  *
 * Author: Yasuhide Sawada                          Revised:           *
 *                                                                     *
 * Function:                                                           *
 *                                                                     *
 ***********************************************************************/
{
  static G4ThreadLocal G4int iFace = 1;
  static G4ThreadLocal G4int iQVertex = 0;
  G4int vIndex = pF[iFace].edge[iQVertex].v;

  edgeFlag = (vIndex > 0) ? 1 : 0;
  index = std::abs(vIndex);

  if (iQVertex >= 3 || pF[iFace].edge[iQVertex+1].v == 0) {
    iQVertex = 0;
    if (++iFace > nface) iFace = 1;
    return false;  // Last Edge
  }
  
  ++iQVertex;
  return true;  // not Last Edge
}

G4Point3D HepPolyhedron::GetVertex(G4int index) const
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
    return G4Point3D();
  }
  return pV[index];
}

G4bool
HepPolyhedron::GetNextVertex(G4Point3D &vertex, G4int &edgeFlag) const
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
  G4int index;
  G4bool rep = GetNextVertexIndex(index, edgeFlag);
  vertex = pV[index];
  return rep;
}

G4bool HepPolyhedron::GetNextVertex(G4Point3D &vertex, G4int &edgeFlag,
                                  G4Normal3D &normal) const
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
  static G4ThreadLocal G4int iFace = 1;
  static G4ThreadLocal G4int iNode = 0;

  if (nface == 0) return false;  // empty polyhedron

  G4int k = pF[iFace].edge[iNode].v;
  if (k > 0) { edgeFlag = 1; } else { edgeFlag = -1; k = -k; }
  vertex = pV[k];
  normal = FindNodeNormal(iFace,k);
  if (iNode >= 3 || pF[iFace].edge[iNode+1].v == 0) {
    iNode = 0;
    if (++iFace > nface) iFace = 1;
    return false;                // last node
  }
  ++iNode;
  return true;                 // not last node
}

G4bool HepPolyhedron::GetNextEdgeIndices(G4int &i1, G4int &i2, G4int &edgeFlag,
                                       G4int &iface1, G4int &iface2) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdgeIndices          Date:    30.09.96  *
 * Author: E.Chernyaev                              Revised: 17.11.99  *
 *                                                                     *
 * Function: Get indices of the next edge together with indices of     *
 *           of the faces which share the edge.                        *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  static G4ThreadLocal G4int iFace    = 1;
  static G4ThreadLocal G4int iQVertex = 0;
  static G4ThreadLocal G4int iOrder   = 1;
  G4int  k1, k2, kflag, kface1, kface2;

  if (iFace == 1 && iQVertex == 0) {
    k2 = pF[nface].edge[0].v;
    k1 = pF[nface].edge[3].v;
    if (k1 == 0) k1 = pF[nface].edge[2].v;
    if (std::abs(k1) > std::abs(k2)) iOrder = -1;
  }

  do {
    k1     = pF[iFace].edge[iQVertex].v;
    kflag  = k1;
    k1     = std::abs(k1);
    kface1 = iFace;
    kface2 = pF[iFace].edge[iQVertex].f;
    if (iQVertex >= 3 || pF[iFace].edge[iQVertex+1].v == 0) {
      iQVertex = 0;
      k2 = std::abs(pF[iFace].edge[iQVertex].v);
      iFace++;
    }else{
      iQVertex++;
      k2 = std::abs(pF[iFace].edge[iQVertex].v);
    }
  } while (iOrder*k1 > iOrder*k2);

  i1 = k1; i2 = k2; edgeFlag = (kflag > 0) ? 1 : 0;
  iface1 = kface1; iface2 = kface2;

  if (iFace > nface) {
    iFace  = 1; iOrder = 1;
    return false;
  }

  return true;
}

G4bool
HepPolyhedron::GetNextEdgeIndices(G4int &i1, G4int &i2, G4int &edgeFlag) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdgeIndices          Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get indices of the next edge.                             *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  G4int kface1, kface2;
  return GetNextEdgeIndices(i1, i2, edgeFlag, kface1, kface2);
}

G4bool
HepPolyhedron::GetNextEdge(G4Point3D &p1,
                           G4Point3D &p2,
                           G4int &edgeFlag) const
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
  G4int i1,i2;
  G4bool rep = GetNextEdgeIndices(i1,i2,edgeFlag);
  p1 = pV[i1];
  p2 = pV[i2];
  return rep;
}

G4bool
HepPolyhedron::GetNextEdge(G4Point3D &p1, G4Point3D &p2,
                          G4int &edgeFlag, G4int &iface1, G4int &iface2) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetNextEdge                 Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get next edge with indices of the faces which share       *
 *           the edge.                                                 *
 *           Returns false when the last edge.                         *
 *                                                                     *
 ***********************************************************************/
{
  G4int i1,i2;
  G4bool rep = GetNextEdgeIndices(i1,i2,edgeFlag,iface1,iface2);
  p1 = pV[i1];
  p2 = pV[i2];
  return rep;
}

void HepPolyhedron::GetFacet(G4int iFace, G4int &n, G4int *iNodes,
                            G4int *edgeFlags, G4int *iFaces) const
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
    G4int i, k;
    for (i=0; i<4; i++) {
      k = pF[iFace].edge[i].v;
      if (k == 0) break;
      if (iFaces != nullptr) iFaces[i] = pF[iFace].edge[i].f;
      if (k > 0) {
        iNodes[i] = k;
        if (edgeFlags != nullptr) edgeFlags[i] = 1;
      }else{
        iNodes[i] = -k;
        if (edgeFlags != nullptr) edgeFlags[i] = -1;
      }
    }
    n = i;
  }
}

void HepPolyhedron::GetFacet(G4int index, G4int &n, G4Point3D *nodes,
                             G4int *edgeFlags, G4Normal3D *normals) const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetFacet                    Date:    17.11.99  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Get face by index                                         *
 *                                                                     *
 ***********************************************************************/
{
  G4int iNodes[4];
  GetFacet(index, n, iNodes, edgeFlags);
  if (n != 0) {
    for (G4int i=0; i<n; i++) {
      nodes[i] = pV[iNodes[i]];
      if (normals != nullptr) normals[i] = FindNodeNormal(index,iNodes[i]);
    }
  }
}

G4bool
HepPolyhedron::GetNextFacet(G4int &n, G4Point3D *nodes,
                           G4int *edgeFlags, G4Normal3D *normals) const
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
  static G4ThreadLocal G4int iFace = 1;

  if (edgeFlags == nullptr) {
    GetFacet(iFace, n, nodes);
  }else if (normals == nullptr) {
    GetFacet(iFace, n, nodes, edgeFlags);
  }else{
    GetFacet(iFace, n, nodes, edgeFlags, normals);
  }

  if (++iFace > nface) {
    iFace  = 1;
    return false;
  }

  return true;
}

G4Normal3D HepPolyhedron::GetNormal(G4int iFace) const
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
    return G4Normal3D();
  }

  G4int i0  = std::abs(pF[iFace].edge[0].v);
  G4int i1  = std::abs(pF[iFace].edge[1].v);
  G4int i2  = std::abs(pF[iFace].edge[2].v);
  G4int i3  = std::abs(pF[iFace].edge[3].v);
  if (i3 == 0) i3 = i0;
  return (pV[i2] - pV[i0]).cross(pV[i3] - pV[i1]);
}

G4Normal3D HepPolyhedron::GetUnitNormal(G4int iFace) const
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
    return G4Normal3D();
  }

  G4int i0  = std::abs(pF[iFace].edge[0].v);
  G4int i1  = std::abs(pF[iFace].edge[1].v);
  G4int i2  = std::abs(pF[iFace].edge[2].v);
  G4int i3  = std::abs(pF[iFace].edge[3].v);
  if (i3 == 0) i3 = i0;
  return ((pV[i2] - pV[i0]).cross(pV[i3] - pV[i1])).unit();
}

G4bool HepPolyhedron::GetNextNormal(G4Normal3D &normal) const
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
  static G4ThreadLocal G4int iFace = 1;
  normal = GetNormal(iFace);
  if (++iFace > nface) {
    iFace = 1;
    return false;
  }
  return true;
}

G4bool HepPolyhedron::GetNextUnitNormal(G4Normal3D &normal) const
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
  G4bool rep = GetNextNormal(normal);
  normal = normal.unit();
  return rep;
}

G4double HepPolyhedron::GetSurfaceArea() const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetSurfaceArea              Date:    25.05.01  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Returns area of the surface of the polyhedron.            *
 *                                                                     *
 ***********************************************************************/
{
  G4double srf = 0.;
  for (G4int iFace=1; iFace<=nface; iFace++) {
    G4int i0 = std::abs(pF[iFace].edge[0].v);
    G4int i1 = std::abs(pF[iFace].edge[1].v);
    G4int i2 = std::abs(pF[iFace].edge[2].v);
    G4int i3 = std::abs(pF[iFace].edge[3].v);
    if (i3 == 0) i3 = i0;
    srf += ((pV[i2] - pV[i0]).cross(pV[i3] - pV[i1])).mag();
  }
  return srf/2.;
}

G4double HepPolyhedron::GetVolume() const
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::GetVolume                   Date:    25.05.01  *
 * Author: E.Chernyaev                              Revised:           *
 *                                                                     *
 * Function: Returns volume of the polyhedron.                         *
 *                                                                     *
 ***********************************************************************/
{
  G4double v = 0.;
  for (G4int iFace=1; iFace<=nface; iFace++) {
    G4int i0 = std::abs(pF[iFace].edge[0].v);
    G4int i1 = std::abs(pF[iFace].edge[1].v);
    G4int i2 = std::abs(pF[iFace].edge[2].v);
    G4int i3 = std::abs(pF[iFace].edge[3].v);
    G4Point3D pt;
    if (i3 == 0) {
      i3 = i0;
      pt = (pV[i0]+pV[i1]+pV[i2]) * (1./3.);
    }else{
      pt = (pV[i0]+pV[i1]+pV[i2]+pV[i3]) * 0.25;
    }
    v += ((pV[i2] - pV[i0]).cross(pV[i3] - pV[i1])).dot(pt);
  }
  return v/6.;
}

G4int
HepPolyhedron::createTwistedTrap(G4double Dz,
                                 const G4double xy1[][2],
                                 const G4double xy2[][2])
/***********************************************************************
 *                                                                     *
 * Name: createTwistedTrap                           Date:    05.11.02 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Creates polyhedron for twisted trapezoid                  *
 *                                                                     *
 * Input: Dz       - half-length along Z             8----7            *
 *        xy1[2,4] - quadrilateral at Z=-Dz       5----6  !            *
 *        xy2[2,4] - quadrilateral at Z=+Dz       !  4-!--3            *
 *                                                1----2               *
 *                                                                     *
 ***********************************************************************/
{
  AllocateMemory(12,18);

  pV[ 1] = G4Point3D(xy1[0][0],xy1[0][1],-Dz);
  pV[ 2] = G4Point3D(xy1[1][0],xy1[1][1],-Dz);
  pV[ 3] = G4Point3D(xy1[2][0],xy1[2][1],-Dz);
  pV[ 4] = G4Point3D(xy1[3][0],xy1[3][1],-Dz);

  pV[ 5] = G4Point3D(xy2[0][0],xy2[0][1], Dz);
  pV[ 6] = G4Point3D(xy2[1][0],xy2[1][1], Dz);
  pV[ 7] = G4Point3D(xy2[2][0],xy2[2][1], Dz);
  pV[ 8] = G4Point3D(xy2[3][0],xy2[3][1], Dz);

  pV[ 9] = (pV[1]+pV[2]+pV[5]+pV[6])/4.;
  pV[10] = (pV[2]+pV[3]+pV[6]+pV[7])/4.;
  pV[11] = (pV[3]+pV[4]+pV[7]+pV[8])/4.;
  pV[12] = (pV[4]+pV[1]+pV[8]+pV[5])/4.;

  enum {DUMMY, BOTTOM,
        LEFT_BOTTOM,  LEFT_FRONT,   LEFT_TOP,  LEFT_BACK,
        BACK_BOTTOM,  BACK_LEFT,    BACK_TOP,  BACK_RIGHT,
        RIGHT_BOTTOM, RIGHT_BACK,   RIGHT_TOP, RIGHT_FRONT,
        FRONT_BOTTOM, FRONT_RIGHT,  FRONT_TOP, FRONT_LEFT,
        TOP};

  pF[ 1]=G4Facet(1,LEFT_BOTTOM, 4,BACK_BOTTOM, 3,RIGHT_BOTTOM, 2,FRONT_BOTTOM);

  pF[ 2]=G4Facet(4,BOTTOM,     -1,LEFT_FRONT,  -12,LEFT_BACK,    0,0);
  pF[ 3]=G4Facet(1,FRONT_LEFT, -5,LEFT_TOP,    -12,LEFT_BOTTOM,  0,0);
  pF[ 4]=G4Facet(5,TOP,        -8,LEFT_BACK,   -12,LEFT_FRONT,   0,0);
  pF[ 5]=G4Facet(8,BACK_LEFT,  -4,LEFT_BOTTOM, -12,LEFT_TOP,     0,0);

  pF[ 6]=G4Facet(3,BOTTOM,     -4,BACK_LEFT,   -11,BACK_RIGHT,   0,0);
  pF[ 7]=G4Facet(4,LEFT_BACK,  -8,BACK_TOP,    -11,BACK_BOTTOM,  0,0);
  pF[ 8]=G4Facet(8,TOP,        -7,BACK_RIGHT,  -11,BACK_LEFT,    0,0);
  pF[ 9]=G4Facet(7,RIGHT_BACK, -3,BACK_BOTTOM, -11,BACK_TOP,     0,0);

  pF[10]=G4Facet(2,BOTTOM,     -3,RIGHT_BACK,  -10,RIGHT_FRONT,  0,0);
  pF[11]=G4Facet(3,BACK_RIGHT, -7,RIGHT_TOP,   -10,RIGHT_BOTTOM, 0,0);
  pF[12]=G4Facet(7,TOP,        -6,RIGHT_FRONT, -10,RIGHT_BACK,   0,0);
  pF[13]=G4Facet(6,FRONT_RIGHT,-2,RIGHT_BOTTOM,-10,RIGHT_TOP,    0,0);

  pF[14]=G4Facet(1,BOTTOM,     -2,FRONT_RIGHT,  -9,FRONT_LEFT,   0,0);
  pF[15]=G4Facet(2,RIGHT_FRONT,-6,FRONT_TOP,    -9,FRONT_BOTTOM, 0,0);
  pF[16]=G4Facet(6,TOP,        -5,FRONT_LEFT,   -9,FRONT_RIGHT,  0,0);
  pF[17]=G4Facet(5,LEFT_FRONT, -1,FRONT_BOTTOM, -9,FRONT_TOP,    0,0);

  pF[18]=G4Facet(5,FRONT_TOP, 6,RIGHT_TOP, 7,BACK_TOP, 8,LEFT_TOP);

  return 0;
}

G4int
HepPolyhedron::createPolyhedron(G4int Nnodes, G4int Nfaces,
                                const G4double xyz[][3],
                                const G4int  faces[][4])
/***********************************************************************
 *                                                                     *
 * Name: createPolyhedron                            Date:    05.11.02 *
 * Author: E.Chernyaev (IHEP/Protvino)               Revised:          *
 *                                                                     *
 * Function: Creates user defined polyhedron                           *
 *                                                                     *
 * Input: Nnodes  - number of nodes                                    *
 *        Nfaces  - number of faces                                    *
 *        nodes[][3] - node coordinates                                *
 *        faces[][4] - faces                                           *
 *                                                                     *
 ***********************************************************************/
{
  AllocateMemory(Nnodes, Nfaces);
  if (nvert == 0) return 1;

  for (G4int i=0; i<Nnodes; i++) {
    pV[i+1] = G4Point3D(xyz[i][0], xyz[i][1], xyz[i][2]);
  }
  for (G4int k=0; k<Nfaces; k++) {
    pF[k+1] = G4Facet(faces[k][0],0,faces[k][1],0,faces[k][2],0,faces[k][3],0);
  }
  SetReferences();
  return 0;
}

HepPolyhedronTrd2::HepPolyhedronTrd2(G4double Dx1, G4double Dx2,
                                     G4double Dy1, G4double Dy2,
                                     G4double Dz)
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

HepPolyhedronTrd2::~HepPolyhedronTrd2() = default;

HepPolyhedronTrd1::HepPolyhedronTrd1(G4double Dx1, G4double Dx2,
                                     G4double Dy, G4double Dz)
  : HepPolyhedronTrd2(Dx1, Dx2, Dy, Dy, Dz) {}

HepPolyhedronTrd1::~HepPolyhedronTrd1() = default;

HepPolyhedronBox::HepPolyhedronBox(G4double Dx, G4double Dy, G4double Dz)
  : HepPolyhedronTrd2(Dx, Dx, Dy, Dy, Dz) {}

HepPolyhedronBox::~HepPolyhedronBox() = default;

HepPolyhedronTrap::HepPolyhedronTrap(G4double Dz,
                                     G4double Theta,
                                     G4double Phi,
                                     G4double Dy1,
                                     G4double Dx1,
                                     G4double Dx2,
                                     G4double Alp1,
                                     G4double Dy2,
                                     G4double Dx3,
                                     G4double Dx4,
                                     G4double Alp2)
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
  G4double DzTthetaCphi = Dz*std::tan(Theta)*std::cos(Phi);
  G4double DzTthetaSphi = Dz*std::tan(Theta)*std::sin(Phi);
  G4double Dy1Talp1 = Dy1*std::tan(Alp1);
  G4double Dy2Talp2 = Dy2*std::tan(Alp2);

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

HepPolyhedronTrap::~HepPolyhedronTrap() = default;

HepPolyhedronPara::HepPolyhedronPara(G4double Dx, G4double Dy, G4double Dz,
                                     G4double Alpha, G4double Theta,
                                     G4double Phi)
  : HepPolyhedronTrap(Dz, Theta, Phi, Dy, Dx, Dx, Alpha, Dy, Dx, Dx, Alpha) {}

HepPolyhedronPara::~HepPolyhedronPara() = default;

HepPolyhedronParaboloid::HepPolyhedronParaboloid(G4double r1,
                                                 G4double r2,
                                                 G4double dz,
                                                 G4double sPhi,
                                                 G4double dPhi)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronParaboloid                     Date:    28.06.07 *
 * Author: L.Lindroos, T.Nikitina (CERN), July 2007  Revised: 28.06.07 *
 *                                                                     *
 * Function: Constructor for paraboloid                                *
 *                                                                     *
 * Input: r1    - inside and outside radiuses at -Dz                   *
 *        r2    - inside and outside radiuses at +Dz                   *
 *        dz    - half length in Z                                     *
 *        sPhi  - starting angle of the segment                        *
 *        dPhi  - segment range                                        *
 *                                                                     *
 ***********************************************************************/
{
  static const G4double wholeCircle=twopi;

  //   C H E C K   I N P U T   P A R A M E T E R S

  G4int k = 0;
  if (r1 < 0. || r2 <= 0.)        k = 1;

  if (dz <= 0.) k += 2;

  G4double phi1, phi2, dphi;

  if(dPhi < 0.)
  {
    phi2 = sPhi; phi1 = phi2 + dPhi;
  }
  else if(dPhi == 0.)
  {
    phi1 = sPhi; phi2 = phi1 + wholeCircle;
  }
  else
  {
    phi1 = sPhi; phi2 = phi1 + dPhi;
  }
  dphi  = phi2 - phi1;

  if (std::abs(dphi-wholeCircle) < perMillion) dphi = wholeCircle;
  if (dphi > wholeCircle) k += 4;

  if (k != 0) {
    std::cerr << "HepPolyhedronParaboloid: error in input parameters";
    if ((k & 1) != 0) std::cerr << " (radiuses)";
    if ((k & 2) != 0) std::cerr << " (half-length)";
    if ((k & 4) != 0) std::cerr << " (angles)";
    std::cerr << std::endl;
    std::cerr << " r1=" << r1;
    std::cerr << " r2=" << r2;
    std::cerr << " dz=" << dz << " sPhi=" << sPhi << " dPhi=" << dPhi
              << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4int n = GetNumberOfRotationSteps();
  G4double dl = (r2 - r1) / n;
  G4double k1 = (r2*r2 - r1*r1) / 2 / dz;
  G4double k2 = (r2*r2 + r1*r1) / 2;

  auto zz = new G4double[n + 2], rr = new G4double[n + 2];

  zz[0] = dz;
  rr[0] = r2;

  for(G4int i = 1; i < n - 1; i++)
  {
    rr[i] = rr[i-1] - dl;
    zz[i] = (rr[i]*rr[i] - k2) / k1;
    if(rr[i] < 0)
    {
      rr[i] = 0;
      zz[i] = 0;
    }
  }

  zz[n-1] = -dz;
  rr[n-1] = r1;

  zz[n] = dz;
  rr[n] = 0;

  zz[n+1] = -dz;
  rr[n+1] = 0;

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, phi1, dphi, n, 2, zz, rr, -1, -1);
  SetReferences();

  delete [] zz;
  delete [] rr;
}

HepPolyhedronParaboloid::~HepPolyhedronParaboloid() = default;

HepPolyhedronHype::HepPolyhedronHype(G4double r1,
                                     G4double r2,
                                     G4double sqrtan1,
                                     G4double sqrtan2,
                                     G4double halfZ)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronHype                           Date:    14.04.08 *
 * Author: Tatiana Nikitina (CERN)                   Revised: 14.04.08 *
 *         Evgueni Tcherniaev                                 01.12.20 *
 *                                                                     *
 * Function: Constructor for Hype                                      *
 *                                                                     *
 * Input: r1       - inside radius at z=0                              *
 *        r2       - outside radiuses at z=0                           *
 *        sqrtan1  - sqr of tan of Inner Stereo Angle                  *
 *        sqrtan2  - sqr of tan of Outer Stereo Angle                  *
 *        halfZ    - half length in Z                                  *
 *                                                                     *
 ***********************************************************************/
{
  static const G4double wholeCircle = twopi;

  //   C H E C K   I N P U T   P A R A M E T E R S

  G4int k = 0;
  if (r1 < 0. || r2 < 0. || r1 >= r2) k = 1;
  if (halfZ <= 0.) k += 2;
  if (sqrtan1 < 0.|| sqrtan2 < 0.) k += 4;

  if (k != 0)
  {
    std::cerr << "HepPolyhedronHype: error in input parameters";
    if ((k & 1) != 0) std::cerr << " (radiuses)";
    if ((k & 2) != 0) std::cerr << " (half-length)";
    if ((k & 4) != 0) std::cerr << " (angles)";
    std::cerr << std::endl;
    std::cerr << " r1=" << r1 << " r2=" << r2;
    std::cerr << " halfZ=" << halfZ << " sqrTan1=" << sqrtan1
              << " sqrTan2=" << sqrtan2
              << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4int ns = std::max(3, GetNumberOfRotationSteps()/4);
  G4int nz1 = (sqrtan1 == 0.) ? 2 : ns + 1;
  G4int nz2 = (sqrtan2 == 0.) ? 2 : ns + 1;
  auto  zz = new G4double[nz1 + nz2];
  auto  rr = new G4double[nz1 + nz2];

  // external polyline
  G4double dz2 = 2.*halfZ/(nz2 - 1);
  for(G4int i = 0; i < nz2; ++i)
  {
    zz[i] = halfZ - dz2*i;
    rr[i] = std::sqrt(sqrtan2*zz[i]*zz[i] + r2*r2);
  }

  // internal polyline
  G4double dz1 = 2.*halfZ/(nz1 - 1);
  for(G4int i = 0; i < nz1; ++i)
  {
    G4int j = nz2 + i;
    zz[j] = halfZ - dz1*i;
    rr[j] = std::sqrt(sqrtan1*zz[j]*zz[j] + r1*r1);
  }

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, 0., wholeCircle, nz2, nz1, zz, rr, -1, -1);
  SetReferences();

  delete [] zz;
  delete [] rr;
}

HepPolyhedronHype::~HepPolyhedronHype() = default;

HepPolyhedronCons::HepPolyhedronCons(G4double Rmn1,
                                     G4double Rmx1,
                                     G4double Rmn2,
                                     G4double Rmx2,
                                     G4double Dz,
                                     G4double Phi1,
                                     G4double Dphi)
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
  static const G4double wholeCircle=twopi;

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
  if (std::abs(dphi-wholeCircle) < perMillion) dphi = wholeCircle;
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

  RotateAroundZ(0, phi1, dphi, 2, 2, zz, rr, -1, -1);
  SetReferences();
}

HepPolyhedronCons::~HepPolyhedronCons() = default;

HepPolyhedronCone::HepPolyhedronCone(G4double Rmn1, G4double Rmx1,
                                     G4double Rmn2, G4double Rmx2,
                                     G4double Dz) :
  HepPolyhedronCons(Rmn1, Rmx1, Rmn2, Rmx2, Dz, 0*deg, 360*deg) {}

HepPolyhedronCone::~HepPolyhedronCone() = default;

HepPolyhedronTubs::HepPolyhedronTubs(G4double Rmin, G4double Rmax,
                                     G4double Dz,
                                     G4double Phi1, G4double Dphi)
  :   HepPolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, Phi1, Dphi) {}

HepPolyhedronTubs::~HepPolyhedronTubs() = default;

HepPolyhedronTube::HepPolyhedronTube (G4double Rmin, G4double Rmax,
                                      G4double Dz)
  : HepPolyhedronCons(Rmin, Rmax, Rmin, Rmax, Dz, 0*deg, 360*deg) {}

HepPolyhedronTube::~HepPolyhedronTube () = default;

HepPolyhedronPgon::HepPolyhedronPgon(G4double phi,
                                     G4double dphi,
                                     G4int npdv,
                                     G4int nz,
                                     const G4double *z,
                                     const G4double *rmin,
                                     const G4double *rmax)
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

  if (dphi <= 0. || dphi > twopi) {
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

  G4int i;
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

  G4int nodeVis = 1;
  G4int edgeVis = (npdv == 0) ? -1 : 1;
  RotateAroundZ(npdv, phi, dphi, nz, nz, zz, rr, nodeVis, edgeVis);
  SetReferences();

  delete [] zz;
  delete [] rr;
}

HepPolyhedronPgon::HepPolyhedronPgon(G4double phi,
                                     G4double dphi,
                                     G4int npdv,
                                     const std::vector<G4TwoVector> &rz)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronPgon                           Date:    12.05.21 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Constructor of polyhedron for PGON, PCON                  *
 *                                                                     *
 * Input: phi  - initial phi                                           *
 *        dphi - delta phi                                             *
 *        npdv - number of steps along phi                             *
 *        rz   - rz-contour                                            *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  if (dphi <= 0. || dphi > twopi) {
    std::cerr
      << "HepPolyhedronPgon/Pcon: wrong delta phi = " << dphi
      << std::endl;
    return;
  }

  if (npdv < 0) {
    std::cerr
      << "HepPolyhedronPgon/Pcon: error in number of phi-steps = " << npdv
      << std::endl;
    return;
  }

  G4int nrz = (G4int)rz.size();
  if (nrz < 3) {
    std::cerr
      << "HepPolyhedronPgon/Pcon: invalid number of nodes in rz-contour = " << nrz
      << std::endl;
    return;
  }

  //   R O T A T E    P O L Y L I N E

  G4int nodeVis = 1;
  G4int edgeVis = (npdv == 0) ? -1 : 1;
  RotateContourAroundZ(npdv, phi, dphi, rz, nodeVis, edgeVis);
  SetReferences();
}

HepPolyhedronPgon::~HepPolyhedronPgon() = default;

HepPolyhedronPcon::HepPolyhedronPcon(G4double phi, G4double dphi, G4int nz,
                                     const G4double *z,
                                     const G4double *rmin,
                                     const G4double *rmax)
  : HepPolyhedronPgon(phi, dphi, 0, nz, z, rmin, rmax) {}

HepPolyhedronPcon::HepPolyhedronPcon(G4double phi, G4double dphi,
                                     const std::vector<G4TwoVector> &rz)
  : HepPolyhedronPgon(phi, dphi, 0, rz) {}

HepPolyhedronPcon::~HepPolyhedronPcon() = default;

HepPolyhedronSphere::HepPolyhedronSphere(G4double rmin, G4double rmax,
                                         G4double phi, G4double dphi,
                                         G4double the, G4double dthe)
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

  if (dphi <= 0. || dphi > twopi) {
    std::cerr
      << "HepPolyhedronSphere: wrong delta phi = " << dphi
      << std::endl;
    return;
  }

  if (the < 0. || the > pi) {
    std::cerr
      << "HepPolyhedronSphere: wrong theta = " << the
      << std::endl;
    return;
  }

  if (dthe <= 0. || dthe > pi) {
    std::cerr
      << "HepPolyhedronSphere: wrong delta theta = " << dthe
      << std::endl;
    return;
  }

  if (the+dthe > pi) {
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

  G4int nds = (GetNumberOfRotationSteps() + 1) / 2;
  G4int np1 = G4int(dthe*nds/pi+.5) + 1;
  if (np1 <= 1) np1 = 2;
  G4int np2 = rmin < spatialTolerance ? 1 : np1;

  G4double *zz, *rr;
  zz = new G4double[np1+np2];
  rr = new G4double[np1+np2];

  G4double a = dthe/(np1-1);
  G4double cosa, sina;
  for (G4int i=0; i<np1; i++) {
    cosa  = std::cos(the+i*a);
    sina  = std::sin(the+i*a);
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

HepPolyhedronSphere::~HepPolyhedronSphere() = default;

HepPolyhedronTorus::HepPolyhedronTorus(G4double rmin,
                                       G4double rmax,
                                       G4double rtor,
                                       G4double phi,
                                       G4double dphi)
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

  if (dphi <= 0. || dphi > twopi) {
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

  G4int np1 = GetNumberOfRotationSteps();
  G4int np2 = rmin < spatialTolerance ? 1 : np1;

  G4double *zz, *rr;
  zz = new G4double[np1+np2];
  rr = new G4double[np1+np2];

  G4double a = twopi/np1;
  G4double cosa, sina;
  for (G4int i=0; i<np1; i++) {
    cosa  = std::cos(i*a);
    sina  = std::sin(i*a);
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

HepPolyhedronTorus::~HepPolyhedronTorus() = default;

HepPolyhedronTet::HepPolyhedronTet(const G4double p0[3],
                                   const G4double p1[3],
                                   const G4double p2[3],
                                   const G4double p3[3])
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronTet                            Date:  21.02.2020 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Constructor of polyhedron for TETrahedron                 *
 *                                                                     *
 * Input: p0,p1,p2,p3 - vertices                                       *
 *                                                                     *
 ***********************************************************************/
{
  AllocateMemory(4,4);

  pV[1].set(p0[0], p0[1], p0[2]);
  pV[2].set(p1[0], p1[1], p1[2]);
  pV[3].set(p2[0], p2[1], p2[2]);
  pV[4].set(p3[0], p3[1], p3[2]);

  G4Vector3D v1(pV[2] - pV[1]);
  G4Vector3D v2(pV[3] - pV[1]);
  G4Vector3D v3(pV[4] - pV[1]);

  if (v1.cross(v2).dot(v3) < 0.)
  {
    pV[3].set(p3[0], p3[1], p3[2]);
    pV[4].set(p2[0], p2[1], p2[2]);
  }

  pF[1] = G4Facet(1,2,  3,4,  2,3);
  pF[2] = G4Facet(1,3,  4,4,  3,1);
  pF[3] = G4Facet(1,1,  2,4,  4,2);
  pF[4] = G4Facet(2,1,  3,2,  4,3);
}

HepPolyhedronTet::~HepPolyhedronTet() = default;

HepPolyhedronEllipsoid::HepPolyhedronEllipsoid(G4double ax, G4double by,
                                               G4double cz, G4double zCut1,
                                               G4double zCut2)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronEllipsoid                      Date:    25.02.05 *
 * Author: G.Guerrieri                               Revised:          *
 *         Evgueni Tcherniaev                                 20.01.22 *
 *                                                                     *
 * Function: Constructor of polyhedron for ELLIPSOID                   *
 *                                                                     *
 * Input: ax - semiaxis x                                              *
 *        by - semiaxis y                                              *
 *        cz - semiaxis z                                              *
 *        zCut1 - lower cut plane level (solid lies above this plane)  *
 *        zCut2 - upper cut plane level (solid lies below this plane)  *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  if (zCut1 >= cz || zCut2 <= -cz || zCut1 > zCut2) {
    std::cerr << "HepPolyhedronEllipsoid: wrong zCut1 = " << zCut1
           << " zCut2 = " << zCut2
           << " for given cz = " << cz << std::endl;
    return;
  }
  if (cz <= 0.0) {
    std::cerr << "HepPolyhedronEllipsoid: bad z semi-axis: cz = " << cz
      << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S
  //   generate sphere of radius cz first, then rescale x and y later

  G4double sthe = std::acos(zCut2/cz);
  G4double dthe = std::acos(zCut1/cz) - sthe;
  G4int nds = (GetNumberOfRotationSteps() + 1)/2;
  G4int np1 = G4int(dthe*nds/pi + 0.5) + 1;
  if (np1 <= 1) np1 = 2;
  G4int np2 = 2;

  G4double *zz, *rr;
  zz = new G4double[np1 + np2];
  rr = new G4double[np1 + np2];
  if ((zz == nullptr) || (rr == nullptr))
  {
    G4Exception("HepPolyhedronEllipsoid::HepPolyhedronEllipsoid",
                "greps1002", FatalException, "Out of memory");
  }

  G4double a = dthe/(np1 - 1);
  G4double cosa, sina;
  for (G4int i = 0; i < np1; ++i)
  {
    cosa  = std::cos(sthe + i*a);
    sina  = std::sin(sthe + i*a);
    zz[i] = cz*cosa;
    rr[i] = cz*sina;
  }
  zz[np1 + 0] = zCut2;
  rr[np1 + 0] = 0.;
  zz[np1 + 1] = zCut1;
  rr[np1 + 1] = 0.;

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, 0., twopi, np1, np2, zz, rr, -1, -1);
  SetReferences();

  delete [] zz;
  delete [] rr;

  // rescale x and y vertex coordinates
  G4double kx = ax/cz;
  G4double ky = by/cz;
  G4Point3D* p = pV;
  for (G4int i = 0; i < nvert; ++i, ++p)
  {
    p->setX(p->x()*kx);
    p->setY(p->y()*ky);
  }
}

HepPolyhedronEllipsoid::~HepPolyhedronEllipsoid() = default;

HepPolyhedronEllipticalCone::HepPolyhedronEllipticalCone(G4double ax,
                                                         G4double ay,
                                                         G4double h,
                                                         G4double zTopCut)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronEllipticalCone                 Date:    8.9.2005 *
 * Author: D.Anninos                                 Revised: 9.9.2005 *
 *                                                                     *
 * Function: Constructor for EllipticalCone                            *
 *                                                                     *
 * Input: ax, ay     - X & Y semi axes at z = 0                        *
 *        h          - height of full cone                             *
 *        zTopCut    - Top Cut in Z Axis                               *
 *                                                                     *
 ***********************************************************************/
{
  //   C H E C K   I N P U T   P A R A M E T E R S

  G4int k = 0;
  if ( (ax <= 0.) || (ay <= 0.) || (h <= 0.) || (zTopCut <= 0.) ) { k = 1; }

  if (k != 0) {
    std::cerr << "HepPolyhedronCone: error in input parameters";
    std::cerr << std::endl;
    return;
  }

  //   P R E P A R E   T W O   P O L Y L I N E S

  zTopCut = (h >= zTopCut ? zTopCut : h);

  G4double *zz, *rr;
  zz = new G4double[4];
  rr = new G4double[4];
  zz[0] =   zTopCut;
  zz[1] =  -zTopCut;
  zz[2] =   zTopCut;
  zz[3] =  -zTopCut;
  rr[0] =  (h-zTopCut);
  rr[1] =  (h+zTopCut);
  rr[2] =  0.;
  rr[3] =  0.;

  //   R O T A T E    P O L Y L I N E S

  RotateAroundZ(0, 0., twopi, 2, 2, zz, rr, -1, -1);
  SetReferences();

  delete [] zz;
  delete [] rr;

  // rescale x and y vertex coordinates
 {
   G4Point3D * p= pV;
   for (G4int i=0; i<nvert; i++, p++) {
     p->setX( p->x() * ax );
     p->setY( p->y() * ay );
   }
 }
}

HepPolyhedronEllipticalCone::~HepPolyhedronEllipticalCone() = default;

HepPolyhedronHyperbolicMirror::HepPolyhedronHyperbolicMirror(G4double a,
                                                             G4double h,
                                                             G4double r)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronHyperbolicMirror               Date:  22.02.2020 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Create polyhedron for Hyperbolic mirror                   *
 *                                                                     *
 * Input: a - half-separation                                          *
 *        h - height                                                   *
 *        r - radius                                                   *
 *                                                                     *
 ***********************************************************************/
{
  G4double H = std::abs(h);
  G4double R = std::abs(r);
  G4double A = std::abs(a);
  G4double B = A*R/std::sqrt(2*A*H + H*H);

  //   P R E P A R E   T W O   P O L Y L I N E S

  G4int np1 = (A == 0.) ? 2 : std::max(3, GetNumberOfRotationSteps()/4) + 1;
  G4int np2 = 2;
  G4double maxAng = (A == 0.) ? 0. : std::acosh(1. + H/A);
  G4double delAng = maxAng/(np1 - 1);

  auto zz = new G4double[np1 + np2];
  auto rr = new G4double[np1 + np2];

  // 1st polyline
  zz[0] = H;
  rr[0] = R;
  for (G4int iz = 1; iz < np1 - 1; ++iz)
  {
    G4double ang = maxAng - iz*delAng;
    zz[iz] = A*std::cosh(ang) - A;
    rr[iz] = B*std::sinh(ang);
  }
  zz[np1 - 1] = 0.;
  rr[np1 - 1] = 0.;

  // 2nd polyline
  zz[np1] = H;
  rr[np1] = 0.;
  zz[np1 + 1] = 0.;
  rr[np1 + 1] = 0.;

  //   R O T A T E    P O L Y L I N E S

  G4double phi  = 0.;
  G4double dphi = CLHEP::twopi;
  RotateAroundZ(0, phi, dphi, np1, np2, zz, rr, -1, -1);
  SetReferences();

  delete [] zz;
  delete [] rr;
}

HepPolyhedronHyperbolicMirror::~HepPolyhedronHyperbolicMirror() = default;

HepPolyhedronTetMesh::
HepPolyhedronTetMesh(const std::vector<G4ThreeVector>& tetrahedra)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronTetMesh                        Date:  26.03.2022 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Create polyhedron for tetrahedron mesh                    *
 *                                                                     *
 * Input: tetrahedra - array of tetrahedron vertices, four vertices    *
 *                     per tetrahedron                                 *
 *                                                                     *
 ***********************************************************************/
{
  // Check size of input vector
  G4int nnodes = (G4int)tetrahedra.size();
  if (nnodes == 0)
  {
    std::cerr
      << "HepPolyhedronTetMesh: Empty tetrahedron mesh" << std::endl;
    return;
  }
  G4int ntet = nnodes/4;
  if (nnodes != ntet*4)
  {
    std::cerr << "HepPolyhedronTetMesh: Number of nodes = " << nnodes
              << " in tetrahedron mesh is NOT multiple of 4"
              << std::endl;
    return;
  }

  // Find coincident vertices using hash table techniques.
  // This could be done using std::unordered_map, but the code
  // below runs faster.
  std::vector<G4int> iheads(nnodes, -1);
  std::vector<std::pair<G4int,G4int>> ipairs(nnodes,std::pair(-1,-1));
  for (G4int i = 0; i < nnodes; ++i)
  {
    // Generate hash key
    G4ThreeVector point = tetrahedra[i];
    auto key = std::hash<G4double>()(point.x());
    key ^= std::hash<G4double>()(point.y());
    key ^= std::hash<G4double>()(point.z());
    key %= nnodes;
    // Check head of the list
    if (iheads[key] < 0)
    {
      iheads[key] = i;
      ipairs[i].first = i;
      continue;
    }
    // Loop along the list
    for (G4int icur = iheads[key], iprev = 0;;)
    {
      G4int icheck = ipairs[icur].first;
      if (tetrahedra[icheck] == point)
      {
        ipairs[i].first = icheck; // coincident vertex
        break;
      }
      iprev = icur;
      icur = ipairs[icur].second;
      // Append vertex to the list
      if (icur < 0)
      {
        ipairs[i].first = i;
        ipairs[iprev].second = i;
        break;
      }
    }
  }

  // Create vector of original facets
  struct facet
  {
    G4int i1, i2, i3;
    facet() : i1(0), i2(0), i3(0) {};
    facet(G4int k1, G4int k2, G4int k3) : i1(k1), i2(k2), i3(k3) {};
  };
  G4int nfacets = nnodes;
  std::vector<facet> ifacets(nfacets);
  for (G4int i = 0; i < nfacets; i += 4)
  {
    G4int i0 = ipairs[i + 0].first;
    G4int i1 = ipairs[i + 1].first;
    G4int i2 = ipairs[i + 2].first;
    G4int i3 = ipairs[i + 3].first;
    if (i0 > i1) std::swap(i0, i1);
    if (i0 > i2) std::swap(i0, i2);
    if (i0 > i3) std::swap(i0, i3);
    if (i1 > i2) std::swap(i1, i2);
    if (i1 > i3) std::swap(i1, i3);
    G4ThreeVector e1 = tetrahedra[i1] - tetrahedra[i0];
    G4ThreeVector e2 = tetrahedra[i2] - tetrahedra[i0];
    G4ThreeVector e3 = tetrahedra[i3] - tetrahedra[i0];
    G4double volume = (e1.cross(e2)).dot(e3);
    if (volume > 0.) std::swap(i2, i3);
    ifacets[i + 0] = facet(i0, i1, i2);
    ifacets[i + 1] = facet(i0, i2, i3);
    ifacets[i + 2] = facet(i0, i3, i1);
    ifacets[i + 3] = facet(i1, i3, i2);
  }

  // Find shared facets
  std::fill(iheads.begin(), iheads.end(), -1);
  std::fill(ipairs.begin(), ipairs.end(), std::pair(-1,-1));
  for (G4int i = 0; i < nfacets; ++i)
  {
    // Check head of the list
    G4int key = ifacets[i].i1;
    if (iheads[key] < 0)
    {
      iheads[key] = i;
      ipairs[i].first = i;
      continue;
    }
    // Loop along the list
    G4int i2 = ifacets[i].i2, i3 = ifacets[i].i3;
    for (G4int icur = iheads[key], iprev = -1;;)
    {
      G4int icheck = ipairs[icur].first;
      if (ifacets[icheck].i2 == i3 && ifacets[icheck].i3 == i2)
      {
        if (iprev < 0)
        {
          iheads[key] = ipairs[icur].second;
        }
        else
        {
          ipairs[iprev].second = ipairs[icur].second;
        }
        ipairs[icur].first = -1; // shared facet
        ipairs[icur].second = -1;
        break;
      }
      iprev = icur;
      icur = ipairs[icur].second;
      // Append facet to the list
      if (icur < 0)
      {
        ipairs[i].first = i;
        ipairs[iprev].second = i;
        break;
      }
    }
  }

  // Count vertices and facets skipping shared facets
  std::fill(iheads.begin(), iheads.end(), -1);
  G4int nver = 0, nfac = 0;
  for (G4int i = 0; i < nfacets; ++i)
  {
    if (ipairs[i].first < 0) continue;
    G4int i1 = ifacets[i].i1;
    G4int i2 = ifacets[i].i2;
    G4int i3 = ifacets[i].i3;
    if (iheads[i1] < 0) iheads[i1] = nver++;
    if (iheads[i2] < 0) iheads[i2] = nver++;
    if (iheads[i3] < 0) iheads[i3] = nver++;
    nfac++;
  }

  // Construct polyhedron
  AllocateMemory(nver, nfac);
  for (G4int i = 0; i < nnodes; ++i)
  {
    G4int k = iheads[i];
    if (k >= 0) SetVertex(k + 1, tetrahedra[i]);
  }
  for (G4int i = 0, k = 0; i < nfacets; ++i)
  {
    if (ipairs[i].first < 0) continue;
    G4int i1 = iheads[ifacets[i].i1] + 1;
    G4int i2 = iheads[ifacets[i].i2] + 1;
    G4int i3 = iheads[ifacets[i].i3] + 1;
    SetFacet(++k, i1, i2, i3);
  }
  SetReferences();
}

HepPolyhedronTetMesh::~HepPolyhedronTetMesh() = default;

HepPolyhedronBoxMesh::
HepPolyhedronBoxMesh(G4double sizeX, G4double sizeY, G4double sizeZ,
                     const std::vector<G4ThreeVector>& positions)
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedronBoxMesh                        Date:  07.04.2022 *
 * Author: E.Tcherniaev (E.Chernyaev)                Revised:          *
 *                                                                     *
 * Function: Create polyhedron for box mesh                            *
 *                                                                     *
 * Input: sizeX, sizeY, sizeZ - dimensions of the mesh cell            *
 *        positions - vector of cell centres                           *
 *                                                                     *
 ***********************************************************************/
{
  G4int nbox = (G4int)positions.size();
  if (nbox == 0)
  {
    std::cerr << "HepPolyhedronBoxMesh: Empty box mesh" << std::endl;
    return;
  }
  // compute inverse dimensions
  G4double invx = 1./sizeX, invy = 1./sizeY, invz = 1./sizeZ;
  // find mesh bounding box
  G4ThreeVector pmin = positions[0], pmax = positions[0];
  for (const auto& p: positions)
  {
    if (pmin.x() > p.x()) pmin.setX(p.x());
    if (pmin.y() > p.y()) pmin.setY(p.y());
    if (pmin.z() > p.z()) pmin.setZ(p.z());
    if (pmax.x() < p.x()) pmax.setX(p.x());
    if (pmax.y() < p.y()) pmax.setY(p.y());
    if (pmax.z() < p.z()) pmax.setZ(p.z());
  }
  // find number of voxels
  G4int nx = (pmax.x() - pmin.x())*invx + 1.5;
  G4int ny = (pmax.y() - pmin.y())*invy + 1.5;
  G4int nz = (pmax.z() - pmin.z())*invz + 1.5;
  // create structures for voxels and node indices
  std::vector<char> voxels(nx*ny*nz, 0);
  std::vector<G4int> indices((nx+1)*(ny+1)*(nz+1), 0);
  // mark voxels listed in positions
  G4int kx =  ny*nz, ky = nz;
  for (const auto& p: positions)
  {
    G4int ix = (p.x() - pmin.x())*invx + 0.5;
    G4int iy = (p.y() - pmin.y())*invy + 0.5;
    G4int iz = (p.z() - pmin.z())*invz + 0.5;
    G4int i = ix*kx + iy*ky + iz;
    voxels[i] = 1;
  }
  // count number of vertices and facets
  // set indices
  G4int kvx = (ny + 1)*(nz + 1), kvy = nz + 1;
  G4int nver = 0, nfac = 0;
  for (const auto& p: positions)
  {
    G4int ix = (p.x() - pmin.x())*invx + 0.5;
    G4int iy = (p.y() - pmin.y())*invy + 0.5;
    G4int iz = (p.z() - pmin.z())*invz + 0.5;
    //
    //    011       111
    //      +---–---+
    //      | 001   |   101
    //      |   +---–---+
    //      |   |   |   |
    //      +---|---+   |
    //    010   |   110 |
    //          +-------+
    //        000       100
    //
    G4int vcheck = 0;
    // check (ix - 1) side
    vcheck = (ix == 0) ? 0 : voxels[(ix-1)*kx + iy*ky + iz];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      G4int i3 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      G4int i4 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
    // check (ix + 1) side
    vcheck = (ix == nx - 1) ? 0 : voxels[(ix+1)*kx + iy*ky + iz];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      G4int i2 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i3 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i4 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
    // check (iy - 1) side
    vcheck = (iy == 0) ? 0 : voxels[ix*kx + (iy-1)*ky + iz];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      G4int i3 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i4 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
    // check (iy + 1) side
    vcheck = (iy == ny - 1) ? 0 : voxels[ix*kx + (iy+1)*ky + iz];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      G4int i2 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i4 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
    // check (iz - 1) side
    vcheck = (iz == 0) ? 0 : voxels[ix*kx + iy*ky + (iz-1)];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      G4int i4 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
    // check (iz + 1) side
    vcheck = (iz == nz - 1) ? 0 : voxels[ix*kx + iy*ky + (iz+1)];
    if (vcheck == 0)
    {
      nfac++;
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      G4int i2 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i4 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      if (indices[i1] == 0) indices[i1] = ++nver;
      if (indices[i2] == 0) indices[i2] = ++nver;
      if (indices[i3] == 0) indices[i3] = ++nver;
      if (indices[i4] == 0) indices[i4] = ++nver;
    }
  }
  // Construct polyhedron
  AllocateMemory(nver, nfac);
  G4ThreeVector p0(pmin.x() - 0.5*sizeX, pmin.y() - 0.5*sizeY, pmin.z() - 0.5*sizeZ);
  for (G4int ix = 0; ix <= nx; ++ix)
  {
    for (G4int iy = 0; iy <= ny; ++iy)
    {
      for (G4int iz = 0; iz <= nz; ++iz)
      {
	G4int i = ix*kvx + iy*kvy + iz;
	if (indices[i] == 0) continue;
	SetVertex(indices[i], p0 + G4ThreeVector(ix*sizeX, iy*sizeY, iz*sizeZ));
      }
    }
  }
  nfac = 0;
  for (const auto& p: positions)
  {
    G4int ix = (p.x() - pmin.x())*invx + 0.5;
    G4int iy = (p.y() - pmin.y())*invy + 0.5;
    G4int iz = (p.z() - pmin.z())*invz + 0.5;
    G4int vcheck = 0;
    // check (ix - 1) side
    vcheck = (ix == 0) ? 0 : voxels[(ix-1)*kx + iy*ky + iz];
    if (vcheck == 0)
    {
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      G4int i3 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      G4int i4 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);
    }
    // check (ix + 1) side
    vcheck = (ix == nx - 1) ? 0 : voxels[(ix+1)*kx + iy*ky + iz];
    if (vcheck == 0)
    {
      G4int i1 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      G4int i2 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i3 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i4 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);

    }
    // check (iy - 1) side
    vcheck = (iy == 0) ? 0 : voxels[ix*kx + (iy-1)*ky + iz];
    if (vcheck == 0)
    {
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      G4int i3 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i4 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);
    }
    // check (iy + 1) side
    vcheck = (iy == ny - 1) ? 0 : voxels[ix*kx + (iy+1)*ky + iz];
    if (vcheck == 0)
    {
      G4int i1 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      G4int i2 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i4 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);
    }
    // check (iz - 1) side
    vcheck = (iz == 0) ? 0 : voxels[ix*kx + iy*ky + (iz-1)];
    if (vcheck == 0)
    {
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+0); // 000
      G4int i2 = (ix+0)*kvx + (iy+1)*kvy + (iz+0); // 010
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+0); // 110
      G4int i4 = (ix+1)*kvx + (iy+0)*kvy + (iz+0); // 100
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);
    }
    // check (iz + 1) side
    vcheck = (iz == nz - 1) ? 0 : voxels[ix*kx + iy*ky + (iz+1)];
    if (vcheck == 0)
    {
      G4int i1 = (ix+0)*kvx + (iy+0)*kvy + (iz+1); // 001
      G4int i2 = (ix+1)*kvx + (iy+0)*kvy + (iz+1); // 101
      G4int i3 = (ix+1)*kvx + (iy+1)*kvy + (iz+1); // 111
      G4int i4 = (ix+0)*kvx + (iy+1)*kvy + (iz+1); // 011
      SetFacet(++nfac, indices[i1], indices[i2], indices[i3], indices[i4]);
    }
  }
  SetReferences();
}

HepPolyhedronBoxMesh::~HepPolyhedronBoxMesh() = default;

G4ThreadLocal
G4int HepPolyhedron::fNumberOfRotationSteps = DEFAULT_NUMBER_OF_STEPS;
/***********************************************************************
 *                                                                     *
 * Name: HepPolyhedron::fNumberOfRotationSteps       Date:    24.06.97 *
 * Author: J.Allison (Manchester University)         Revised:          *
 *                                                                     *
 * Function: Number of steps for whole circle                          *
 *                                                                     *
 ***********************************************************************/

#include "BooleanProcessor.src"

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
  G4int ierr;
  BooleanProcessor processor;
  return processor.execute(OP_UNION, *this, p,ierr);
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
  G4int ierr;
  BooleanProcessor processor;
  return processor.execute(OP_INTERSECTION, *this, p,ierr);
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
  G4int ierr;
  BooleanProcessor processor;
  return processor.execute(OP_SUBTRACTION, *this, p,ierr);
}

//NOTE : include the code of HepPolyhedronProcessor here
//       since there is no BooleanProcessor.h

#undef INTERSECTION

#include "HepPolyhedronProcessor.src"
