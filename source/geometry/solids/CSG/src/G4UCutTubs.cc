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
// Implementation for G4UCutTubs wrapper class
// --------------------------------------------------------------------

#include "G4CutTubs.hh"
#include "G4UCutTubs.hh"

#if ( defined(G4GEOM_USE_USOLIDS) || defined(G4GEOM_USE_PARTIAL_USOLIDS) )

#include "G4GeomTools.hh"
#include "G4AffineTransform.hh"
#include "G4VPVParameterisation.hh"
#include "G4BoundingEnvelope.hh"

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4UCutTubs::G4UCutTubs( const G4String& pName,
                              G4double pRMin, G4double pRMax,
                              G4double pDz,
                              G4double pSPhi, G4double pDPhi,
                              G4ThreeVector pLowNorm,
                              G4ThreeVector pHighNorm )
  : Base_t(pName, pRMin, pRMax, pDz, pSPhi, pDPhi,
           pLowNorm.x(), pLowNorm.y(), pLowNorm.z(),
           pHighNorm.x(), pHighNorm.y(), pHighNorm.z())
{
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4UCutTubs::G4UCutTubs( __void__& a )
  : Base_t(a)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4UCutTubs::~G4UCutTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4UCutTubs::G4UCutTubs(const G4UCutTubs& rhs)
  : Base_t(rhs)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4UCutTubs& G4UCutTubs::operator = (const G4UCutTubs& rhs) 
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   Base_t::operator=(rhs);

   return *this;
}

/////////////////////////////////////////////////////////////////////////
//
// Accessors and modifiers

G4double G4UCutTubs::GetInnerRadius() const
{
  return rmin();
}
G4double G4UCutTubs::GetOuterRadius() const
{
  return rmax();
}
G4double G4UCutTubs::GetZHalfLength() const
{
  return z();
}
G4double G4UCutTubs::GetStartPhiAngle() const
{
  return sphi();
}
G4double G4UCutTubs::GetDeltaPhiAngle() const
{
  return dphi();
}
G4double G4UCutTubs::GetSinStartPhi() const
{
  return std::sin(GetStartPhiAngle());
}
G4double G4UCutTubs::GetCosStartPhi() const
{
  return std::cos(GetStartPhiAngle());
}
G4double G4UCutTubs::GetSinEndPhi() const
{
  return std::sin(GetStartPhiAngle()+GetDeltaPhiAngle());
}
G4double G4UCutTubs::GetCosEndPhi() const
{
  return std::cos(GetStartPhiAngle()+GetDeltaPhiAngle());
}
G4ThreeVector G4UCutTubs::GetLowNorm  () const
{
  U3Vector lc = BottomNormal();
  return G4ThreeVector(lc.x(), lc.y(), lc.z());
}
G4ThreeVector G4UCutTubs::GetHighNorm () const
{
  U3Vector hc = TopNormal();
  return G4ThreeVector(hc.x(), hc.y(), hc.z());
} 

void G4UCutTubs::SetInnerRadius(G4double newRMin)
{
  SetRMin(newRMin);
  fRebuildPolyhedron = true;
}
void G4UCutTubs::SetOuterRadius(G4double newRMax)
{
  SetRMax(newRMax);
  fRebuildPolyhedron = true;
}
void G4UCutTubs::SetZHalfLength(G4double newDz)
{
  SetDz(newDz);
  fRebuildPolyhedron = true;
}
void G4UCutTubs::SetStartPhiAngle(G4double newSPhi, G4bool)
{
  SetSPhi(newSPhi);
  fRebuildPolyhedron = true;
}
void G4UCutTubs::SetDeltaPhiAngle(G4double newDPhi)
{
  SetDPhi(newDPhi);
  fRebuildPolyhedron = true;
}

/////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object

G4VSolid* G4UCutTubs::Clone() const
{
  return new G4UCutTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4UCutTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
  static G4bool checkBBox = true;

  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();
  G4double dz   = GetZHalfLength();
  G4double dphi = GetDeltaPhiAngle();

  G4double sinSphi = GetSinStartPhi(); 
  G4double cosSphi = GetCosStartPhi(); 
  G4double sinEphi = GetSinEndPhi(); 
  G4double cosEphi = GetCosEndPhi(); 

  G4ThreeVector norm;
  G4double mag, topx, topy, dists, diste; 
  G4bool iftop;

  // Find Zmin
  //
  G4double zmin;
  norm = GetLowNorm();
  mag  = std::sqrt(norm.x()*norm.x() + norm.y()*norm.y());
  topx = (mag == 0) ? 0 : -rmax*norm.x()/mag; 
  topy = (mag == 0) ? 0 : -rmax*norm.y()/mag; 
  dists =  sinSphi*topx - cosSphi*topy;
  diste = -sinEphi*topx + cosEphi*topy;
  if (dphi > pi)
  {
    iftop = true;
    if (dists > 0 && diste > 0)iftop = false;
  }
  else
  {
    iftop = false;
    if (dists <= 0 && diste <= 0) iftop = true;
  }
  if (iftop)
  {
    zmin = -(norm.x()*topx + norm.y()*topy)/norm.z() - dz;
  }
  else
  {
    G4double z1 = -rmin*(norm.x()*cosSphi + norm.y()*sinSphi)/norm.z() - dz;  
    G4double z2 = -rmin*(norm.x()*cosEphi + norm.y()*sinEphi)/norm.z() - dz;  
    G4double z3 = -rmax*(norm.x()*cosSphi + norm.y()*sinSphi)/norm.z() - dz;  
    G4double z4 = -rmax*(norm.x()*cosEphi + norm.y()*sinEphi)/norm.z() - dz;  
    zmin = std::min(std::min(std::min(z1,z2),z3),z4);
  }

  // Find Zmax
  //
  G4double zmax;
  norm = GetHighNorm();
  mag  = std::sqrt(norm.x()*norm.x() + norm.y()*norm.y());
  topx = (mag == 0) ? 0 : -rmax*norm.x()/mag; 
  topy = (mag == 0) ? 0 : -rmax*norm.y()/mag; 
  dists =  sinSphi*topx - cosSphi*topy;
  diste = -sinEphi*topx + cosEphi*topy;
  if (dphi > pi)
  {
    iftop = true;
    if (dists > 0 && diste > 0) iftop = false;
  }
  else
  {
    iftop = false;
    if (dists <= 0 && diste <= 0) iftop = true;
  }
  if (iftop)
  {
    zmax = -(norm.x()*topx + norm.y()*topy)/norm.z() + dz;
  }
  else
  {
    G4double z1 = -rmin*(norm.x()*cosSphi + norm.y()*sinSphi)/norm.z() + dz;  
    G4double z2 = -rmin*(norm.x()*cosEphi + norm.y()*sinEphi)/norm.z() + dz;  
    G4double z3 = -rmax*(norm.x()*cosSphi + norm.y()*sinSphi)/norm.z() + dz;  
    G4double z4 = -rmax*(norm.x()*cosEphi + norm.y()*sinEphi)/norm.z() + dz;  
    zmax = std::max(std::max(std::max(z1,z2),z3),z4);
  }

  // Find bounding box
  //
  if (GetDeltaPhiAngle() < twopi)
  {
    G4TwoVector vmin,vmax;
    G4GeomTools::DiskExtent(rmin,rmax,
                            GetSinStartPhi(),GetCosStartPhi(),
                            GetSinEndPhi(),GetCosEndPhi(),
                            vmin,vmax);
    pMin.set(vmin.x(),vmin.y(), zmin);
    pMax.set(vmax.x(),vmax.y(), zmax);
  }
  else
  {
    pMin.set(-rmax,-rmax, zmin);
    pMax.set( rmax, rmax, zmax);
  }

  // Check correctness of the bounding box
  //
  if (pMin.x() >= pMax.x() || pMin.y() >= pMax.y() || pMin.z() >= pMax.z())
  {
    std::ostringstream message;
    message << "Bad bounding box (min >= max) for solid: "
            << GetName() << " !"
            << "\npMin = " << pMin
            << "\npMax = " << pMax;
    G4Exception("G4CUutTubs::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    StreamInfo(G4cout);
  }

  // Check consistency of bounding boxes
  //
  if (checkBBox)
  {
    U3Vector vmin, vmax;
    Extent(vmin,vmax);
    if (std::abs(pMin.x()-vmin.x()) > kCarTolerance ||
        std::abs(pMin.y()-vmin.y()) > kCarTolerance ||
        std::abs(pMin.z()-vmin.z()) > kCarTolerance ||
        std::abs(pMax.x()-vmax.x()) > kCarTolerance ||
        std::abs(pMax.y()-vmax.y()) > kCarTolerance ||
        std::abs(pMax.z()-vmax.z()) > kCarTolerance)
    {
      std::ostringstream message;
      message << "Inconsistency in bounding boxes for solid: "
              << GetName() << " !"
              << "\nBBox min: wrapper = " << pMin << " solid = " << vmin
              << "\nBBox max: wrapper = " << pMax << " solid = " << vmax;
      G4Exception("G4UCutTubs::BoundingLimits()", "GeomMgt0001",
                  JustWarning, message);
      checkBBox = false;
    }
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool
G4UCutTubs::CalculateExtent(const EAxis pAxis,
                            const G4VoxelLimits& pVoxelLimit,
                            const G4AffineTransform& pTransform,
                                  G4double& pMin, G4double& pMax) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Check bounding box
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  if (true) return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
#endif
  if (bbox.BoundingBoxVsVoxelLimits(pAxis,pVoxelLimit,pTransform,pMin,pMax))
  {
    return exist = (pMin < pMax) ? true : false;
  }

  // Get parameters of the solid
  G4double rmin = GetInnerRadius();
  G4double rmax = GetOuterRadius();
  G4double dphi = GetDeltaPhiAngle();
  G4double zmin = bmin.z();
  G4double zmax = bmax.z();

  // Find bounding envelope and calculate extent
  //
  const G4int NSTEPS = 24;            // number of steps for whole circle
  G4double astep  = twopi/NSTEPS;     // max angle for one step
  G4int    ksteps = (dphi <= astep) ? 1 : (G4int)((dphi-deg)/astep) + 1;
  G4double ang    = dphi/ksteps;

  G4double sinHalf = std::sin(0.5*ang);
  G4double cosHalf = std::cos(0.5*ang);
  G4double sinStep = 2.*sinHalf*cosHalf;
  G4double cosStep = 1. - 2.*sinHalf*sinHalf;
  G4double rext    = rmax/cosHalf;

  // bounding envelope for full cylinder consists of two polygons,
  // in other cases it is a sequence of quadrilaterals
  if (rmin == 0 && dphi == twopi)
  {
    G4double sinCur = sinHalf;
    G4double cosCur = cosHalf;

    G4ThreeVectorList baseA(NSTEPS),baseB(NSTEPS);
    for (G4int k=0; k<NSTEPS; ++k)
    {
      baseA[k].set(rext*cosCur,rext*sinCur,zmin);
      baseB[k].set(rext*cosCur,rext*sinCur,zmax);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    std::vector<const G4ThreeVectorList *> polygons(2);
    polygons[0] = &baseA;
    polygons[1] = &baseB;
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }
  else
  {
    G4double sinStart = GetSinStartPhi();
    G4double cosStart = GetCosStartPhi();
    G4double sinEnd   = GetSinEndPhi();
    G4double cosEnd   = GetCosEndPhi();
    G4double sinCur   = sinStart*cosHalf + cosStart*sinHalf;
    G4double cosCur   = cosStart*cosHalf - sinStart*sinHalf;

    // set quadrilaterals
    G4ThreeVectorList pols[NSTEPS+2];
    for (G4int k=0; k<ksteps+2; ++k) pols[k].resize(4);
    pols[0][0].set(rmin*cosStart,rmin*sinStart,zmax);
    pols[0][1].set(rmin*cosStart,rmin*sinStart,zmin);
    pols[0][2].set(rmax*cosStart,rmax*sinStart,zmin);
    pols[0][3].set(rmax*cosStart,rmax*sinStart,zmax);
    for (G4int k=1; k<ksteps+1; ++k)
    {
      pols[k][0].set(rmin*cosCur,rmin*sinCur,zmax);
      pols[k][1].set(rmin*cosCur,rmin*sinCur,zmin);
      pols[k][2].set(rext*cosCur,rext*sinCur,zmin);
      pols[k][3].set(rext*cosCur,rext*sinCur,zmax);

      G4double sinTmp = sinCur;
      sinCur = sinCur*cosStep + cosCur*sinStep;
      cosCur = cosCur*cosStep - sinTmp*sinStep;
    }
    pols[ksteps+1][0].set(rmin*cosEnd,rmin*sinEnd,zmax);
    pols[ksteps+1][1].set(rmin*cosEnd,rmin*sinEnd,zmin);
    pols[ksteps+1][2].set(rmax*cosEnd,rmax*sinEnd,zmin);
    pols[ksteps+1][3].set(rmax*cosEnd,rmax*sinEnd,zmax);

    // set envelope and calculate extent
    std::vector<const G4ThreeVectorList *> polygons;
    polygons.resize(ksteps+2);
    for (G4int k=0; k<ksteps+2; ++k) polygons[k] = &pols[k];
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }
  return exist;
}

///////////////////////////////////////////////////////////////////////////
//
// Return real Z coordinate of point on Cutted +/- fDZ plane

G4double G4UCutTubs::GetCutZ(const G4ThreeVector& p) const
{
  G4double newz = p.z();  // p.z() should be either +fDz or -fDz
  G4ThreeVector fLowNorm = GetLowNorm();
  G4ThreeVector fHighNorm = GetHighNorm();
  
  if (p.z()<0)
  {
    if(fLowNorm.z()!=0.)
    {
       newz = -GetZHalfLength()
            - (p.x()*fLowNorm.x()+p.y()*fLowNorm.y())/fLowNorm.z();
    }
  }
  else
  {
    if(fHighNorm.z()!=0.)
    {
       newz = GetZHalfLength()
            - (p.x()*fHighNorm.x()+p.y()*fHighNorm.y())/fHighNorm.z();
    }
  }
  return newz;
}

//////////////////////////////////////////////////////////////////////////
//
// Create polyhedron for visualization
//
G4Polyhedron* G4UCutTubs::CreatePolyhedron() const
{
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];

  G4Polyhedron *ph  = new G4Polyhedron;
  G4Polyhedron *ph1 = new G4PolyhedronTubs(GetInnerRadius(),
                                           GetOuterRadius(),
                                           GetZHalfLength(),
                                           GetStartPhiAngle(),
                                           GetDeltaPhiAngle());
  G4int nn=ph1->GetNoVertices();
  G4int nf=ph1->GetNoFacets();
  G4double3* xyz = new G4double3[nn];  // number of nodes 
  G4int4*  faces = new G4int4[nf] ;    // number of faces
  G4double fDz = GetZHalfLength();

  for(G4int i=0;i<nn;++i)
  {
    xyz[i][0]=ph1->GetVertex(i+1).x();
    xyz[i][1]=ph1->GetVertex(i+1).y();
    G4double tmpZ=ph1->GetVertex(i+1).z();
    if(tmpZ>=fDz-kCarTolerance)
    {
      xyz[i][2]=GetCutZ(G4ThreeVector(xyz[i][0],xyz[i][1],fDz));
    }
    else if(tmpZ<=-fDz+kCarTolerance)
    {
      xyz[i][2]=GetCutZ(G4ThreeVector(xyz[i][0],xyz[i][1],-fDz));
    }
    else
    {
      xyz[i][2]=tmpZ;
    }
  }
  G4int iNodes[4];
  G4int *iEdge=0;
  G4int n;
  for(G4int i=0;i<nf;++i)
  {
    ph1->GetFacet(i+1,n,iNodes,iEdge);
    for(G4int k=0;k<n;++k)
    {
      faces[i][k]=iNodes[k];
    }
    for(G4int k=n;k<4;++k)
    {
      faces[i][k]=0;
    }
  }
  ph->createPolyhedron(nn,nf,xyz,faces);

  delete [] xyz;
  delete [] faces;
  delete ph1;

  return ph;
}

#endif  // G4GEOM_USE_USOLIDS
