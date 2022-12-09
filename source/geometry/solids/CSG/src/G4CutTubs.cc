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
// G4CutTubs implementation
//
// 01.06.11 T.Nikitina - Derived from G4Tubs
// 30.10.16 E.Tcherniaev - reimplemented CalculateExtent(),
//                         removed CreateRotatedVetices()
// --------------------------------------------------------------------

#include "G4CutTubs.hh"

#if !defined(G4GEOM_USE_UCTUBS)

#include "G4GeomTools.hh"
#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"
#include "G4GeometryTolerance.hh"
#include "G4BoundingEnvelope.hh"

#include "G4VPVParameterisation.hh"
#include "G4QuickRand.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"

#include "G4AutoLock.hh"

namespace
{
  G4Mutex zminmaxMutex = G4MUTEX_INITIALIZER;
}

using namespace CLHEP;

/////////////////////////////////////////////////////////////////////////
//
// Constructor - check parameters, convert angles so 0<sphi+dpshi<=2_PI
//             - note if pdphi>2PI then reset to 2PI

G4CutTubs::G4CutTubs( const G4String &pName,
                      G4double pRMin, G4double pRMax,
                      G4double pDz,
                      G4double pSPhi, G4double pDPhi,
                      G4ThreeVector pLowNorm,G4ThreeVector pHighNorm )
  : G4CSGSolid(pName), fRMin(pRMin), fRMax(pRMax), fDz(pDz),
    fSPhi(0.), fDPhi(0.), fZMin(0.), fZMax(0.)
{
  kRadTolerance = G4GeometryTolerance::GetInstance()->GetRadialTolerance();
  kAngTolerance = G4GeometryTolerance::GetInstance()->GetAngularTolerance();

  halfCarTolerance = kCarTolerance*0.5;
  halfRadTolerance = kRadTolerance*0.5;
  halfAngTolerance = kAngTolerance*0.5;

  if (pDz<=0) // Check z-len
  {
    std::ostringstream message;
    message << "Negative Z half-length (" << pDz << ") in solid: " << GetName();
    G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002", FatalException, message);
  }
  if ( (pRMin >= pRMax) || (pRMin < 0) ) // Check radii
  {
    std::ostringstream message;
    message << "Invalid values for radii in solid: " << GetName()
            << G4endl
            << "        pRMin = " << pRMin << ", pRMax = " << pRMax;
    G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002", FatalException, message);
  }

  // Check angles
  //
  CheckPhiAngles(pSPhi, pDPhi);

  // Check on Cutted Planes Normals
  // If there is NO CUT, propose to use G4Tubs instead
  //
  if ( ( !pLowNorm.x()) && ( !pLowNorm.y())
    && ( !pHighNorm.x()) && (!pHighNorm.y()) )
  {
    std::ostringstream message;
    message << "Inexisting Low/High Normal to Z plane or Parallel to Z."
            << G4endl
            << "Normals to Z plane are " << pLowNorm << " and "
            << pHighNorm << " in solid: " << GetName() << " \n";
    G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids1001",
                JustWarning, message, "Should use G4Tubs!");
  }

  // If Normal is (0,0,0),means parallel to R, give it value of (0,0,+/-1)
  //
  if (pLowNorm.mag2() == 0.)  { pLowNorm.setZ(-1.); }
  if (pHighNorm.mag2()== 0.)  { pHighNorm.setZ(1.); }

  // Given Normals to Cut Planes have to be an unit vectors.
  // Normalize if it is needed.
  //
  if (pLowNorm.mag2() != 1.)  { pLowNorm  = pLowNorm.unit();  }
  if (pHighNorm.mag2()!= 1.)  { pHighNorm = pHighNorm.unit(); }

  // Normals to cutted planes have to point outside Solid
  //
  if( (pLowNorm.mag2() != 0.) && (pHighNorm.mag2()!= 0. ) )
  {
    if( ( pLowNorm.z()>= 0. ) || ( pHighNorm.z() <= 0.))
    {
      std::ostringstream message;
      message << "Invalid Low or High Normal to Z plane; "
                 "has to point outside Solid." << G4endl
              << "Invalid Norm to Z plane (" << pLowNorm << " or  "
              << pHighNorm << ") in solid: " << GetName();
      G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002",
                  FatalException, message);
    }
  }
  fLowNorm  = pLowNorm;
  fHighNorm = pHighNorm;

  // Check intersection of cut planes, they MUST NOT intersect
  // each other inside the lateral surface
  //
  if(IsCrossingCutPlanes())
  {
    std::ostringstream message;
    message << "Invalid normals to Z plane in solid : " << GetName() << G4endl
            << "Cut planes are crossing inside lateral surface !!!\n"
            << " Solid type: G4CutTubs\n"
            << " Parameters: \n"
            << "    inner radius : " << fRMin/mm << " mm \n"
            << "    outer radius : " << fRMax/mm << " mm \n"
            << "    half length Z: " << fDz/mm << " mm \n"
            << "    starting phi : " << fSPhi/degree << " degrees \n"
            << "    delta phi    : " << fDPhi/degree << " degrees \n"
            << "    low Norm     : " << fLowNorm << "  \n"
            << "    high Norm    : " << fHighNorm;
    G4Exception("G4CutTubs::G4CutTubs()", "GeomSolids0002",
                FatalException, message);
  }
}

///////////////////////////////////////////////////////////////////////
//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//
G4CutTubs::G4CutTubs( __void__& a )
  : G4CSGSolid(a), kRadTolerance(0.), kAngTolerance(0.),
    fRMin(0.), fRMax(0.), fDz(0.), fSPhi(0.), fDPhi(0.), fZMin(0.), fZMax(0.),
    sinCPhi(0.), cosCPhi(0.), cosHDPhi(0.), cosHDPhiOT(0.), cosHDPhiIT(0.),
    sinSPhi(0.), cosSPhi(0.), sinEPhi(0.), cosEPhi(0.),
    halfCarTolerance(0.), halfRadTolerance(0.), halfAngTolerance(0.),
    fLowNorm(G4ThreeVector()), fHighNorm(G4ThreeVector())
{
}

//////////////////////////////////////////////////////////////////////////
//
// Destructor

G4CutTubs::~G4CutTubs()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

G4CutTubs::G4CutTubs(const G4CutTubs& rhs)
  : G4CSGSolid(rhs),
    kRadTolerance(rhs.kRadTolerance), kAngTolerance(rhs.kAngTolerance),
    fRMin(rhs.fRMin), fRMax(rhs.fRMax), fDz(rhs.fDz),
    fSPhi(rhs.fSPhi), fDPhi(rhs.fDPhi),
    fZMin(rhs.fZMin), fZMax(rhs.fZMax),
    sinCPhi(rhs.sinCPhi), cosCPhi(rhs.cosCPhi), cosHDPhi(rhs.cosHDPhi),
    cosHDPhiOT(rhs.cosHDPhiOT), cosHDPhiIT(rhs.cosHDPhiIT),
    sinSPhi(rhs.sinSPhi), cosSPhi(rhs.cosSPhi),
    sinEPhi(rhs.sinEPhi), cosEPhi(rhs.cosEPhi),
    fPhiFullCutTube(rhs.fPhiFullCutTube),
    halfCarTolerance(rhs.halfCarTolerance),
    halfRadTolerance(rhs.halfRadTolerance),
    halfAngTolerance(rhs.halfAngTolerance),
    fLowNorm(rhs.fLowNorm), fHighNorm(rhs.fHighNorm)
{
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

G4CutTubs& G4CutTubs::operator = (const G4CutTubs& rhs)
{
   // Check assignment to self
   //
   if (this == &rhs)  { return *this; }

   // Copy base class data
   //
   G4CSGSolid::operator=(rhs);

   // Copy data
   //
   kRadTolerance = rhs.kRadTolerance; kAngTolerance = rhs.kAngTolerance;
   fRMin = rhs.fRMin; fRMax = rhs.fRMax; fDz = rhs.fDz;
   fSPhi = rhs.fSPhi; fDPhi = rhs.fDPhi;
   fZMin = rhs.fZMin; fZMax = rhs.fZMax;
   sinCPhi = rhs.sinCPhi; cosCPhi = rhs.cosCPhi;
   cosHDPhiOT = rhs.cosHDPhiOT; cosHDPhiIT = rhs.cosHDPhiIT;
   sinSPhi = rhs.sinSPhi; cosSPhi = rhs.cosSPhi;
   sinEPhi = rhs.sinEPhi; cosEPhi = rhs.cosEPhi;
   fPhiFullCutTube = rhs.fPhiFullCutTube;
   halfCarTolerance = rhs.halfCarTolerance;
   halfRadTolerance = rhs.halfRadTolerance;
   halfAngTolerance = rhs.halfAngTolerance;
   fLowNorm = rhs.fLowNorm; fHighNorm = rhs.fHighNorm;

   return *this;
}

//////////////////////////////////////////////////////////////////////////
//
// Get volume

G4double G4CutTubs::GetCubicVolume()
{
  constexpr G4int nphi = 200, nrho = 100;

  if (fCubicVolume == 0.)
  {
    // get parameters
    G4double rmin = GetInnerRadius();
    G4double rmax = GetOuterRadius();
    G4double dz   = GetZHalfLength();
    G4double sphi = GetStartPhiAngle();
    G4double dphi = GetDeltaPhiAngle();

    // calculate volume
    G4double volume = dz*dphi*(rmax*rmax - rmin*rmin);
    if (dphi < twopi) // make recalculation
    {
      // set values for calculation of h - distance between
      // opposite points on bases
      G4ThreeVector nbot = GetLowNorm();
      G4ThreeVector ntop = GetHighNorm();
      G4double nx = nbot.x()/nbot.z() - ntop.x()/ntop.z();
      G4double ny = nbot.y()/nbot.z() - ntop.y()/ntop.z();

      // compute volume by integration
      G4double delrho = (rmax - rmin)/nrho;
      G4double delphi = dphi/nphi;
      volume = 0.;
      for (G4int irho=0; irho<nrho; ++irho)
      {
        G4double r1  = rmin + delrho*irho;
        G4double r2  = rmin + delrho*(irho + 1);
        G4double rho = 0.5*(r1 + r2);
        G4double sector = 0.5*delphi*(r2*r2 - r1*r1);
        for (G4int iphi=0; iphi<nphi; ++iphi)
        {
          G4double phi = sphi + delphi*(iphi + 0.5);
          G4double h = nx*rho*std::cos(phi) + ny*rho*std::sin(phi) + 2.*dz;
          volume += sector*h;
        }
      }
    }
    fCubicVolume = volume;
  }
  return fCubicVolume;
}

//////////////////////////////////////////////////////////////////////////
//
// Get surface area

G4double G4CutTubs::GetSurfaceArea()
{
  constexpr G4int nphi = 400;

  if (fSurfaceArea == 0.)
  {
    // get parameters
    G4double rmin = GetInnerRadius();
    G4double rmax = GetOuterRadius();
    G4double dz   = GetZHalfLength();
    G4double sphi = GetStartPhiAngle();
    G4double dphi = GetDeltaPhiAngle();
    G4ThreeVector nbot = GetLowNorm();
    G4ThreeVector ntop = GetHighNorm();

    // calculate lateral surface area
    G4double sinner = 2.*dz*dphi*rmin;
    G4double souter = 2.*dz*dphi*rmax;
    if (dphi < twopi) // make recalculation
    {
      // set values for calculation of h - distance between
      // opposite points on bases
      G4double nx = nbot.x()/nbot.z() - ntop.x()/ntop.z();
      G4double ny = nbot.y()/nbot.z() - ntop.y()/ntop.z();

      // compute lateral surface area by integration
      G4double delphi = dphi/nphi;
      sinner = 0.;
      souter = 0.;
      for (G4int iphi=0; iphi<nphi; ++iphi)
      {
        G4double phi = sphi + delphi*(iphi + 0.5);
        G4double cosphi = std::cos(phi);
        G4double sinphi = std::sin(phi);
        sinner += rmin*(nx*cosphi + ny*sinphi) + 2.*dz;
        souter += rmax*(nx*cosphi + ny*sinphi) + 2.*dz;
      }
      sinner *= delphi*rmin;
      souter *= delphi*rmax;
    }
    // set surface area
    G4double scut  = (dphi == twopi) ? 0. : 2.*dz*(rmax - rmin);
    G4double szero = 0.5*dphi*(rmax*rmax - rmin*rmin);
    G4double slow  = szero/std::abs(nbot.z());
    G4double shigh = szero/std::abs(ntop.z());
    fSurfaceArea = sinner + souter + 2.*scut + slow + shigh;
  }
  return fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Get bounding box

void G4CutTubs::BoundingLimits(G4ThreeVector& pMin, G4ThreeVector& pMax) const
{
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
  if (dphi < twopi)
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
    G4Exception("G4CutTubs::BoundingLimits()", "GeomMgt0001",
                JustWarning, message);
    DumpInfo();
  }
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate extent under transform and specified limit

G4bool G4CutTubs::CalculateExtent( const EAxis              pAxis,
                                   const G4VoxelLimits&     pVoxelLimit,
                                   const G4AffineTransform& pTransform,
                                         G4double&          pMin,
                                         G4double&          pMax    ) const
{
  G4ThreeVector bmin, bmax;
  G4bool exist;

  // Get bounding box
  BoundingLimits(bmin,bmax);

  // Check bounding box
  G4BoundingEnvelope bbox(bmin,bmax);
#ifdef G4BBOX_EXTENT
  return bbox.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
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
    for (G4int k=0; k<ksteps+2; ++k) { polygons[k] = &pols[k]; }
    G4BoundingEnvelope benv(bmin,bmax,polygons);
    exist = benv.CalculateExtent(pAxis,pVoxelLimit,pTransform,pMin,pMax);
  }
  return exist;
}

//////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface

EInside G4CutTubs::Inside( const G4ThreeVector& p ) const
{
  G4ThreeVector vZ = G4ThreeVector(0,0,fDz);
  EInside in = kInside;

  // Check the lower cut plane
  //
  G4double zinLow =(p+vZ).dot(fLowNorm);
  if (zinLow > halfCarTolerance)  { return kOutside; }

  // Check the higher cut plane
  //
  G4double zinHigh = (p-vZ).dot(fHighNorm);
  if (zinHigh > halfCarTolerance)  { return kOutside; }

  // Check radius
  //
  G4double r2 = p.x()*p.x() + p.y()*p.y() ;

  G4double tolRMin = fRMin - halfRadTolerance;
  G4double tolRMax = fRMax + halfRadTolerance;
  if ( tolRMin < 0 )  { tolRMin = 0; }

  if (r2 < tolRMin*tolRMin || r2 > tolRMax*tolRMax) { return kOutside; }

  // Check Phi cut
  //
  if(!fPhiFullCutTube)
  {
    if ((tolRMin == 0) && (std::fabs(p.x()) <= halfCarTolerance)
                       && (std::fabs(p.y()) <= halfCarTolerance))
    {
      return kSurface;
    }

    G4double phi0 = std::atan2(p.y(),p.x());
    G4double phi1 = phi0 - twopi;
    G4double phi2 = phi0 + twopi;

    in = kOutside;
    G4double sphi = fSPhi - halfAngTolerance;
    G4double ephi = sphi + fDPhi + kAngTolerance;
    if ((phi0  >= sphi && phi0  <= ephi) ||
        (phi1  >= sphi && phi1  <= ephi) ||
        (phi2  >= sphi && phi2  <= ephi)) in = kSurface;
    if (in == kOutside)  { return kOutside; }

    sphi += kAngTolerance;
    ephi -= kAngTolerance;
    if ((phi0  >= sphi && phi0  <= ephi) ||
        (phi1  >= sphi && phi1  <= ephi) ||
        (phi2  >= sphi && phi2  <= ephi)) in = kInside;
    if (in == kSurface)  { return kSurface; }
  }

  // Check on the Surface for Z
  //
  if ((zinLow >= -halfCarTolerance) || (zinHigh >= -halfCarTolerance))
  {
    return kSurface;
  }

  // Check on the Surface for R
  //
  if (fRMin) { tolRMin = fRMin + halfRadTolerance; }
  else       { tolRMin = 0; }
  tolRMax = fRMax - halfRadTolerance;
  if (((r2 <= tolRMin*tolRMin) || (r2 >= tolRMax*tolRMax)) &&
       (r2 >= halfRadTolerance*halfRadTolerance))
  {
    return kSurface;
  }

  return in;
}

///////////////////////////////////////////////////////////////////////////
//
// Return unit normal of surface closest to p
// - note if point on z axis, ignore phi divided sides
// - unsafe if point close to z axis a rmin=0 - no explicit checks

G4ThreeVector G4CutTubs::SurfaceNormal( const G4ThreeVector& p ) const
{
  G4int noSurfaces = 0;
  G4double rho, pPhi;
  G4double distZLow,distZHigh, distRMin, distRMax;
  G4double distSPhi = kInfinity, distEPhi = kInfinity;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  G4ThreeVector norm, sumnorm(0.,0.,0.);
  G4ThreeVector nZ = G4ThreeVector(0, 0, 1.0);
  G4ThreeVector nR, nPs, nPe;

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y());

  distRMin = std::fabs(rho - fRMin);
  distRMax = std::fabs(rho - fRMax);

  // dist to Low Cut
  //
  distZLow =std::fabs((p+vZ).dot(fLowNorm));

  // dist to High Cut
  //
  distZHigh = std::fabs((p-vZ).dot(fHighNorm));

  if (!fPhiFullCutTube)    // Protected against (0,0,z)
  {
    if ( rho > halfCarTolerance )
    {
      pPhi = std::atan2(p.y(),p.x());

      if(pPhi  < fSPhi- halfCarTolerance)           { pPhi += twopi; }
      else if(pPhi > fSPhi+fDPhi+ halfCarTolerance) { pPhi -= twopi; }

      distSPhi = std::fabs(pPhi - fSPhi);
      distEPhi = std::fabs(pPhi - fSPhi - fDPhi);
    }
    else if( !fRMin )
    {
      distSPhi = 0.;
      distEPhi = 0.;
    }
    nPs = G4ThreeVector( sinSPhi, -cosSPhi, 0 );
    nPe = G4ThreeVector( -sinEPhi, cosEPhi, 0 );
  }
  if ( rho > halfCarTolerance ) { nR = G4ThreeVector(p.x()/rho,p.y()/rho,0); }

  if( distRMax <= halfCarTolerance )
  {
    ++noSurfaces;
    sumnorm += nR;
  }
  if( fRMin && (distRMin <= halfCarTolerance) )
  {
    ++noSurfaces;
    sumnorm -= nR;
  }
  if( fDPhi < twopi )
  {
    if (distSPhi <= halfAngTolerance)
    {
      ++noSurfaces;
      sumnorm += nPs;
    }
    if (distEPhi <= halfAngTolerance)
    {
      ++noSurfaces;
      sumnorm += nPe;
    }
  }
  if (distZLow <= halfCarTolerance)
  {
    ++noSurfaces;
    sumnorm += fLowNorm;
  }
  if (distZHigh <= halfCarTolerance)
  {
    ++noSurfaces;
    sumnorm += fHighNorm;
  }
  if ( noSurfaces == 0 )
  {
#ifdef G4CSGDEBUG
    G4Exception("G4CutTubs::SurfaceNormal(p)", "GeomSolids1002",
                JustWarning, "Point p is not on surface !?" );
    G4int oldprc = G4cout.precision(20);
    G4cout<< "G4CutTubs::SN ( "<<p.x()<<", "<<p.y()<<", "<<p.z()<<" ); "
          << G4endl << G4endl;
    G4cout.precision(oldprc) ;
#endif
     norm = ApproxSurfaceNormal(p);
  }
  else if ( noSurfaces == 1 )  { norm = sumnorm; }
  else                         { norm = sumnorm.unit(); }

  return norm;
}

/////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

G4ThreeVector G4CutTubs::ApproxSurfaceNormal( const G4ThreeVector& p ) const
{
  enum ENorm {kNRMin,kNRMax,kNSPhi,kNEPhi,kNZ};

  ENorm side ;
  G4ThreeVector norm ;
  G4double rho, phi ;
  G4double distZLow,distZHigh,distZ;
  G4double distRMin, distRMax, distSPhi, distEPhi, distMin ;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  distRMin = std::fabs(rho - fRMin) ;
  distRMax = std::fabs(rho - fRMax) ;

  //dist to Low Cut
  //
  distZLow =std::fabs((p+vZ).dot(fLowNorm));

  //dist to High Cut
  //
  distZHigh = std::fabs((p-vZ).dot(fHighNorm));
  distZ=std::min(distZLow,distZHigh);

  if (distRMin < distRMax) // First minimum
  {
    if ( distZ < distRMin )
    {
       distMin = distZ ;
       side    = kNZ ;
    }
    else
    {
      distMin = distRMin ;
      side    = kNRMin   ;
    }
  }
  else
  {
    if ( distZ < distRMax )
    {
      distMin = distZ ;
      side    = kNZ   ;
    }
    else
    {
      distMin = distRMax ;
      side    = kNRMax   ;
    }
  }
  if (!fPhiFullCutTube  &&  rho ) // Protected against (0,0,z)
  {
    phi = std::atan2(p.y(),p.x()) ;

    if ( phi < 0 )  { phi += twopi; }

    if ( fSPhi < 0 )
    {
      distSPhi = std::fabs(phi - (fSPhi + twopi))*rho ;
    }
    else
    {
      distSPhi = std::fabs(phi - fSPhi)*rho ;
    }
    distEPhi = std::fabs(phi - fSPhi - fDPhi)*rho ;

    if (distSPhi < distEPhi) // Find new minimum
    {
      if ( distSPhi < distMin )
      {
        side = kNSPhi ;
      }
    }
    else
    {
      if ( distEPhi < distMin )
      {
        side = kNEPhi ;
      }
    }
  }
  switch ( side )
  {
    case kNRMin : // Inner radius
    {
      norm = G4ThreeVector(-p.x()/rho, -p.y()/rho, 0) ;
      break ;
    }
    case kNRMax : // Outer radius
    {
      norm = G4ThreeVector(p.x()/rho, p.y()/rho, 0) ;
      break ;
    }
    case kNZ :    // + or - dz
    {
      if ( distZHigh > distZLow )  { norm = fHighNorm ; }
      else                         { norm = fLowNorm; }
      break ;
    }
    case kNSPhi:
    {
      norm = G4ThreeVector(sinSPhi, -cosSPhi, 0) ;
      break ;
    }
    case kNEPhi:
    {
      norm = G4ThreeVector(-sinEPhi, cosEPhi, 0) ;
      break;
    }
    default:      // Should never reach this case ...
    {
      DumpInfo();
      G4Exception("G4CutTubs::ApproxSurfaceNormal()",
                  "GeomSolids1002", JustWarning,
                  "Undefined side for valid surface normal to solid.");
      break ;
    }
  }
  return norm;
}

////////////////////////////////////////////////////////////////////
//
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - 'if valid' implies tolerant checking of intersection points

G4double G4CutTubs::DistanceToIn( const G4ThreeVector& p,
                                  const G4ThreeVector& v  ) const
{
  G4double snxt = kInfinity ;      // snxt = default return value
  G4double tolORMin2, tolIRMax2 ;  // 'generous' radii squared
  G4double tolORMax2, tolIRMin2;
  const G4double dRmax = 100.*fRMax;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  // Intersection point variables
  //
  G4double Dist, sd=0, xi, yi, zi, rho2, inum, iden, cosPsi, Comp,calf ;
  G4double t1, t2, t3, b, c, d ;     // Quadratic solver variables
  G4double distZLow,distZHigh;
  // Calculate tolerant rmin and rmax

  if (fRMin > kRadTolerance)
  {
    tolORMin2 = (fRMin - halfRadTolerance)*(fRMin - halfRadTolerance) ;
    tolIRMin2 = (fRMin + halfRadTolerance)*(fRMin + halfRadTolerance) ;
  }
  else
  {
    tolORMin2 = 0.0 ;
    tolIRMin2 = 0.0 ;
  }
  tolORMax2 = (fRMax + halfRadTolerance)*(fRMax + halfRadTolerance) ;
  tolIRMax2 = (fRMax - halfRadTolerance)*(fRMax - halfRadTolerance) ;

  // Intersection with ZCut surfaces

  // dist to Low Cut
  //
  distZLow =(p+vZ).dot(fLowNorm);

  // dist to High Cut
  //
  distZHigh = (p-vZ).dot(fHighNorm);

  if ( distZLow >= -halfCarTolerance )
  {
    calf = v.dot(fLowNorm);
    if (calf<0)
    {
      sd = -distZLow/calf;
      if(sd < 0.0)  { sd = 0.0; }

      xi   = p.x() + sd*v.x() ;                // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullCutTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*cosCPhi + yi*sinCPhi ;
          iden   = std::sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)  { return sd ; }
        }
        else
        {
          return sd ;
        }
      }
    }
    else
    {
      if ( sd<halfCarTolerance )
      {
        if(calf>=0) { sd=kInfinity; }
        return sd ;  // On/outside extent, and heading away
      }              // -> cannot intersect
    }
  }

  if(distZHigh >= -halfCarTolerance )
  {
    calf = v.dot(fHighNorm);
    if (calf<0)
    {
      sd = -distZHigh/calf;

      if(sd < 0.0)  { sd = 0.0; }

      xi   = p.x() + sd*v.x() ;                // Intersection coords
      yi   = p.y() + sd*v.y() ;
      rho2 = xi*xi + yi*yi ;

      // Check validity of intersection

      if ((tolIRMin2 <= rho2) && (rho2 <= tolIRMax2))
      {
        if (!fPhiFullCutTube && rho2)
        {
          // Psi = angle made with central (average) phi of shape
          //
          inum   = xi*cosCPhi + yi*sinCPhi ;
          iden   = std::sqrt(rho2) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)  { return sd ; }
        }
        else
        {
          return sd ;
        }
      }
    }
    else
    {
      if ( sd<halfCarTolerance )
      {
        if(calf>=0) { sd=kInfinity; }
        return sd ;  // On/outside extent, and heading away
      }              // -> cannot intersect
    }
  }

  // -> Can not intersect z surfaces
  //
  // Intersection with rmax (possible return) and rmin (must also check phi)
  //
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //            t1                t2                t3

  t1 = 1.0 - v.z()*v.z() ;
  t2 = p.x()*v.x() + p.y()*v.y() ;
  t3 = p.x()*p.x() + p.y()*p.y() ;
  if ( t1 > 0 )        // Check not || to z axis
  {
    b = t2/t1 ;
    c = t3 - fRMax*fRMax ;

    if ((t3 >= tolORMax2) && (t2<0))   // This also handles the tangent case
    {
      // Try outer cylinder intersection, c=(t3-fRMax*fRMax)/t1;

      c /= t1 ;
      d = b*b - c ;

      if (d >= 0)  // If real root
      {
        sd = c/(-b+std::sqrt(d));
        if (sd >= 0)  // If 'forwards'
        {
          if ( sd>dRmax ) // Avoid rounding errors due to precision issues on
          {               // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          }
          // Check z intersection
          //
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              // Z ok. Check phi intersection if reqd
              //
              if (fPhiFullCutTube)
              {
                return sd ;
              }
              else
              {
                xi     = p.x() + sd*v.x() ;
                yi     = p.y() + sd*v.y() ;
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/fRMax ;
                if (cosPsi >= cosHDPhiIT)  { return sd ; }
              }
            }  //  end if std::fabs(zi)
          }
        }    //  end if (sd>=0)
      }      //  end if (d>=0)
    }        //  end if (r>=fRMax)
    else
    {
      // Inside outer radius :
      // check not inside, and heading through tubs (-> 0 to in)
      if ((t3 > tolIRMin2) && (t2 < 0)
       && (std::fabs(p.z()) <= std::fabs(GetCutZ(p))-halfCarTolerance ))
      {
        // Inside both radii, delta r -ve, inside z extent

        if (!fPhiFullCutTube)
        {
          inum   = p.x()*cosCPhi + p.y()*sinCPhi ;
          iden   = std::sqrt(t3) ;
          cosPsi = inum/iden ;
          if (cosPsi >= cosHDPhiIT)
          {
            // In the old version, the small negative tangent for the point
            // on surface was not taken in account, and returning 0.0 ...
            // New version: check the tangent for the point on surface and
            // if no intersection, return kInfinity, if intersection instead
            // return sd.
            //
            c = t3-fRMax*fRMax;
            if ( c<=0.0 )
            {
              return 0.0;
            }
            else
            {
              c = c/t1 ;
              d = b*b-c;
              if ( d>=0.0 )
              {
                snxt = c/(-b+std::sqrt(d)); // using safe solution
                                            // for quadratic equation
                if ( snxt < halfCarTolerance ) { snxt=0; }
                return snxt ;
              }
              else
              {
                return kInfinity;
              }
            }
          }
        }
        else
        {
          // In the old version, the small negative tangent for the point
          // on surface was not taken in account, and returning 0.0 ...
          // New version: check the tangent for the point on surface and
          // if no intersection, return kInfinity, if intersection instead
          // return sd.
          //
          c = t3 - fRMax*fRMax;
          if ( c<=0.0 )
          {
            return 0.0;
          }
          else
          {
            c = c/t1 ;
            d = b*b-c;
            if ( d>=0.0 )
            {
              snxt= c/(-b+std::sqrt(d)); // using safe solution
                                         // for quadratic equation
              if ( snxt < halfCarTolerance ) { snxt=0; }
              return snxt ;
            }
            else
            {
              return kInfinity;
            }
          }
        } // end if   (!fPhiFullCutTube)
      }   // end if   (t3>tolIRMin2)
    }     // end if   (Inside Outer Radius)

    if ( fRMin )    // Try inner cylinder intersection
    {
      c = (t3 - fRMin*fRMin)/t1 ;
      d = b*b - c ;
      if ( d >= 0.0 )  // If real root
      {
        // Always want 2nd root - we are outside and know rmax Hit was bad
        // - If on surface of rmin also need farthest root

        sd =( b > 0. )? c/(-b - std::sqrt(d)) : (-b + std::sqrt(d));
        if (sd >= -10*halfCarTolerance)  // check forwards
        {
          // Check z intersection
          //
          if (sd < 0.0)  { sd = 0.0; }
          if (sd>dRmax) // Avoid rounding errors due to precision issues seen
          {             // 64 bits systems. Split long distances and recompute
            G4double fTerm = sd-std::fmod(sd,dRmax);
            sd = fTerm + DistanceToIn(p+fTerm*v,v);
          }
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              // Z ok. Check phi
              //
              if ( fPhiFullCutTube )
              {
                return sd ;
              }
              else
              {
                cosPsi = (xi*cosCPhi + yi*sinCPhi)/fRMin ;
                if (cosPsi >= cosHDPhiIT)
                {
                  // Good inner radius isect
                  // - but earlier phi isect still possible
                  //
                  snxt = sd ;
                }
              }
            }      //    end if std::fabs(zi)
          }
        }          //    end if (sd>=0)
      }            //    end if (d>=0)
    }              //    end if (fRMin)
  }

  // Phi segment intersection
  //
  // o Tolerant of points inside phi planes by up to kCarTolerance*0.5
  //
  // o NOTE: Large duplication of code between sphi & ephi checks
  //         -> only diffs: sphi -> ephi, Comp -> -Comp and half-plane
  //            intersection check <=0 -> >=0
  //         -> use some form of loop Construct ?
  //
  if ( !fPhiFullCutTube )
  {
    // First phi surface (Starting phi)
    //
    Comp = v.x()*sinSPhi - v.y()*cosSPhi ;

    if ( Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = (p.y()*cosSPhi - p.x()*sinSPhi) ;

      if ( Dist < halfCarTolerance )
      {
        sd = Dist/Comp ;

        if (sd < snxt)
        {
          if ( sd < 0 )  { sd = 0.0; }
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              rho2 = xi*xi + yi*yi ;
              if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                || ( (rho2 >  tolORMin2) && (rho2 <  tolIRMin2)
                  && ( v.y()*cosSPhi - v.x()*sinSPhi >  0 )
                  && ( v.x()*cosSPhi + v.y()*sinSPhi >= 0 )     )
                || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                  && (v.y()*cosSPhi - v.x()*sinSPhi > 0)
                  && (v.x()*cosSPhi + v.y()*sinSPhi < 0) )    )
              {
                // z and r intersections good
                // - check intersecting with correct half-plane
                //
                if ((yi*cosCPhi-xi*sinCPhi) <= halfCarTolerance) { snxt = sd; }
              }
            }   //two Z conditions
          }
        }
      }
    }

    // Second phi surface (Ending phi)
    //
    Comp = -(v.x()*sinEPhi - v.y()*cosEPhi) ;

    if (Comp < 0 )  // Component in outwards normal dirn
    {
      Dist = -(p.y()*cosEPhi - p.x()*sinEPhi) ;

      if ( Dist < halfCarTolerance )
      {
        sd = Dist/Comp ;

        if (sd < snxt)
        {
          if ( sd < 0 )  { sd = 0; }
          zi = p.z() + sd*v.z() ;
          xi = p.x() + sd*v.x() ;
          yi = p.y() + sd*v.y() ;
          if ((-xi*fLowNorm.x()-yi*fLowNorm.y()
               -(zi+fDz)*fLowNorm.z())>-halfCarTolerance)
          {
            if ((-xi*fHighNorm.x()-yi*fHighNorm.y()
                 +(fDz-zi)*fHighNorm.z())>-halfCarTolerance)
            {
              xi   = p.x() + sd*v.x() ;
              yi   = p.y() + sd*v.y() ;
              rho2 = xi*xi + yi*yi ;
              if ( ( (rho2 >= tolIRMin2) && (rho2 <= tolIRMax2) )
                  || ( (rho2 > tolORMin2)  && (rho2 < tolIRMin2)
                    && (v.x()*sinEPhi - v.y()*cosEPhi >  0)
                    && (v.x()*cosEPhi + v.y()*sinEPhi >= 0) )
                  || ( (rho2 > tolIRMax2) && (rho2 < tolORMax2)
                    && (v.x()*sinEPhi - v.y()*cosEPhi > 0)
                    && (v.x()*cosEPhi + v.y()*sinEPhi < 0) ) )
              {
                // z and r intersections good
                // - check intersecting with correct half-plane
                //
                if ( (yi*cosCPhi-xi*sinCPhi) >= -halfCarTolerance )
                {
                  snxt = sd;
                }
              }    //?? >=-halfCarTolerance
            }
          }  // two Z conditions
        }
      }
    }         //  Comp < 0
  }           //  !fPhiFullTube
  if ( snxt<halfCarTolerance )  { snxt=0; }

  return snxt ;
}

//////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance
//
// - Compute the intersection with the z planes
//        - if at valid r, phi, return
//
// -> If point is outer outer radius, compute intersection with rmax
//        - if at valid phi,z return
//
// -> Compute intersection with inner radius, taking largest +ve root
//        - if valid (in z,phi), save intersction
//
//    -> If phi segmented, compute intersections with phi half planes
//        - return smallest of valid phi intersections and
//          inner radius intersection
//
// NOTE:
// - Precalculations for phi trigonometry are Done `just in time'
// - `if valid' implies tolerant checking of intersection points
//   Calculate distance (<= actual) to closest surface of shape from outside
// - Calculate distance to z, radial planes
// - Only to phi planes if outside phi extent
// - Return 0 if point inside

G4double G4CutTubs::DistanceToIn( const G4ThreeVector& p ) const
{
  G4double safRMin,safRMax,safZLow,safZHigh,safePhi,safe,rho,cosPsi;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  // Distance to R
  //
  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;

  safRMin = fRMin- rho ;
  safRMax = rho - fRMax ;

  // Distances to ZCut(Low/High)

  // Dist to Low Cut
  //
  safZLow = (p+vZ).dot(fLowNorm);

  // Dist to High Cut
  //
  safZHigh = (p-vZ).dot(fHighNorm);

  safe = std::max(safZLow,safZHigh);

  if ( safRMin > safe ) { safe = safRMin; }
  if ( safRMax> safe )  { safe = safRMax; }

  // Distance to Phi
  //
  if ( (!fPhiFullCutTube) && (rho) )
   {
     // Psi=angle from central phi to point
     //
     cosPsi = (p.x()*cosCPhi + p.y()*sinCPhi)/rho ;

     if ( cosPsi < cosHDPhi )
     {
       // Point lies outside phi range

       if ( (p.y()*cosCPhi - p.x()*sinCPhi) <= 0 )
       {
         safePhi = std::fabs(p.x()*sinSPhi - p.y()*cosSPhi) ;
       }
       else
       {
         safePhi = std::fabs(p.x()*sinEPhi - p.y()*cosEPhi) ;
       }
       if ( safePhi > safe )  { safe = safePhi; }
     }
   }
   if ( safe < 0 )  { safe = 0; }

   return safe ;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from `inside', allowing for tolerance
// - Only Calc rmax intersection if no valid rmin intersection

G4double G4CutTubs::DistanceToOut( const G4ThreeVector& p,
                                   const G4ThreeVector& v,
                                   const G4bool calcNorm,
                                         G4bool* validNorm,
                                         G4ThreeVector* n ) const
{
  enum ESide {kNull,kRMin,kRMax,kSPhi,kEPhi,kPZ,kMZ};

  ESide side=kNull , sider=kNull, sidephi=kNull ;
  G4double snxt=kInfinity, srd=kInfinity,sz=kInfinity, sphi=kInfinity ;
  G4double deltaR, t1, t2, t3, b, c, d2, roMin2 ;
  G4double distZLow,distZHigh,calfH,calfL;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  // Vars for phi intersection:
  //
  G4double pDistS, compS, pDistE, compE, sphi2, xi, yi, vphi, roi2 ;

  // Z plane intersection
  // Distances to ZCut(Low/High)

  // dist to Low Cut
  //
  distZLow =(p+vZ).dot(fLowNorm);

  // dist to High Cut
  //
  distZHigh = (p-vZ).dot(fHighNorm);

  calfH = v.dot(fHighNorm);
  calfL = v.dot(fLowNorm);

  if (calfH > 0 )
  {
    if ( distZHigh < halfCarTolerance )
    {
      snxt = -distZHigh/calfH ;
      side = kPZ ;
    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector(0,0,1) ;
        *validNorm = true ;
      }
      return snxt = 0 ;
    }
 }
  if ( calfL>0)
  {

    if ( distZLow < halfCarTolerance )
    {
      sz = -distZLow/calfL ;
      if(sz<snxt){
      snxt=sz;
      side = kMZ ;
      }

    }
    else
    {
      if (calcNorm)
      {
        *n         = G4ThreeVector(0,0,-1) ;
        *validNorm = true ;
      }
      return snxt = 0.0 ;
    }
  }
  if((calfH<=0)&&(calfL<=0))
  {
    snxt = kInfinity ;    // Travel perpendicular to z axis
    side = kNull;
  }
  // Radial Intersections
  //
  // Find intersection with cylinders at rmax/rmin
  // Intersection point (xi,yi,zi) on line x=p.x+t*v.x etc.
  //
  // Intersects with x^2+y^2=R^2
  //
  // Hence (v.x^2+v.y^2)t^2+ 2t(p.x*v.x+p.y*v.y)+p.x^2+p.y^2-R^2=0
  //
  //            t1                t2                    t3

  t1   = 1.0 - v.z()*v.z() ;      // since v normalised
  t2   = p.x()*v.x() + p.y()*v.y() ;
  t3   = p.x()*p.x() + p.y()*p.y() ;

  if ( snxt > 10*(fDz+fRMax) )  { roi2 = 2*fRMax*fRMax; }
  else  { roi2 = snxt*snxt*t1 + 2*snxt*t2 + t3; }        // radius^2 on +-fDz

  if ( t1 > 0 ) // Check not parallel
  {
    // Calculate srd, r exit distance

    if ( (t2 >= 0.0) && (roi2 > fRMax*(fRMax + kRadTolerance)) )
    {
      // Delta r not negative => leaving via rmax

      deltaR = t3 - fRMax*fRMax ;

      // NOTE: Should use rho-fRMax<-kRadTolerance*0.5
      // - avoid sqrt for efficiency

      if ( deltaR < -kRadTolerance*fRMax )
      {
        b     = t2/t1 ;
        c     = deltaR/t1 ;
        d2    = b*b-c;
        if( d2 >= 0 ) { srd = c/( -b - std::sqrt(d2)); }
        else          { srd = 0.; }
        sider = kRMax ;
      }
      else
      {
        // On tolerant boundary & heading outwards (or perpendicular to)
        // outer radial surface -> leaving immediately

        if ( calcNorm )
        {
          *n         = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
          *validNorm = true ;
        }
        return snxt = 0 ; // Leaving by rmax immediately
      }
    }
    else if ( t2 < 0. ) // i.e.  t2 < 0; Possible rmin intersection
    {
      roMin2 = t3 - t2*t2/t1 ; // min ro2 of the plane of movement

      if ( fRMin && (roMin2 < fRMin*(fRMin - kRadTolerance)) )
      {
        deltaR = t3 - fRMin*fRMin ;
        b      = t2/t1 ;
        c      = deltaR/t1 ;
        d2     = b*b - c ;

        if ( d2 >= 0 )   // Leaving via rmin
        {
          // NOTE: SHould use rho-rmin>kRadTolerance*0.5
          // - avoid sqrt for efficiency

          if (deltaR > kRadTolerance*fRMin)
          {
            srd = c/(-b+std::sqrt(d2));
            sider = kRMin ;
          }
          else
          {
            if ( calcNorm ) { *validNorm = false; }  // Concave side
            return snxt = 0.0;
          }
        }
        else    // No rmin intersect -> must be rmax intersect
        {
          deltaR = t3 - fRMax*fRMax ;
          c     = deltaR/t1 ;
          d2    = b*b-c;
          if( d2 >=0. )
          {
            srd    = -b + std::sqrt(d2) ;
            sider  = kRMax ;
          }
          else // Case: On the border+t2<kRadTolerance
               //       (v is perpendicular to the surface)
          {
            if (calcNorm)
            {
              *n = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
              *validNorm = true ;
            }
            return snxt = 0.0;
          }
        }
      }
      else if ( roi2 > fRMax*(fRMax + kRadTolerance) )
           // No rmin intersect -> must be rmax intersect
      {
        deltaR = t3 - fRMax*fRMax ;
        b      = t2/t1 ;
        c      = deltaR/t1;
        d2     = b*b-c;
        if( d2 >= 0 )
        {
          srd    = -b + std::sqrt(d2) ;
          sider  = kRMax ;
        }
        else // Case: On the border+t2<kRadTolerance
             //       (v is perpendicular to the surface)
        {
          if (calcNorm)
          {
            *n = G4ThreeVector(p.x()/fRMax,p.y()/fRMax,0) ;
            *validNorm = true ;
          }
          return snxt = 0.0;
        }
      }
    }
    // Phi Intersection

    if ( !fPhiFullCutTube )
    {
      // add angle calculation with correction
      // of the difference in domain of atan2 and Sphi
      //
      vphi = std::atan2(v.y(),v.x()) ;

      if ( vphi < fSPhi - halfAngTolerance  )             { vphi += twopi; }
      else if ( vphi > fSPhi + fDPhi + halfAngTolerance ) { vphi -= twopi; }


      if ( p.x() || p.y() )  // Check if on z axis (rho not needed later)
      {
        // pDist -ve when inside

        pDistS = p.x()*sinSPhi - p.y()*cosSPhi ;
        pDistE = -p.x()*sinEPhi + p.y()*cosEPhi ;

        // Comp -ve when in direction of outwards normal

        compS   = -sinSPhi*v.x() + cosSPhi*v.y() ;
        compE   =  sinEPhi*v.x() - cosEPhi*v.y() ;

        sidephi = kNull;

        if( ( (fDPhi <= pi) && ( (pDistS <= halfCarTolerance)
                              && (pDistE <= halfCarTolerance) ) )
         || ( (fDPhi >  pi) && !((pDistS >  halfCarTolerance)
                              && (pDistE >  halfCarTolerance) ) )  )
        {
          // Inside both phi *full* planes

          if ( compS < 0 )
          {
            sphi = pDistS/compS ;

            if (sphi >= -halfCarTolerance)
            {
              xi = p.x() + sphi*v.x() ;
              yi = p.y() + sphi*v.y() ;

              // Check intersecting with correct half-plane
              // (if not -> no intersect)
              //
              if( (std::fabs(xi)<=kCarTolerance)
               && (std::fabs(yi)<=kCarTolerance) )
              {
                sidephi = kSPhi;
                if (((fSPhi-halfAngTolerance)<=vphi)
                   &&((fSPhi+fDPhi+halfAngTolerance)>=vphi))
                {
                  sphi = kInfinity;
                }
              }
              else if ( yi*cosCPhi-xi*sinCPhi >=0 )
              {
                sphi = kInfinity ;
              }
              else
              {
                sidephi = kSPhi ;
                if ( pDistS > -halfCarTolerance )
                {
                  sphi = 0.0 ; // Leave by sphi immediately
                }
              }
            }
            else
            {
              sphi = kInfinity ;
            }
          }
          else
          {
            sphi = kInfinity ;
          }

          if ( compE < 0 )
          {
            sphi2 = pDistE/compE ;

            // Only check further if < starting phi intersection
            //
            if ( (sphi2 > -halfCarTolerance) && (sphi2 < sphi) )
            {
              xi = p.x() + sphi2*v.x() ;
              yi = p.y() + sphi2*v.y() ;

              if ((std::fabs(xi)<=kCarTolerance)&&(std::fabs(yi)<=kCarTolerance))
              {
                // Leaving via ending phi
                //
                if( !((fSPhi-halfAngTolerance <= vphi)
                     &&(fSPhi+fDPhi+halfAngTolerance >= vphi)) )
                {
                  sidephi = kEPhi ;
                  if ( pDistE <= -halfCarTolerance )  { sphi = sphi2 ; }
                  else                                { sphi = 0.0 ;   }
                }
              }
              else    // Check intersecting with correct half-plane

              if ( (yi*cosCPhi-xi*sinCPhi) >= 0)
              {
                // Leaving via ending phi
                //
                sidephi = kEPhi ;
                if ( pDistE <= -halfCarTolerance ) { sphi = sphi2 ; }
                else                               { sphi = 0.0 ;   }
              }
            }
          }
        }
        else
        {
          sphi = kInfinity ;
        }
      }
      else
      {
        // On z axis + travel not || to z axis -> if phi of vector direction
        // within phi of shape, Step limited by rmax, else Step =0

        if ( (fSPhi - halfAngTolerance <= vphi)
           && (vphi <= fSPhi + fDPhi + halfAngTolerance ) )
        {
          sphi = kInfinity ;
        }
        else
        {
          sidephi = kSPhi ; // arbitrary
          sphi    = 0.0 ;
        }
      }
      if (sphi < snxt)  // Order intersecttions
      {
        snxt = sphi ;
        side = sidephi ;
      }
    }
    if (srd < snxt)  // Order intersections
    {
      snxt = srd ;
      side = sider ;
    }
  }
  if (calcNorm)
  {
    switch(side)
    {
      case kRMax:
        // Note: returned vector not normalised
        // (divide by fRMax for unit vector)
        //
        xi = p.x() + snxt*v.x() ;
        yi = p.y() + snxt*v.y() ;
        *n = G4ThreeVector(xi/fRMax,yi/fRMax,0) ;
        *validNorm = true ;
        break ;

      case kRMin:
        *validNorm = false ;  // Rmin is inconvex
        break ;

      case kSPhi:
        if ( fDPhi <= pi )
        {
          *n         = G4ThreeVector(sinSPhi,-cosSPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kEPhi:
        if (fDPhi <= pi)
        {
          *n = G4ThreeVector(-sinEPhi,cosEPhi,0) ;
          *validNorm = true ;
        }
        else
        {
          *validNorm = false ;
        }
        break ;

      case kPZ:
        *n         = fHighNorm ;
        *validNorm = true ;
        break ;

      case kMZ:
        *n         = fLowNorm ;
        *validNorm = true ;
        break ;

      default:
        G4cout << G4endl ;
        DumpInfo();
        std::ostringstream message;
        G4long oldprc = message.precision(16);
        message << "Undefined side for valid surface normal to solid."
                << G4endl
                << "Position:"  << G4endl << G4endl
                << "p.x() = "   << p.x()/mm << " mm" << G4endl
                << "p.y() = "   << p.y()/mm << " mm" << G4endl
                << "p.z() = "   << p.z()/mm << " mm" << G4endl << G4endl
                << "Direction:" << G4endl << G4endl
                << "v.x() = "   << v.x() << G4endl
                << "v.y() = "   << v.y() << G4endl
                << "v.z() = "   << v.z() << G4endl << G4endl
                << "Proposed distance :" << G4endl << G4endl
                << "snxt = "    << snxt/mm << " mm" << G4endl ;
        message.precision(oldprc) ;
        G4Exception("G4CutTubs::DistanceToOut(p,v,..)", "GeomSolids1002",
                    JustWarning, message);
        break ;
    }
  }
  if ( snxt<halfCarTolerance )  { snxt=0 ; }
  return snxt ;
}

//////////////////////////////////////////////////////////////////////////
//
// Calculate distance (<=actual) to closest surface of shape from inside

G4double G4CutTubs::DistanceToOut( const G4ThreeVector& p ) const
{
  G4double safRMin,safRMax,safZLow,safZHigh,safePhi,safe,rho;
  G4ThreeVector vZ=G4ThreeVector(0,0,fDz);

  rho = std::sqrt(p.x()*p.x() + p.y()*p.y()) ;  // Distance to R

  safRMin =  rho - fRMin ;
  safRMax =  fRMax - rho ;

  // Distances to ZCut(Low/High)

  // Dist to Low Cut
  //
  safZLow = std::fabs((p+vZ).dot(fLowNorm));

  // Dist to High Cut
  //
  safZHigh = std::fabs((p-vZ).dot(fHighNorm));
  safe = std::min(safZLow,safZHigh);

  if ( safRMin < safe ) { safe = safRMin; }
  if ( safRMax< safe )  { safe = safRMax; }

  // Check if phi divided, Calc distances closest phi plane
  //
  if ( !fPhiFullCutTube )
  {
    if ( p.y()*cosCPhi-p.x()*sinCPhi <= 0 )
    {
      safePhi = -(p.x()*sinSPhi - p.y()*cosSPhi) ;
    }
    else
    {
      safePhi = (p.x()*sinEPhi - p.y()*cosEPhi) ;
    }
    if (safePhi < safe)  { safe = safePhi ; }
  }
  if ( safe < 0 )  { safe = 0; }

  return safe ;
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

G4GeometryType G4CutTubs::GetEntityType() const
{
  return G4String("G4CutTubs");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
G4VSolid* G4CutTubs::Clone() const
{
  return new G4CutTubs(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& G4CutTubs::StreamInfo( std::ostream& os ) const
{
  G4long oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: G4CutTubs\n"
     << " Parameters: \n"
     << "    inner radius : " << fRMin/mm << " mm \n"
     << "    outer radius : " << fRMax/mm << " mm \n"
     << "    half length Z: " << fDz/mm << " mm \n"
     << "    starting phi : " << fSPhi/degree << " degrees \n"
     << "    delta phi    : " << fDPhi/degree << " degrees \n"
     << "    low Norm     : " << fLowNorm     << "  \n"
     << "    high Norm    : "  <<fHighNorm    << "  \n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnSurface

G4ThreeVector G4CutTubs::GetPointOnSurface() const
{
  // Set min and max z
  if (fZMin == 0. && fZMax == 0.)
  {
    G4AutoLock l(&zminmaxMutex);
    G4ThreeVector bmin, bmax;
    BoundingLimits(bmin,bmax);
    fZMin = bmin.z();
    fZMax = bmax.z();
    l.unlock();
  }

  // Set parameters
  G4double hmax = fZMax - fZMin;
  G4double sphi = fSPhi;
  G4double dphi = fDPhi;
  G4double rmin = fRMin;
  G4double rmax = fRMax;
  G4double rrmax = rmax*rmax;
  G4double rrmin = rmin*rmin;

  G4ThreeVector nbot = GetLowNorm();
  G4ThreeVector ntop = GetHighNorm();

  // Set array of surface areas
  G4double sbase = 0.5*dphi*(rrmax - rrmin);
  G4double sbot = sbase/std::abs(nbot.z());
  G4double stop = sbase/std::abs(ntop.z());
  G4double scut = (dphi == twopi) ? 0. : hmax*(rmax - rmin);
  G4double ssurf[6] = { scut, scut, sbot, stop, dphi*rmax*hmax, dphi*rmin*hmax };
  ssurf[1] += ssurf[0];
  ssurf[2] += ssurf[1];
  ssurf[3] += ssurf[2];
  ssurf[4] += ssurf[3];
  ssurf[5] += ssurf[4];

  constexpr G4int ntry = 100000;
  for (G4int i=0; i<ntry; ++i)
  {
    // Select surface
    G4double select = ssurf[5]*G4QuickRand();
    G4int k = 5;
    k -= (select <= ssurf[4]);
    k -= (select <= ssurf[3]);
    k -= (select <= ssurf[2]);
    k -= (select <= ssurf[1]);
    k -= (select <= ssurf[0]);

    // Generate point on selected surface (rejection sampling)
    G4ThreeVector p(0,0,0);
    switch(k)
    {
      case 0: // cut at start phi
      {
        G4double r = rmin + (rmax - rmin)*G4QuickRand();
        p.set(r*cosSPhi, r*sinSPhi, fZMin + hmax*G4QuickRand());
        break;
      }
      case 1: // cut at end phi
      {
        G4double r = rmin + (rmax - rmin)*G4QuickRand();
        p.set(r*cosEPhi, r*sinEPhi, fZMin + hmax*G4QuickRand());
        break;
      }
      case 2: // base at low z
      {
        G4double r = std::sqrt(rrmin + (rrmax - rrmin)*G4QuickRand());
        G4double phi = sphi + dphi*G4QuickRand();
        G4double x = r*std::cos(phi);
        G4double y = r*std::sin(phi);
        G4double z = -fDz - (x*nbot.x() + y*nbot.y())/nbot.z();
        return G4ThreeVector(x, y, z);
      }
      case 3: // base at high z
      {
        G4double r = std::sqrt(rrmin + (rrmax - rrmin)*G4QuickRand());
        G4double phi = sphi + dphi*G4QuickRand();
        G4double x = r*std::cos(phi);
        G4double y = r*std::sin(phi);
        G4double z = fDz - (x*ntop.x() + y*ntop.y())/ntop.z();
        return G4ThreeVector(x, y, z);
      }
      case 4: // external lateral surface
      {
        G4double phi = sphi + dphi*G4QuickRand();
        G4double z = fZMin + hmax*G4QuickRand();
        G4double x = rmax*std::cos(phi);
        G4double y = rmax*std::sin(phi);
        p.set(x, y, z);
        break;
      }
      case 5: // internal lateral surface
      {
        G4double phi = sphi + dphi*G4QuickRand();
        G4double z = fZMin + hmax*G4QuickRand();
        G4double x = rmin*std::cos(phi);
        G4double y = rmin*std::sin(phi);
        p.set(x, y, z);
        break;
      }
    }
    if ((ntop.dot(p) - fDz*ntop.z()) > 0.) continue;
    if ((nbot.dot(p) + fDz*nbot.z()) > 0.) continue;
    return p;
  }
  // Just in case, if all attempts to generate a point have failed
  // Normally should never happen
  G4double x = rmax*std::cos(sphi + 0.5*dphi);
  G4double y = rmax*std::sin(sphi + 0.5*dphi);
  G4double z = fDz - (x*ntop.x() + y*ntop.y())/ntop.z();
  return G4ThreeVector(x, y, z);
}

///////////////////////////////////////////////////////////////////////////
//
// Methods for visualisation

void G4CutTubs::DescribeYourselfTo ( G4VGraphicsScene& scene ) const
{
  scene.AddSolid (*this) ;
}

G4Polyhedron* G4CutTubs::CreatePolyhedron () const
{
  typedef G4double G4double3[3];
  typedef G4int G4int4[4];

  G4Polyhedron *ph  = new G4Polyhedron;
  G4Polyhedron *ph1 = new G4PolyhedronTubs (fRMin, fRMax, fDz, fSPhi, fDPhi);
  G4int nn=ph1->GetNoVertices();
  G4int nf=ph1->GetNoFacets();
  G4double3* xyz = new G4double3[nn];  // number of nodes
  G4int4*  faces = new G4int4[nf] ;    // number of faces

  for(G4int i=0; i<nn; ++i)
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
  for(G4int i=0; i<nf ; ++i)
  {
    ph1->GetFacet(i+1,n,iNodes,iEdge);
    for(G4int k=0; k<n; ++k)
    {
      faces[i][k]=iNodes[k];
    }
    for(G4int k=n; k<4; ++k)
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

// Auxilary Methods for Solid

//////////////////////////////////////////////////////////////////////////
//
// Check set of points on the outer lateral surface and return true
// if the cut planes are crossing inside the surface
//

G4bool G4CutTubs::IsCrossingCutPlanes() const
{
  constexpr G4int npoints = 30;

  // set values for calculation of h - distance between
  // opposite points on bases
  G4ThreeVector nbot = GetLowNorm();
  G4ThreeVector ntop = GetHighNorm();
  if (std::abs(nbot.z()) < kCarTolerance) return true;
  if (std::abs(ntop.z()) < kCarTolerance) return true;
  G4double nx = nbot.x()/nbot.z() - ntop.x()/ntop.z();
  G4double ny = nbot.y()/nbot.z() - ntop.y()/ntop.z();

  // check points
  G4double cosphi = GetCosStartPhi();
  G4double sinphi = GetSinStartPhi();
  G4double delphi = GetDeltaPhiAngle()/npoints;
  G4double cosdel = std::cos(delphi);
  G4double sindel = std::sin(delphi);
  G4double hzero = 2.*GetZHalfLength()/GetOuterRadius();
  for (G4int i=0; i<npoints+1; ++i)
  {
    G4double h = nx*cosphi + ny*sinphi + hzero;
    if (h < 0.) return true;
    G4double sintmp = sinphi;
    sinphi = sintmp*cosdel + cosphi*sindel;
    cosphi = cosphi*cosdel - sintmp*sindel;
  }
  return false;
}

///////////////////////////////////////////////////////////////////////////
//
// Return real Z coordinate of point on Cutted +/- fDZ plane

G4double G4CutTubs::GetCutZ(const G4ThreeVector& p) const
{
  G4double newz = p.z();  // p.z() should be either +fDz or -fDz
  if (p.z()<0)
  {
    if(fLowNorm.z()!=0.)
    {
       newz = -fDz-(p.x()*fLowNorm.x()+p.y()*fLowNorm.y())/fLowNorm.z();
    }
  }
  else
  {
    if(fHighNorm.z()!=0.)
    {
       newz = fDz-(p.x()*fHighNorm.x()+p.y()*fHighNorm.y())/fHighNorm.z();
    }
  }
  return newz;
}
#endif
