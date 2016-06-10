//
// ********************************************************************
// * This Software is part of the AIDA Unified Solids Library package *
// * See: https://aidasoft.web.cern.ch/USolids                        *
// ********************************************************************
//
// $Id:$
//
// --------------------------------------------------------------------
//
// UTrap
//
// 12.02.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UTrap.hh"

using namespace std;

////////////////////////////////////////////////////////////////////////
//
// Accuracy of coplanarity

const double kCoplanar_Tolerance = 1E-4 ;

//////////////////////////////////////////////////////////////////////////
//
// Private enum: Not for external use

enum Eside {kUndef, ks0, ks1, ks2, ks3, kPZ, kMZ};

//////////////////////////////////////////////////////////////////////////
//
// Constructor - check and Set half-widths as well as angles:
// final check of coplanarity

UTrap::UTrap(const std::string& pName,
             double pDz,
             double pTheta, double pPhi,
             double pDy1, double pDx1, double pDx2,
             double pAlp1,
             double pDy2, double pDx3, double pDx4,
             double pAlp2)
  : VUSolid(pName)
{
  if (pDz <= 0 || pDy1 <= 0 || pDx1 <= 0 ||
      pDx2 <= 0 || pDy2 <= 0 || pDx3 <= 0 || pDx4 <= 0)
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName() << std::endl
            << "				X - "
            << pDx1 << ", " << pDx2 << ", " << pDx3 << ", " << pDx4 << std::endl
            << "					Y - " << pDy1 << ", " << pDy2 << std::endl
            << "					Z - " << pDz;
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  fDz = pDz;
  fTthetaCphi = std::tan(pTheta) * std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta) * std::sin(pPhi);

  fDy1 = pDy1;
  fDx1 = pDx1;
  fDx2 = pDx2;
  fTalpha1 = std::tan(pAlp1);

  fDy2 = pDy2;
  fDx3 = pDx3;
  fDx4 = pDx4;
  fTalpha2 = std::tan(pAlp2);

  MakePlanes();
  fCubicVolume=0.;
  fSurfaceArea=0.;
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor - Design of trapezoid based on 8 UVector3 parameters,
// which are its vertices. Checking of planarity with preparation of
// fPlanes[] and than calculation of other members

UTrap::UTrap(const std::string& pName,
             const UVector3 pt[8])
  : VUSolid(pName)
{
  SetPlanes(pt);
  fCubicVolume=0.;
  fSurfaceArea=0.;
}

//////////////////////////////////////////////////////////////////////////////
//
// Constructor for Right Angular Wedge from STEP

UTrap::UTrap(const std::string& pName,
             double pZ,
             double pY,
             double pX, double pLTX)
  : VUSolid(pName)
{
  bool good;

  if (pZ <= 0 || pY <= 0 || pX <= 0 || pLTX <= 0 || pLTX > pX)
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  fDz = 0.5 * pZ ;
  fTthetaCphi = 0 ;
  fTthetaSphi = 0 ;

  fDy1 = 0.5 * pY;
  fDx1 = 0.5 * pX ;
  fDx2 = 0.5 * pLTX;
  fTalpha1 =  0.5 * (pLTX - pX) / pY;

  fDy2 = fDy1 ;
  fDx3 = fDx1;
  fDx4 = fDx2 ;
  fTalpha2 = fTalpha1 ;

  UVector3 pt[8] ;

  pt[0] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 - fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[1] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 + fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[2] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 - fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[3] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 + fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[4] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 - fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[5] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 + fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[6] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 - fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);
  pt[7] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 + fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);

  // Bottom side with normal approx. -Y
  //
  good = MakePlane(pt[0], pt[4], pt[5], pt[1], fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Top side with normal approx. +Y
  //
  good = MakePlane(pt[2], pt[3], pt[7], pt[6], fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Front side with normal approx. -X
  //
  good = MakePlane(pt[0], pt[2], pt[6], pt[4], fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Back side iwth normal approx. +X
  //
  good = MakePlane(pt[1], pt[5], pt[7], pt[3], fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }
  fCubicVolume=0.;
  fSurfaceArea=0.;
}

///////////////////////////////////////////////////////////////////////////////
//
// Constructor for UTrd

UTrap::UTrap(const std::string& pName,
             double pDx1,  double pDx2,
             double pDy1,  double pDy2,
             double pDz)
  : VUSolid(pName)
{
  bool good;

  if (pDz <= 0 || pDy1 <= 0 || pDx1 <= 0 || pDx2 <= 0 || pDy2 <= 0)
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  fDz = pDz;
  fTthetaCphi = 0 ;
  fTthetaSphi = 0 ;

  fDy1 = pDy1 ;
  fDx1 = pDx1 ;
  fDx2 = pDx1 ;
  fTalpha1 = 0 ;

  fDy2 = pDy2 ;
  fDx3 = pDx2 ;
  fDx4 = pDx2 ;
  fTalpha2 = 0 ;

  UVector3 pt[8] ;

  pt[0] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 - fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[1] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 + fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[2] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 - fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[3] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 + fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[4] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 - fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[5] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 + fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[6] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 - fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);
  pt[7] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 + fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);

  // Bottom side with normal approx. -Y
  //
  good = MakePlane(pt[0], pt[4], pt[5], pt[1], fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Top side with normal approx. +Y
  //
  good = MakePlane(pt[2], pt[3], pt[7], pt[6], fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Front side with normal approx. -X
  //
  good = MakePlane(pt[0], pt[2], pt[6], pt[4], fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Back side iwth normal approx. +X
  //
  good = MakePlane(pt[1], pt[5], pt[7], pt[3], fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume=0.;
  fSurfaceArea=0.;
}

////////////////////////////////////////////////////////////////////////////
//
// Constructor for UPara

UTrap::UTrap(const std::string& pName,
             double pDx, double pDy,
             double pDz,
             double pAlpha,
             double pTheta, double pPhi)
  : VUSolid(pName)
{
  bool good;

  if (pDz <= 0 || pDy <= 0 || pDx <= 0)
  {
    std::ostringstream message;
    message << "Invalid length parameters for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  fDz = pDz ;
  fTthetaCphi = std::tan(pTheta) * std::cos(pPhi) ;
  fTthetaSphi = std::tan(pTheta) * std::sin(pPhi) ;

  fDy1 = pDy ;
  fDx1 = pDx ;
  fDx2 = pDx ;
  fTalpha1 = std::tan(pAlpha) ;

  fDy2 = pDy ;
  fDx3 = pDx ;
  fDx4 = pDx ;
  fTalpha2 = fTalpha1 ;

  UVector3 pt[8] ;

  pt[0] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 - fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[1] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 + fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[2] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 - fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[3] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 + fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[4] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 - fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[5] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 + fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[6] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 - fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);
  pt[7] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 + fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);

  // Bottom side with normal approx. -Y
  //
  good = MakePlane(pt[0], pt[4], pt[5], pt[1], fPlanes[0]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Top side with normal approx. +Y
  //
  good = MakePlane(pt[2], pt[3], pt[7], pt[6], fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Front side with normal approx. -X
  //
  good = MakePlane(pt[0], pt[2], pt[6], pt[4], fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Back side iwth normal approx. +X
  //
  good = MakePlane(pt[1], pt[5], pt[7], pt[3], fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume=0.;
  fSurfaceArea=0.;
}

///////////////////////////////////////////////////////////////////////////
//
// Nominal constructor for UTrap whose parameters are to be Set by
// a UVParamaterisation later.  Check and Set half-widths as well as
// angles: final check of coplanarity

UTrap::UTrap(const std::string& pName)
  : VUSolid(pName), fDz(1.), fTthetaCphi(0.), fTthetaSphi(0.),
    fDy1(1.), fDx1(1.), fDx2(1.), fTalpha1(0.),
    fDy2(1.), fDx3(1.), fDx4(1.), fTalpha2(0.)
{
  MakePlanes();
  fCubicVolume=0.;
  fSurfaceArea=0.;
}


////////////////////////////////////////////////////////////////////////
//
// Destructor

UTrap::~UTrap()
{
}

//////////////////////////////////////////////////////////////////////////
//
// Copy constructor

UTrap::UTrap(const UTrap& rhs)
  : VUSolid(rhs), fDz(rhs.fDz),
    fTthetaCphi(rhs.fTthetaCphi), fTthetaSphi(rhs.fTthetaSphi),
    fDy1(rhs.fDy1), fDx1(rhs.fDx1), fDx2(rhs.fDx2), fTalpha1(rhs.fTalpha1),
    fDy2(rhs.fDy2), fDx3(rhs.fDx3), fDx4(rhs.fDx4), fTalpha2(rhs.fTalpha2)
{
  for (size_t i = 0; i < 4; ++i)
  {
    fPlanes[i].a = rhs.fPlanes[i].a;
    fPlanes[i].b = rhs.fPlanes[i].b;
    fPlanes[i].c = rhs.fPlanes[i].c;
    fPlanes[i].d = rhs.fPlanes[i].d;
  }
  fCubicVolume=rhs.fCubicVolume;
  fSurfaceArea=rhs.fSurfaceArea;
}

//////////////////////////////////////////////////////////////////////////
//
// Assignment operator

UTrap& UTrap::operator = (const UTrap& rhs)
{
  // Check assignment to self
  //
  if (this == &rhs)
  {
    return *this;
  }

  // Copy base class data
  //
  VUSolid::operator=(rhs);

  // Copy data
  //
  fDz = rhs.fDz;
  fTthetaCphi = rhs.fTthetaCphi;
  fTthetaSphi = rhs.fTthetaSphi;
  fDy1 = rhs.fDy1;
  fDx1 = rhs.fDx1;
  fDx2 = rhs.fDx2;
  fTalpha1 = rhs.fTalpha1;
  fDy2 = rhs.fDy2;
  fDx3 = rhs.fDx3;
  fDx4 = rhs.fDx4;
  fTalpha2 = rhs.fTalpha2;
  for (size_t i = 0; i < 4; ++i)
  {
    fPlanes[i].a = rhs.fPlanes[i].a;
    fPlanes[i].b = rhs.fPlanes[i].b;
    fPlanes[i].c = rhs.fPlanes[i].c;
    fPlanes[i].d = rhs.fPlanes[i].d;
  }
  fCubicVolume=rhs.fCubicVolume;
  fSurfaceArea=rhs.fSurfaceArea;
  return *this;
}

///////////////////////////////////////////////////////////////////////
//
// Set all parameters, as for constructor - check and Set half-widths
// as well as angles: final check of coplanarity

void UTrap::SetAllParameters(double pDz,
                             double pTheta,
                             double pPhi,
                             double pDy1,
                             double pDx1,
                             double pDx2,
                             double pAlp1,
                             double pDy2,
                             double pDx3,
                             double pDx4,
                             double pAlp2)
{
  if (pDz <= 0 || pDy1 <= 0 || pDx1 <= 0 || pDx2 <= 0 || pDy2 <= 0 || pDx3 <= 0 || pDx4 <= 0)
  {
    std::ostringstream message;
    message << "Invalid Length Parameters for Solid: " << GetName() << std::endl
            << "				X - "
            << pDx1 << ", " << pDx2 << ", " << pDx3 << ", " << pDx4 << std::endl
            << "					Y - " << pDy1 << ", " << pDy2 << std::endl
            << "					Z - " << pDz;
    UUtils::Exception("UTrap::SetAllParameters()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fCubicVolume = 0.;
  fSurfaceArea = 0.;
  fDz = pDz;
  fTthetaCphi = std::tan(pTheta) * std::cos(pPhi);
  fTthetaSphi = std::tan(pTheta) * std::sin(pPhi);

  fDy1 = pDy1;
  fDx1 = pDx1;
  fDx2 = pDx2;
  fTalpha1 = std::tan(pAlp1);

  fDy2 = pDy2;
  fDx3 = pDx3;
  fDx4 = pDx4;
  fTalpha2 = std::tan(pAlp2);

  MakePlanes();
}

void UTrap::SetPlanes(const UVector3 pt[8])
{
  bool good;

  // Start with check of centering - the center of gravity trap line
  // should Cross the origin of frame

  if (!(pt[0].z() < 0
        && pt[0].z() == pt[1].z() && pt[0].z() == pt[2].z()
        && pt[0].z() == pt[3].z()
        && pt[4].z() > 0
        && pt[4].z() == pt[5].z() && pt[4].z() == pt[6].z()
        && pt[4].z() == pt[7].z()
        && std::fabs(pt[0].z() + pt[4].z()) < VUSolid::Tolerance()
        && pt[0].y() == pt[1].y() && pt[2].y() == pt[3].y()
        && pt[4].y() == pt[5].y() && pt[6].y() == pt[7].y()
        && std::fabs(pt[0].y() + pt[2].y() + pt[4].y() + pt[6].y()) < VUSolid::Tolerance()
        && std::fabs(pt[0].x() + pt[1].x() + pt[4].x() + pt[5].x() +
                     pt[2].x() + pt[3].x() + pt[6].x() + pt[7].x()) < VUSolid::Tolerance()))
  {
    std::ostringstream message;
    message << "Invalid vertice coordinates for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Bottom side with normal approx. -Y

  good = MakePlane(pt[0], pt[4], pt[5], pt[1], fPlanes[0]);

  if (!good)
  {

    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002", UFatalError, 1,
                      "Face at ~-Y not planar.");
  }

  // Top side with normal approx. +Y

  good = MakePlane(pt[2], pt[3], pt[7], pt[6], fPlanes[1]);

  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Front side with normal approx. -X

  good = MakePlane(pt[0], pt[2], pt[6], pt[4], fPlanes[2]);

  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  // Back side iwth normal approx. +X

  good = MakePlane(pt[1], pt[5], pt[7], pt[3], fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::UTrap()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }
  fDz = (pt[7]).z() ;

  fDy1     = ((pt[2]).y() - (pt[1]).y()) * 0.5;
  fDx1     = ((pt[1]).x() - (pt[0]).x()) * 0.5;
  fDx2     = ((pt[3]).x() - (pt[2]).x()) * 0.5;
  fTalpha1 = ((pt[2]).x() + (pt[3]).x() - (pt[1]).x() - (pt[0]).x()) * 0.25 / fDy1;

  fDy2     = ((pt[6]).y() - (pt[5]).y()) * 0.5;
  fDx3     = ((pt[5]).x() - (pt[4]).x()) * 0.5;
  fDx4     = ((pt[7]).x() - (pt[6]).x()) * 0.5;
  fTalpha2 = ((pt[6]).x() + (pt[7]).x() - (pt[5]).x() - (pt[4]).x()) * 0.25 / fDy2;

  fTthetaCphi = ((pt[4]).x() + fDy2 * fTalpha2 + fDx3) / fDz;
  fTthetaSphi = ((pt[4]).y() + fDy2) / fDz;
}

//////////////////////////////////////////////////////////////////////////
//
// Checking of coplanarity

bool UTrap::MakePlanes()
{
  bool good = true;

  UVector3 pt[8] ;

  pt[0] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 - fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[1] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 + fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[2] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 - fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[3] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 + fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[4] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 - fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[5] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 + fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[6] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 - fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);
  pt[7] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 + fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);

  // Bottom side with normal approx. -Y
  //
  good = MakePlane(pt[0], pt[4], pt[5], pt[1], fPlanes[0]) ;
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::MakePlanes()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  // Top side with normal approx. +Y
  //
  good = MakePlane(pt[2], pt[3], pt[7], pt[6], fPlanes[1]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+Y not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::MakePlanes()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  // Front side with normal approx. -X
  //
  good = MakePlane(pt[0], pt[2], pt[6], pt[4], fPlanes[2]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~-X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::MakePlanes()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  // Back side iwth normal approx. +X
  //
  good = MakePlane(pt[1], pt[5], pt[7], pt[3], fPlanes[3]);
  if (!good)
  {
    std::ostringstream message;
    message << "Face at ~+X not planar for Solid: " << GetName();
    UUtils::Exception("UTrap::MakePlanes()", "GeomSolids0002",
                      UFatalError, 1, message.str().c_str());
  }

  return good;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate the coef's of the plane p1->p2->p3->p4->p1
// where the ThreeVectors 1-4 are in anti-clockwise order when viewed from
// infront of the plane (i.e. from normal direction).
//
// Return true if the ThreeVectors are coplanar + Set coef;s
//        false if ThreeVectors are not coplanar

bool UTrap::MakePlane(const UVector3& p1,
                      const UVector3& p2,
                      const UVector3& p3,
                      const UVector3& p4,
                      UTrapSidePlane& plane)
{
  double a, b, c, sd;
  UVector3 v12, v13, v14, Vcross;

  bool good;

  v12   = p2 - p1;
  v13   = p3 - p1;
  v14   = p4 - p1;
  Vcross = v12.Cross(v13);

  if (std::fabs(Vcross.Dot(v14) / (Vcross.Mag()*v14.Mag())) > kCoplanar_Tolerance)
  {
    good = false;
  }
  else
  {
    // a,b,c correspond to the x/y/z components of the
    // normal vector to the plane

    // a  = (p2.y()-p1.y())*(p1.z()+p2.z())+(p3.y()-p2.y())*(p2.z()+p3.z());
    // a += (p4.y()-p3.y())*(p3.z()+p4.z())+(p1.y()-p4.y())*(p4.z()+p1.z()); // ?
    // b  = (p2.z()-p1.z())*(p1.x()+p2.x())+(p3.z()-p2.z())*(p2.x()+p3.x());
    // b += (p4.z()-p3.z())*(p3.x()+p4.x())+(p1.z()-p4.z())*(p4.x()+p1.x()); // ?
    // c  = (p2.x()-p1.x())*(p1.y()+p2.y())+(p3.x()-p2.x())*(p2.y()+p3.y());
    // c += (p4.x()-p3.x())*(p3.y()+p4.y())+(p1.x()-p4.x())*(p4.y()+p1.y()); // ?

    // Let create diagonals 4-2 and 3-1 than (4-2)x(3-1) provides
    // vector perpendicular to the plane directed to outside !!!
    // and a,b,c, = f(1,2,3,4) external relative to trap normal

    a = +(p4.y() - p2.y()) * (p3.z() - p1.z())
        - (p3.y() - p1.y()) * (p4.z() - p2.z());

    b = -(p4.x() - p2.x()) * (p3.z() - p1.z())
        + (p3.x() - p1.x()) * (p4.z() - p2.z());

    c = +(p4.x() - p2.x()) * (p3.y() - p1.y())
        - (p3.x() - p1.x()) * (p4.y() - p2.y());

    sd = std::sqrt(a * a + b * b + c * c); // so now vector plane.(a,b,c) is Unit

    if (sd > 0)
    {
      plane.a = a / sd;
      plane.b = b / sd;
      plane.c = c / sd;
    }
    else
    {
      std::ostringstream message;
      message << "Invalid parameters: norm.mod() <= 0, for Solid: "
              << GetName();
      UUtils::Exception("UTrap::MakePlanes()", "GeomSolids0002",
                        UFatalError, 1, message.str().c_str()) ;
    }
    // Calculate D: p1 in in plane so D=-n.p1.Vect()

    plane.d = -(plane.a * p1.x() + plane.b * p1.y() + plane.c * p1.z());

    good = true;
  }
  return good;
}




////////////////////////////////////////////////////////////////////////
//
// Return whether point inside/outside/on surface, using tolerance

VUSolid::EnumInside UTrap::Inside(const UVector3& p) const
{
  VUSolid::EnumInside in;
  double Dist;
  int i;
  if (std::fabs(p.z()) <= fDz - VUSolid::Tolerance() * 0.5)
  {
    in = eInside;

    for (i = 0; i < 4; i++)
    {
      Dist = fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
             + fPlanes[i].c * p.z() + fPlanes[i].d;

      if (Dist > VUSolid::Tolerance() * 0.5) return in = eOutside;
      else if (Dist > -VUSolid::Tolerance() * 0.5)  in = eSurface;

    }
  }
  else if (std::fabs(p.z()) <= fDz + VUSolid::Tolerance() * 0.5)
  {
    in = eSurface;

    for (i = 0; i < 4; i++)
    {
      Dist =  fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
              + fPlanes[i].c * p.z() + fPlanes[i].d;

      if (Dist > VUSolid::Tolerance() * 0.5) return in = eOutside;
    }
  }
  else  in = eOutside;

  return in;
}

/////////////////////////////////////////////////////////////////////////////
//
// Calculate side nearest to p, and return normal
// If 2+ sides equidistant, first side's normal returned (arbitrarily)

bool UTrap::Normal(const UVector3& p, UVector3& aNormal) const
{
  int i, noSurfaces = 0;
  double dist, distz, distx, disty, distmx, distmy, safe = UUtils::kInfinity;
  double delta    = 0.5 * VUSolid::Tolerance();
  UVector3 norm, sumnorm(0., 0., 0.);

  for (i = 0; i < 4; i++)
  {
    dist =  std::fabs(fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
                      + fPlanes[i].c * p.z() + fPlanes[i].d);
    if (dist < safe)
    {
      safe = dist;
    }
  }
  distz = std::fabs(std::fabs(p.z()) - fDz);

  distmy = std::fabs(fPlanes[0].a * p.x() + fPlanes[0].b * p.y()
                     + fPlanes[0].c * p.z() + fPlanes[0].d);

  disty = std::fabs(fPlanes[1].a * p.x() + fPlanes[1].b * p.y()
                    + fPlanes[1].c * p.z() + fPlanes[1].d);

  distmx = std::fabs(fPlanes[2].a * p.x() + fPlanes[2].b * p.y()
                     + fPlanes[2].c * p.z() + fPlanes[2].d);

  distx = std::fabs(fPlanes[3].a * p.x() + fPlanes[3].b * p.y()
                    + fPlanes[3].c * p.z() + fPlanes[3].d);

  UVector3 nX = UVector3(fPlanes[3].a, fPlanes[3].b, fPlanes[3].c);
  UVector3 nmX = UVector3(fPlanes[2].a, fPlanes[2].b, fPlanes[2].c);
  UVector3 nY = UVector3(fPlanes[1].a, fPlanes[1].b, fPlanes[1].c);
  UVector3 nmY = UVector3(fPlanes[0].a, fPlanes[0].b, fPlanes[0].c);
  UVector3 nZ = UVector3(0., 0., 1.0);

  if (distx <= delta)
  {
    noSurfaces ++;
    sumnorm += nX;
  }
  if (distmx <= delta)
  {
    noSurfaces ++;
    sumnorm += nmX;
  }
  if (disty <= delta)
  {
    noSurfaces ++;
    sumnorm += nY;
  }
  if (distmy <= delta)
  {
    noSurfaces ++;
    sumnorm += nmY;
  }
  if (distz <= delta)
  {
    noSurfaces ++;
    if (p.z() >= 0.)  sumnorm += nZ;
    else               sumnorm -= nZ;
  }
  if (noSurfaces == 0)
  {
#ifdef UDEBUG
    UUtils::Exception("UTrap::SurfaceNormal(p)", "GeomSolids1002",
                      UWarning, 1, "Point p is not on surface !?");
#endif
    norm = ApproxSurfaceNormal(p);
  }
  else if (noSurfaces == 1) norm = sumnorm;
  else                  norm = sumnorm.Unit();
  aNormal = norm;
  return noSurfaces != 0;
}

////////////////////////////////////////////////////////////////////////////////////
//
// Algorithm for SurfaceNormal() following the original specification
// for points not on the surface

UVector3 UTrap::ApproxSurfaceNormal(const UVector3& p) const
{
  double safe = UUtils::kInfinity, Dist, safez;
  int i, imin = 0;
  for (i = 0; i < 4; i++)
  {
    Dist = std::fabs(fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
                     + fPlanes[i].c * p.z() + fPlanes[i].d);
    if (Dist < safe)
    {
      safe = Dist;
      imin = i;
    }
  }
  safez = std::fabs(std::fabs(p.z()) - fDz);
  if (safe < safez)
  {
    return UVector3(fPlanes[imin].a, fPlanes[imin].b, fPlanes[imin].c);
  }
  else
  {
    if (p.z() > 0)
    {
      return UVector3(0, 0, 1);
    }
    else
    {
      return UVector3(0, 0, -1);
    }
  }
}

////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to shape from outside - return UUtils::kInfinity if no intersection
//
// ALGORITHM:
// For each component, calculate pair of minimum and maximum intersection
// values for which the particle is in the extent of the shape
// - The smallest (MAX minimum) allowed distance of the pairs is intersect

double UTrap::DistanceToIn(const UVector3& p,
                           const UVector3& v, double) const
{

  double snxt;    // snxt = default return value
  double max, smax, smin;
  double pdist, Comp, vdist;
  int i;
  //
  // Z Intersection range
  //
  if (v.z() > 0)
  {
    max = fDz - p.z() ;
    if (max > 0.5 * VUSolid::Tolerance())
    {
      smax = max / v.z();
      smin = (-fDz - p.z()) / v.z();
    }
    else
    {
      return snxt = UUtils::kInfinity;
    }
  }
  else if (v.z() < 0)
  {
    max = - fDz - p.z() ;
    if (max < -0.5 * VUSolid::Tolerance())
    {
      smax = max / v.z();
      smin = (fDz - p.z()) / v.z();
    }
    else
    {
      return snxt = UUtils::kInfinity;
    }
  }
  else
  {
    if (std::fabs(p.z()) < fDz - 0.5 * VUSolid::Tolerance()) // Inside was <=fDz
    {
      smin = 0;
      smax = UUtils::kInfinity;
    }
    else
    {
      return snxt = UUtils::kInfinity;
    }
  }

  for (i = 0; i < 4; i++)
  {
    pdist = fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
            + fPlanes[i].c * p.z() + fPlanes[i].d;
    Comp = fPlanes[i].a * v.x() + fPlanes[i].b * v.y() + fPlanes[i].c * v.z();
    if (pdist >= -0.5 * VUSolid::Tolerance())      // was >0
    {
      //
      // Outside the plane -> this is an extent entry distance
      //
      if (Comp >= 0)   // was >0
      {
        return snxt = UUtils::kInfinity ;
      }
      else
      {
        vdist = -pdist / Comp;
        if (vdist > smin)
        {
          if (vdist < smax)
          {
            smin = vdist;
          }
          else
          {
            return snxt = UUtils::kInfinity;
          }
        }
      }
    }
    else
    {
      //
      // Inside the plane -> couble be an extent exit distance (smax)
      //
      if (Comp > 0) // Will leave extent
      {
        vdist = -pdist / Comp;
        if (vdist < smax)
        {
          if (vdist > smin)
          {
            smax = vdist;
          }
          else
          {
            return snxt = UUtils::kInfinity;
          }
        }
      }
    }
  }
  //
  // Checks in non z plane intersections ensure smin<smax
  //
  if (smin >= 0)
  {
    snxt = smin ;
  }
  else
  {
    snxt = 0 ;
  }
  return snxt;
}

///////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from outside
// This is the best fast estimation of the shortest distance to trap
// - Returns 0 is ThreeVector inside

double UTrap::SafetyFromOutside(const UVector3& p, bool) const
{
  double safe = 0.0, Dist;
  int i;
  safe = std::fabs(p.z()) - fDz;
  for (i = 0; i < 4; i++)
  {
    Dist = fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
           + fPlanes[i].c * p.z() + fPlanes[i].d;
    if (Dist > safe) safe = Dist;
  }
  if (safe < 0) safe = 0;
  return safe;
}

/////////////////////////////////////////////////////////////////////////////////
//
// Calculate distance to surface of shape from inside
// Calculate distance to x/y/z planes - smallest is exiting distance

double UTrap::DistanceToOut(const UVector3& p, const UVector3&  v, UVector3&       aNormalVector, bool&           aConvex, double) const
{
  Eside side = kUndef;
  double snxt;    // snxt = return value
  double pdist, Comp, vdist, max;
  //
  // Z Intersections
  //
  if (v.z() > 0)
  {
    max = fDz - p.z();
    if (max > VUSolid::Tolerance() / 2)
    {
      snxt = max / v.z();
      side = kPZ;
    }
    else
    {
      aConvex = true;
      aNormalVector = UVector3(0, 0, 1);
      return snxt = 0;
    }
  }
  else if (v.z() < 0)
  {
    max = -fDz - p.z();
    if (max < -VUSolid::Tolerance() / 2)
    {
      snxt = max / v.z();
      side = kMZ;
    }
    else
    {
      aConvex = true;
      aNormalVector = UVector3(0, 0, -1);
      return snxt = 0;
    }
  }
  else
  {
    snxt = UUtils::kInfinity;
  }

  //
  // Intersections with planes[0] (expanded because of setting enum)
  //
  pdist = fPlanes[0].a * p.x() + fPlanes[0].b * p.y() + fPlanes[0].c * p.z() + fPlanes[0].d;
  Comp = fPlanes[0].a * v.x() + fPlanes[0].b * v.y() + fPlanes[0].c * v.z();
  if (pdist > 0)
  {
    // Outside the plane
    if (Comp > 0)
    {
      // Leaving immediately
      aConvex = true;
      aNormalVector = UVector3(fPlanes[0].a, fPlanes[0].b, fPlanes[0].c);
      return snxt = 0;
    }
  }
  else if (pdist < -VUSolid::Tolerance() / 2)
  {
    // Inside the plane
    if (Comp > 0)
    {
      // Will leave extent
      vdist = -pdist / Comp;
      if (vdist < snxt)
      {
        snxt = vdist;
        side = ks0;
      }
    }
  }
  else
  {
    // On surface
    if (Comp > 0)
    {
      aConvex = true;
      aNormalVector = UVector3(fPlanes[0].a, fPlanes[0].b, fPlanes[0].c);
      return snxt = 0;
    }
  }

  //
  // Intersections with planes[1] (expanded because of setting enum)
  //
  pdist = fPlanes[1].a * p.x() + fPlanes[1].b * p.y() + fPlanes[1].c * p.z() + fPlanes[1].d;
  Comp = fPlanes[1].a * v.x() + fPlanes[1].b * v.y() + fPlanes[1].c * v.z();
  if (pdist > 0)
  {
    // Outside the plane
    if (Comp > 0)
    {
      // Leaving immediately
      aConvex = true;
      aNormalVector = UVector3(fPlanes[1].a, fPlanes[1].b, fPlanes[1].c);
      return snxt = 0;
    }
  }
  else if (pdist < -VUSolid::Tolerance() / 2)
  {
    // Inside the plane
    if (Comp > 0)
    {
      // Will leave extent
      vdist = -pdist / Comp;
      if (vdist < snxt)
      {
        snxt = vdist;
        side = ks1;
      }
    }
  }
  else
  {
    // On surface
    if (Comp > 0)
    {
      aConvex = true;
      aNormalVector = UVector3(fPlanes[1].a, fPlanes[1].b, fPlanes[1].c);
      return snxt = 0;
    }
  }

  //
  // Intersections with planes[2] (expanded because of setting enum)
  //
  pdist = fPlanes[2].a * p.x() + fPlanes[2].b * p.y() + fPlanes[2].c * p.z() + fPlanes[2].d;
  Comp = fPlanes[2].a * v.x() + fPlanes[2].b * v.y() + fPlanes[2].c * v.z();
  if (pdist > 0)
  {
    // Outside the plane
    if (Comp > 0)
    {
      // Leaving immediately
      aConvex = true;
      aNormalVector = UVector3(fPlanes[2].a, fPlanes[2].b, fPlanes[2].c);
      return snxt = 0;
    }
  }
  else if (pdist < -VUSolid::Tolerance() / 2)
  {
    // Inside the plane
    if (Comp > 0)
    {
      // Will leave extent
      vdist = -pdist / Comp;
      if (vdist < snxt)
      {
        snxt = vdist;
        side = ks2;
      }
    }
  }
  else
  {
    // On surface
    if (Comp > 0)
    {
      aConvex = true;
      aNormalVector = UVector3(fPlanes[2].a, fPlanes[2].b, fPlanes[2].c);
      return snxt = 0;
    }
  }

  //
  // Intersections with planes[3] (expanded because of setting enum)
  //
  pdist = fPlanes[3].a * p.x() + fPlanes[3].b * p.y() + fPlanes[3].c * p.z() + fPlanes[3].d;
  Comp = fPlanes[3].a * v.x() + fPlanes[3].b * v.y() + fPlanes[3].c * v.z();
  if (pdist > 0)
  {
    // Outside the plane
    if (Comp > 0)
    {
      // Leaving immediately
      aConvex = true;
      aNormalVector = UVector3(fPlanes[3].a, fPlanes[3].b, fPlanes[3].c);
      return snxt = 0;
    }
  }
  else if (pdist < -VUSolid::Tolerance() / 2)
  {
    // Inside the plane
    if (Comp > 0)
    {
      // Will leave extent
      vdist = -pdist / Comp;
      if (vdist < snxt)
      {
        snxt = vdist;
        side = ks3;
      }
    }
  }
  else
  {
    // On surface
    if (Comp > 0)
    {
      aConvex = true;
      aNormalVector = UVector3(fPlanes[3].a, fPlanes[3].b, fPlanes[3].c);
      return snxt = 0;
    }
  }

  // Set normal
  aConvex = true;
  switch (side)
  {
    case ks0:
      aNormalVector = UVector3(fPlanes[0].a, fPlanes[0].b, fPlanes[0].c);
      break;
    case ks1:
      aNormalVector = UVector3(fPlanes[1].a, fPlanes[1].b, fPlanes[1].c);
      break;
    case ks2:
      aNormalVector = UVector3(fPlanes[2].a, fPlanes[2].b, fPlanes[2].c);
      break;
    case ks3:
      aNormalVector = UVector3(fPlanes[3].a, fPlanes[3].b, fPlanes[3].c);
      break;
    case kMZ:
      aNormalVector = UVector3(0, 0, -1);
      break;
    case kPZ:
      aNormalVector = UVector3(0, 0, 1);
      break;
    default:
      cout << std::endl;

      std::ostringstream message;
      int oldprc = message.precision(16);
      message << "Undefined side for valid surface normal to solid."
              << std::endl
              << "Position:"  << std::endl << std::endl
              << "p.x = "  << p.x() << " mm" << std::endl
              << "p.y = "  << p.y() << " mm" << std::endl
              << "p.z = "  << p.z() << " mm" << std::endl << std::endl
              << "Direction:" << std::endl << std::endl
              << "v.x = "  << v.x() << std::endl
              << "v.y = "  << v.y() << std::endl
              << "v.z = "  << v.z() << std::endl << std::endl
              << "Proposed distance :" << std::endl << std::endl
              << "snxt = "    << snxt << " mm" << std::endl;
      message.precision(oldprc);
      UUtils::Exception("UTrap::DistanceToOut(p,v,..)", "GeomSolids1002",
                        UWarning, 1, message.str().c_str());
      break;
  }
  return snxt;
}

//////////////////////////////////////////////////////////////////////////////
//
// Calculate exact shortest distance to any boundary from inside
// - Returns 0 is ThreeVector outside

double UTrap::SafetyFromInside(const UVector3& p, bool /*precise*/) const
{
  double safe = 0.0, Dist;
  int i;

#ifdef UDEBUG
  if (Inside(p) == eOutside)
  {
    int oldprc = cout.precision(16) ;
    cout << std::endl ;

    cout << "Position:"  << std::endl << std::endl ;
    cout << "p.x = "  << p.x() << " mm" << std::endl ;
    cout << "p.y = "  << p.y() << " mm" << std::endl ;
    cout << "p.z = "  << p.z() << " mm" << std::endl << std::endl ;
    cout.precision(oldprc) ;
    UUtils::Exception("UTrap::DistanceToOut(p)",
                      "GeomSolids1002", UWarning, 1, "Point p is outside !?");
  }
#endif

  safe = fDz - std::fabs(p.z());
  if (safe < 0) safe = 0;
  else
  {
    for (i = 0; i < 4; i++)
    {
      Dist = -(fPlanes[i].a * p.x() + fPlanes[i].b * p.y()
               + fPlanes[i].c * p.z() + fPlanes[i].d);
      if (Dist < safe) safe = Dist;
    }
    if (safe < 0) safe = 0;
  }
  return safe;
}


//////////////////////////////////////////////////////////////////////////
//
// GetEntityType

UGeometryType UTrap::GetEntityType() const
{
  return std::string("Trap");
}

//////////////////////////////////////////////////////////////////////////
//
// Make a clone of the object
//
VUSolid* UTrap::Clone() const
{
  return new UTrap(*this);
}

//////////////////////////////////////////////////////////////////////////
//
// Stream object contents to an output stream

std::ostream& UTrap::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "		*** Dump for solid - " << GetName() << " ***\n"
     << "		===================================================\n"
     << " Solid type: UTrap\n"
     << " Parameters: \n"
     << "		half length Z: " << fDz << " mm \n"
     << "		half length Y of face -fDz: " << fDy1 << " mm \n"
     << "		half length X of side -fDy1, face -fDz: " << fDx1 << " mm \n"
     << "		half length X of side +fDy1, face -fDz: " << fDx2 << " mm \n"
     << "		half length Y of face +fDz: " << fDy2 << " mm \n"
     << "		half length X of side -fDy2, face +fDz: " << fDx3 << " mm \n"
     << "		half length X of side +fDy2, face +fDz: " << fDx4 << " mm \n"
     << "		std::tan(theta)*std::cos(phi): " << fTthetaCphi / (UUtils::kPi / 180.0) << " degrees \n"
     << "		std::tan(theta)*std::sin(phi): " << fTthetaSphi / (UUtils::kPi / 180.0) << " degrees \n"
     << "		std::tan(alpha), -fDz: " << fTalpha1 / (UUtils::kPi / 180.0) << " degrees \n"
     << "		std::tan(alpha), +fDz: " << fTalpha2 / (UUtils::kPi / 180.0) << " degrees \n"
     << "		trap side plane equations:\n"
     << "				" << fPlanes[0].a << " X + " << fPlanes[0].b << " Y + "
     << fPlanes[0].c << " Z + " << fPlanes[0].d << " = 0\n"
     << "				" << fPlanes[1].a << " X + " << fPlanes[1].b << " Y + "
     << fPlanes[1].c << " Z + " << fPlanes[1].d << " = 0\n"
     << "				" << fPlanes[2].a << " X + " << fPlanes[2].b << " Y + "
     << fPlanes[2].c << " Z + " << fPlanes[2].d << " = 0\n"
     << "				" << fPlanes[3].a << " X + " << fPlanes[3].b << " Y + "
     << fPlanes[3].c << " Z + " << fPlanes[3].d << " = 0\n"
     << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

/////////////////////////////////////////////////////////////////////////
//
// GetPointOnPlane
//
// Auxiliary method for Get Point on Surface

UVector3 UTrap::GetPointOnPlane(UVector3 p0, UVector3 p1,
                                UVector3 p2, UVector3 p3,
                                double& area) const
{
  double lambda1, lambda2, chose, aOne, aTwo;
  UVector3 t, u, v, w, Area, normal;

  t = p1 - p0;
  u = p2 - p1;
  v = p3 - p2;
  w = p0 - p3;

  Area = UVector3(w.y() * v.z() - w.z() * v.y(),
                  w.z() * v.x() - w.x() * v.z(),
                  w.x() * v.y() - w.y() * v.x());

  aOne = 0.5 * Area.Mag();

  Area = UVector3(t.y() * u.z() - t.z() * u.y(),
                  t.z() * u.x() - t.x() * u.z(),
                  t.x() * u.y() - t.y() * u.x());

  aTwo = 0.5 * Area.Mag();

  area = aOne + aTwo;

  chose = UUtils::Random(0., aOne + aTwo);

  if ((chose >= 0.) && (chose < aOne))
  {
    lambda1 = UUtils::Random(0., 1.);
    lambda2 = UUtils::Random(0., lambda1);
    return (p2 + lambda1 * v + lambda2 * w);
  }

  // else

  lambda1 = UUtils::Random(0., 1.);
  lambda2 = UUtils::Random(0., lambda1);

  return (p0 + lambda1 * t + lambda2 * u);
}

///////////////////////////////////////////////////////////////
//
// GetPointOnSurface

UVector3 UTrap::GetPointOnSurface() const
{
  double aOne, aTwo, aThree, aFour, aFive, aSix, chose;
  UVector3 One, Two, Three, Four, Five, Six, test;
  UVector3 pt[8];

  pt[0] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 - fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[1] = UVector3(-fDz * fTthetaCphi - fDy1 * fTalpha1 + fDx1,
                   -fDz * fTthetaSphi - fDy1, -fDz);
  pt[2] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 - fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[3] = UVector3(-fDz * fTthetaCphi + fDy1 * fTalpha1 + fDx2,
                   -fDz * fTthetaSphi + fDy1, -fDz);
  pt[4] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 - fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[5] = UVector3(+fDz * fTthetaCphi - fDy2 * fTalpha2 + fDx3,
                   +fDz * fTthetaSphi - fDy2, +fDz);
  pt[6] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 - fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);
  pt[7] = UVector3(+fDz * fTthetaCphi + fDy2 * fTalpha2 + fDx4,
                   +fDz * fTthetaSphi + fDy2, +fDz);

  // make sure we provide the points in a clockwise fashion

  One  = GetPointOnPlane(pt[0], pt[1], pt[3], pt[2], aOne);
  Two  = GetPointOnPlane(pt[4], pt[5], pt[7], pt[6], aTwo);
  Three = GetPointOnPlane(pt[6], pt[7], pt[3], pt[2], aThree);
  Four  = GetPointOnPlane(pt[4], pt[5], pt[1], pt[0], aFour);
  Five  = GetPointOnPlane(pt[0], pt[2], pt[6], pt[4], aFive);
  Six  = GetPointOnPlane(pt[1], pt[3], pt[7], pt[5], aSix);

  chose = UUtils::Random(0., aOne + aTwo + aThree + aFour + aFive + aSix);
  if ((chose >= 0.) && (chose < aOne))
  {
    return One;
  }
  else if ((chose >= aOne) && (chose < aOne + aTwo))
  {
    return Two;
  }
  else if ((chose >= aOne + aTwo) && (chose < aOne + aTwo + aThree))
  {
    return Three;
  }
  else if ((chose >= aOne + aTwo + aThree) && (chose < aOne + aTwo + aThree + aFour))
  {
    return Four;
  }
  else if ((chose >= aOne + aTwo + aThree + aFour)
           && (chose < aOne + aTwo + aThree + aFour + aFive))
  {
    return Five;
  }
  return Six;
}

void UTrap::Extent(UVector3& aMin, UVector3& aMax) const
{
 //Z axis
  aMin.z() = -fDz;
  aMax.z() = fDz;
  
  double temp[8] ;     // some points for intersection with zMin/zMax
  UVector3 pt[8];   // vertices after translation
    
  //X axis
  pt[0]=UVector3(-fDz*fTthetaCphi-fDy1*fTalpha1-fDx1,
                        -fDz*fTthetaSphi-fDy1,-fDz);
  pt[1]=UVector3(-fDz*fTthetaCphi-fDy1*fTalpha1+fDx1,
                       -fDz*fTthetaSphi-fDy1,-fDz);
  pt[2]=UVector3(-fDz*fTthetaCphi+fDy1*fTalpha1-fDx2,
                        -fDz*fTthetaSphi+fDy1,-fDz);
  pt[3]=UVector3(-fDz*fTthetaCphi+fDy1*fTalpha1+fDx2,
                        -fDz*fTthetaSphi+fDy1,-fDz);
  pt[4]=UVector3(+fDz*fTthetaCphi-fDy2*fTalpha2-fDx3,
                        fDz*fTthetaSphi-fDy2,+fDz);
  pt[5]=UVector3(fDz*fTthetaCphi-fDy2*fTalpha2+fDx3,
                        fDz*fTthetaSphi-fDy2,+fDz);
  pt[6]=UVector3(fDz*fTthetaCphi+fDy2*fTalpha2-fDx4,
                        fDz*fTthetaSphi+fDy2,+fDz);
  pt[7]=UVector3(fDz*fTthetaCphi+fDy2*fTalpha2+fDx4,
                        fDz*fTthetaSphi+fDy2,+fDz);

  temp[0] = pt[0].x()+(pt[4].x()-pt[0].x())
      *(aMin.z()-pt[0].z())/(pt[4].z()-pt[0].z()) ;
  temp[1] = pt[0].x()+(pt[4].x()-pt[0].x())
      *(aMax.z()-pt[0].z())/(pt[4].z()-pt[0].z()) ;
  temp[2] = pt[2].x()+(pt[6].x()-pt[2].x())
      *(aMin.z()-pt[2].z())/(pt[6].z()-pt[2].z()) ;
  temp[3] = pt[2].x()+(pt[6].x()-pt[2].x())
      *(aMax.z()-pt[2].z())/(pt[6].z()-pt[2].z()) ;
  temp[4] = pt[3].x()+(pt[7].x()-pt[3].x())
      *(aMin.z()-pt[3].z())/(pt[7].z()-pt[3].z()) ;
  temp[5] = pt[3].x()+(pt[7].x()-pt[3].x())
      *(aMax.z()-pt[3].z())/(pt[7].z()-pt[3].z()) ;
  temp[6] = pt[1].x()+(pt[5].x()-pt[1].x())
      *(aMin.z()-pt[1].z())/(pt[5].z()-pt[1].z()) ;
  temp[7] = pt[1].x()+(pt[5].x()-pt[1].x())
      *(aMax.z()-pt[1].z())/(pt[5].z()-pt[1].z()) ;
      
  aMax.x() =  - std::fabs(fDz*fTthetaCphi) - fDx1 - fDx2 -fDx3 - fDx4 ;
  aMin.x() = -aMax.x() ;

  for(int i = 0 ; i < 8 ; i++ )
  {
    if( temp[i] > aMax.x()) aMax.x() = temp[i] ;
    if( temp[i] < aMin.x()) aMin.x() = temp[i] ;
  }                                            
  //Y axis
  temp[0] = pt[0].y()+(pt[4].y()-pt[0].y())*(aMin.z()-pt[0].z())
                       /(pt[4].z()-pt[0].z()) ;
  temp[1] = pt[0].y()+(pt[4].y()-pt[0].y())*(aMax.z()-pt[0].z())
                       /(pt[4].z()-pt[0].z()) ;
  temp[2] = pt[2].y()+(pt[6].y()-pt[2].y())*(aMin.z()-pt[2].z())
                       /(pt[6].z()-pt[2].z()) ;
  temp[3] = pt[2].y()+(pt[6].y()-pt[2].y())*(aMax.z()-pt[2].z())
                       /(pt[6].z()-pt[2].z()) ;

  aMax.y() = - std::fabs(fDz*fTthetaSphi) - fDy1 - fDy2 ;
  aMin.y() = -aMax.y() ;
  
  for( int i = 0 ; i < 4 ; i++ )
  {
    if( temp[i] > aMax.y() ) aMax.y() = temp[i] ;
    if( temp[i] < aMin.y() ) aMin.y() = temp[i] ;
  }
}
