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
// UPolyhedra
//
// --------------------------------------------------------------------
//
// To be done:
//    * Cracks: there are probably small cracks in the seams between the
//      phi face (UPolyPhiFace) and sides (UPolyhedraSide) that are not
//      entirely leakproof. Also, I am not sure all vertices are leak proof.
//    * Many optimizations are possible, but not implemented.
//    * Visualization needs to be updated outside of this routine.
//
// Utility classes:
//    * UEnclosingCylinder: I decided a quick check of geometry would be a
//      good idea (for CPU speed). If the quick check fails, the regular
//      full-blown UVCSGfaceted version is invoked.
//    * UReduciblePolygon: Really meant as a check of input parameters,
//      this utility class also "converts" the GEANT3-like PGON/PCON
//      arguments into the newer ones.
// Both these classes are implemented outside this file because they are
// shared with UPolycone.
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>

#include "UPolyhedra.hh"
#include "UPolyhedraSide.hh"
#include "UPolyPhiFace.hh"
#include "UEnclosingCylinder.hh"
#include "UReduciblePolygon.hh"

using namespace std;

UPolyhedra::UPolyhedra(const std::string& name,
                       double phiStart,
                       double thePhiTotal,
                       int thefNumSide,
                       int numZPlanes,
                       const double zPlane[],
                       const double rInner[],
                       const double rOuter[])
  : UVCSGfaceted(name)
{
  Init(phiStart, thePhiTotal, thefNumSide, numZPlanes, zPlane, rInner, rOuter);
}

//
// Constructor (GEANT3 style parameters)
//
// GEANT3 PGON radii are specified in the distance to the norm of each face.
//
void UPolyhedra::Init(
  double phiStart,
  double thePhiTotal,
  int thefNumSide,
  int numZPlanes,
  const double zPlane[],
  const double rInner[],
  const double rOuter[])
{
  fGenericPgon = false;

  if (thefNumSide <= 0)
  {
    std::ostringstream message;
    message << "Solid must have at least one side - " << GetName() << std::endl
            << "        No sides specified !";
    UUtils::Exception("UPolyhedra::UPolyhedra()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  //
  // Calculate conversion factor from G3 radius to U radius
  //
  double phiTotal = thePhiTotal;
  if ((phiTotal <= 0) || (phiTotal >= 2 * UUtils::kPi * (1 - DBL_EPSILON)))
  {
    phiTotal = 2 * UUtils::kPi;
  }
  double convertRad = std::cos(0.5 * phiTotal / thefNumSide);

  //
  // Some historical stuff
  //
//  fOriginalParameters = new UPolyhedraHistorical;

  fOriginalParameters.fNumSide = thefNumSide;
  fOriginalParameters.fStartAngle = phiStart;
  fOriginalParameters.fOpeningAngle = phiTotal;
  fOriginalParameters.fNumZPlanes = numZPlanes;
  fOriginalParameters.fZValues.resize(numZPlanes);
  fOriginalParameters.Rmin.resize(numZPlanes);
  fOriginalParameters.Rmax.resize(numZPlanes);

  int i;
  for (i = 0; i < numZPlanes; i++)
  {
    if ((i < numZPlanes - 1) && (zPlane[i] == zPlane[i + 1]))
    {
      if ((rInner[i]   > rOuter[i + 1])
          || (rInner[i + 1] > rOuter[i]))
      {

        std::ostringstream message;
        message << "Cannot create a Polyhedra with no contiguous segments."
                << std::endl
                << "        Segments are not contiguous !" << std::endl
                << "        rMin[" << i << "] = " << rInner[i]
                << " -- rMax[" << i + 1 << "] = " << rOuter[i + 1] << std::endl
                << "        rMin[" << i + 1 << "] = " << rInner[i + 1]
                << " -- rMax[" << i << "] = " << rOuter[i];
        UUtils::Exception("UPolyhedra::UPolyhedra()", "GeomSolids0002",
                          UFatalErrorInArguments, 1, message.str().c_str());
      }
    }
    fOriginalParameters.fZValues[i] = zPlane[i];
    fOriginalParameters.Rmin[i] = rInner[i] / convertRad;
    fOriginalParameters.Rmax[i] = rOuter[i] / convertRad;
  }


  //
  // Build RZ polygon using special PCON/PGON GEANT3 constructor
  //
  UReduciblePolygon* rz =
    new UReduciblePolygon(rInner, rOuter, zPlane, numZPlanes);
  rz->ScaleA(1 / convertRad);

  //
  // Do the real work
  //
  Create(phiStart, phiTotal, thefNumSide, rz);

  delete rz;
}


//
// Constructor (generic parameters)
//
UPolyhedra::UPolyhedra(const std::string& name,
                       double phiStart,
                       double phiTotal,
                       int    thefNumSide,
                       int    numRZ,
                       const double r[],
                       const double z[])
  : UVCSGfaceted(name), fGenericPgon(true)
{
  UReduciblePolygon* rz = new UReduciblePolygon(r, z, numRZ);

  Create(phiStart, phiTotal, thefNumSide, rz);

  // Set fOriginalParameters struct for consistency
  //
  SetOriginalParameters();

  delete rz;
}


//
// Create
//
// Generic create routine, called by each constructor
// after conversion of arguments
//
void UPolyhedra::Create(double phiStart,
                        double phiTotal,
                        int    thefNumSide,
                        UReduciblePolygon* rz)
{
  //
  // Perform checks of rz values
  //
  if (rz->Amin() < 0.0)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "        All R values must be >= 0 !";
    UUtils::Exception("UPolyhedra::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  double rzArea = rz->Area();
  if (rzArea < -VUSolid::Tolerance())
    rz->ReverseOrder();

  else if (rzArea < -VUSolid::Tolerance())
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "        R/Z Cross section is zero or near zero: " << rzArea;
    UUtils::Exception("UPolyhedra::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  if ((!rz->RemoveDuplicateVertices(VUSolid::Tolerance()))
      || (!rz->RemoveRedundantVertices(VUSolid::Tolerance())))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "        Too few unique R/Z values !";
    UUtils::Exception("UPolyhedra::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  if (rz->CrossesItself(1 / UUtils::kInfinity))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "        R/Z segments Cross !";
    UUtils::Exception("UPolyhedra::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  fNumCorner = rz->NumVertices();

  fStartPhi = phiStart;
  while (fStartPhi < 0) fStartPhi += 2 * UUtils::kPi;
  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if ((phiTotal <= 0) || (phiTotal > 2 * UUtils::kPi * (1 - DBL_EPSILON)))
  {
    fPhiIsOpen = false;
    fEndPhi = phiStart + 2 * UUtils::kPi;
  }
  else
  {
    fPhiIsOpen = true;

    //
    // Convert phi into our convention
    //
    fEndPhi = phiStart + phiTotal;
    while (fEndPhi < fStartPhi) fEndPhi += 2 * UUtils::kPi;
  }

  //
  // Save number sides
  //
  fNumSides = thefNumSide;

  //
  // Allocate corner array.
  //
  fCorners = new UPolyhedraSideRZ[fNumCorner];

  //
  // Copy fCorners
  //
  UReduciblePolygonIterator iterRZ(rz);

  UPolyhedraSideRZ* next = fCorners;
  iterRZ.Begin();
  do
  {
    next->r = iterRZ.GetA();
    next->z = iterRZ.GetB();
  }
  while (++next, iterRZ.Next());

  //
  // Allocate face pointer array
  //
  numFace = fPhiIsOpen ? fNumCorner + 2 : fNumCorner;
  faces = new UVCSGface*[numFace];

  //
  // Construct side faces
  //
  // To do so properly, we need to keep track of four successive RZ
  // fCorners.
  //
  // But! Don't construct a face if both points are at zero radius!

  //
  UPolyhedraSideRZ* corner = fCorners,
                    *prev = fCorners + fNumCorner - 1,
                     *nextNext;
  UVCSGface**   face = faces;
  do
  {
    next = corner + 1;
    if (next >= fCorners + fNumCorner) next = fCorners;
    nextNext = next + 1;
    if (nextNext >= fCorners + fNumCorner) nextNext = fCorners;

    if (corner->r < 1 / UUtils::kInfinity && next->r < 1 / UUtils::kInfinity) continue;
    /*
        // We must decide here if we can dare declare one of our faces
        // as having a "valid" normal (i.e. allBehind = true). This
        // is never possible if the face faces "inward" in r *unless*
        // we have only one side
        //
        bool allBehind;
        if ((corner->z > next->z) && (fNumSides > 1))
        {
          allBehind = false;
        }
        else
        {
          //
          // Otherwise, it is only true if the line passing
          // through the two points of the segment do not
          // split the r/z Cross section
          //
          allBehind = !rz->BisectedBy( corner->r, corner->z,
                                       next->r, next->z, VUSolid::Tolerance() );
        }
    */
    *face++ = new UPolyhedraSide(prev, corner, next, nextNext,
                                 fNumSides, fStartPhi, fEndPhi - fStartPhi, fPhiIsOpen);
  }
  while (prev = corner, corner = next, corner > fCorners);

  if (fPhiIsOpen)
  {
    //
    // Construct phi open edges
    //
    *face++ = new UPolyPhiFace(rz, fStartPhi, phiTotal / fNumSides, fEndPhi);
    *face++ = new UPolyPhiFace(rz, fEndPhi,   phiTotal / fNumSides, fStartPhi);
  }

  //
  // We might have dropped a face or two: recalculate numFace
  //
  numFace = face - faces;

  //
  // Make fEnclosingCylinder
  //

  /*
  double mxy = rz->Amax();
  double alfa = UUtils::kPi / fNumSides;

  double r= rz->Amax();

  if (fNumSides != 0)
  {
    // mxy *= std::sqrt(2.0); // this is old and wrong, works only for n = 4
    double k = std::tan(alfa) * mxy;
    double l = mxy / std::cos(alfa);
    mxy = l;
    r = l;
  }
  mxy += fgTolerance;
  */

  fEnclosingCylinder =
    new UEnclosingCylinder(rz->Amax(), rz->Bmax(), rz->Bmin(), fPhiIsOpen, phiStart, phiTotal);

  InitVoxels(*rz, fEnclosingCylinder->radius);

  fNoVoxels = fMaxSection < 2; // minimally, sections with at least numbers 0,1,2 values required, this corresponds to fMaxSection == 2
}



//
// Destructor
//
UPolyhedra::~UPolyhedra()
{
  delete [] fCorners;
//  if (fOriginalParameters) delete fOriginalParameters;

  delete fEnclosingCylinder;
}


//
// Copy constructor
//
UPolyhedra::UPolyhedra(const UPolyhedra& source)
  : UVCSGfaceted(source)
{
  CopyStuff(source);
}


//
// Assignment operator
//
UPolyhedra& UPolyhedra::operator=(const UPolyhedra& source)
{
  if (this == &source) return *this;

  UVCSGfaceted::operator=(source);

  delete [] fCorners;
//  if (fOriginalParameters) delete fOriginalParameters;

  delete fEnclosingCylinder;

  CopyStuff(source);

  return *this;
}


//
// CopyStuff
//
void UPolyhedra::CopyStuff(const UPolyhedra& source)
{
  //
  // Simple stuff
  //
  fNumSides    = source.fNumSides;
  fStartPhi   = source.fStartPhi;
  fEndPhi     = source.fEndPhi;
  fPhiIsOpen  = source.fPhiIsOpen;
  fNumCorner  = source.fNumCorner;
  fGenericPgon = source.fGenericPgon;

  //
  // The corner array
  //
  fCorners = new UPolyhedraSideRZ[fNumCorner];

  UPolyhedraSideRZ*  corn = fCorners,
                     *sourceCorn = source.fCorners;
  do
  {
    *corn = *sourceCorn;
  }
  while (++sourceCorn, ++corn < fCorners + fNumCorner);

  fOriginalParameters = source.fOriginalParameters;

  //
  // Enclosing cylinder
  //
  fEnclosingCylinder = new UEnclosingCylinder(*source.fEnclosingCylinder);
}


//
// Reset
//
// Recalculates and reshapes the solid, given pre-assigned scaled
// fOriginalParameters.
//
bool UPolyhedra::Reset()
{
  if (fGenericPgon)
  {
    std::ostringstream message;
    message << "Solid " << GetName() << " built using generic construct."
            << std::endl << "Not applicable to the generic construct !";
    UUtils::Exception("UPolyhedra::Reset(,,)", "GeomSolids1001",
                      UWarning, 1,  message.str().c_str());
    return 1;
  }

  //
  // Clear old setup
  //
  UVCSGfaceted::DeleteStuff();
  delete [] fCorners;
  delete fEnclosingCylinder;

  //
  // Rebuild polyhedra
  //
  UReduciblePolygon* rz =
    new UReduciblePolygon(&fOriginalParameters.Rmin[0],
                          &fOriginalParameters.Rmax[0],
                          &fOriginalParameters.fZValues[0],
                          fOriginalParameters.fNumZPlanes);
  Create(fOriginalParameters.fStartAngle,
         fOriginalParameters.fOpeningAngle,
         fOriginalParameters.fNumSide, rz);
  delete rz;

  return 0;
}


//
// Inside
//
// This is an override of UVCSGfaceted::Inside, created in order
// to speed things up by first checking with UEnclosingCylinder.
//
VUSolid::EnumInside UPolyhedra::Inside(const UVector3& p) const
{
  //
  // Quick test
  //
  if (fEnclosingCylinder->MustBeOutside(p)) return eOutside;

  //
  // Long answer
  //
  return UVCSGfaceted::Inside(p);
}


//
// DistanceToIn
//
double UPolyhedra::SafetyFromOutside(const UVector3& aPoint, bool aAccurate) const
{
  return UVCSGfaceted::SafetyFromOutside(aPoint, aAccurate);
}


//
// GetEntityType
//
UGeometryType UPolyhedra::GetEntityType() const
{
  return std::string("Polyhedra");
}


//
// Make a clone of the object
//
VUSolid* UPolyhedra::Clone() const
{
  return new UPolyhedra(*this);
}


//
// Stream object contents to an output stream
//
std::ostream& UPolyhedra::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "    *** Dump for solid - " << GetName() << " ***\n"
     << "    ===================================================\n"
     << " Solid type: UPolyhedra\n"
     << " Parameters: \n"
     << "    starting phi angle : " << fStartPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "    ending phi angle   : " << fEndPhi / (UUtils::kPi / 180.0) << " degrees \n";
  int i = 0;
  if (!fGenericPgon)
  {
    int numPlanes = fOriginalParameters.fNumZPlanes;
    os << "    number of Z planes: " << numPlanes << "\n"
       << "              Z values: \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << fOriginalParameters.fZValues[i] << "\n";
    }
    os << "              Tangent distances to inner surface (Rmin): \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << fOriginalParameters.Rmin[i] << "\n";
    }
    os << "              Tangent distances to outer surface (Rmax): \n";
    for (i = 0; i < numPlanes; i++)
    {
      os << "              Z plane " << i << ": "
         << fOriginalParameters.Rmax[i] << "\n";
    }
  }
  os << "    number of RZ points: " << fNumCorner << "\n"
     << "              RZ values (fCorners): \n";
  for (i = 0; i < fNumCorner; i++)
  {
    os << "                         "
       << fCorners[i].r << ", " << fCorners[i].z << "\n";
  }
  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}


//
// GetPointOnPlane
//
// Auxiliary method for get point on surface
//
UVector3 UPolyhedra::GetPointOnPlane(UVector3 p0, UVector3 p1,
                                     UVector3 p2, UVector3 p3) const
{
  double lambda1, lambda2, chose, aOne, aTwo;
  UVector3 t, u, v, w, Area, normal;
  aOne = 1.;
  aTwo = 1.;

  t = p1 - p0;
  u = p2 - p1;
  v = p3 - p2;
  w = p0 - p3;

  chose = UUtils::Random(0., aOne + aTwo);
  if ((chose >= 0.) && (chose < aOne))
  {
    lambda1 = UUtils::Random(0., 1.);
    lambda2 = UUtils::Random(0., lambda1);
    return (p2 + lambda1 * v + lambda2 * w);
  }

  lambda1 = UUtils::Random(0., 1.);
  lambda2 = UUtils::Random(0., lambda1);
  return (p0 + lambda1 * t + lambda2 * u);
}


//
// GetPointOnTriangle
//
// Auxiliary method for get point on surface
//
UVector3 UPolyhedra::GetPointOnTriangle(UVector3 p1,
                                        UVector3 p2,
                                        UVector3 p3) const
{
  double lambda1, lambda2;
  UVector3 v = p3 - p1, w = p1 - p2;

  lambda1 = UUtils::Random(0., 1.);
  lambda2 = UUtils::Random(0., lambda1);

  return (p2 + lambda1 * w + lambda2 * v);
}


//
// GetPointOnSurface
//
UVector3 UPolyhedra::GetPointOnSurface() const
{
  if (!fGenericPgon)   // Polyhedra by faces
  {
    int j, numPlanes = fOriginalParameters.fNumZPlanes, Flag = 0;
    double chose, totArea = 0., Achose1, Achose2,
                  rad1, rad2, sinphi1, sinphi2, cosphi1, cosphi2;
    double a, b, l2, rang, totalPhi, ksi,
           area, aTop = 0., aBottom = 0., zVal = 0.;

    UVector3 p0, p1, p2, p3;
    std::vector<double> aVector1;
    std::vector<double> aVector2;
    std::vector<double> aVector3;

    totalPhi = (fPhiIsOpen) ? (fEndPhi - fStartPhi) : 2 * UUtils::kPi;
    ksi = totalPhi / fNumSides;
    double cosksi = std::cos(ksi / 2.);

    // Below we generate the areas relevant to our solid
    //
    for (j = 0; j < numPlanes - 1; j++)
    {
      a = fOriginalParameters.Rmax[j + 1];
      b = fOriginalParameters.Rmax[j];
      l2 = UUtils::sqr(fOriginalParameters.fZValues[j]
                       - fOriginalParameters.fZValues[j + 1]) + UUtils::sqr(b - a);
      area = std::sqrt(l2 - UUtils::sqr((a - b) * cosksi)) * (a + b) * cosksi;
      aVector1.push_back(area);
    }

    for (j = 0; j < numPlanes - 1; j++)
    {
      a = fOriginalParameters.Rmin[j + 1]; //*cosksi;
      b = fOriginalParameters.Rmin[j];//*cosksi;
      l2 = UUtils::sqr(fOriginalParameters.fZValues[j]
                       - fOriginalParameters.fZValues[j + 1]) + UUtils::sqr(b - a);
      area = std::sqrt(l2 - UUtils::sqr((a - b) * cosksi)) * (a + b) * cosksi;
      aVector2.push_back(area);
    }

    for (j = 0; j < numPlanes - 1; j++)
    {
      if (fPhiIsOpen == true)
      {
        aVector3.push_back(0.5 * (fOriginalParameters.Rmax[j]
                                  - fOriginalParameters.Rmin[j]
                                  + fOriginalParameters.Rmax[j + 1]
                                  - fOriginalParameters.Rmin[j + 1])
                           *std::fabs(fOriginalParameters.fZValues[j + 1]
                                      - fOriginalParameters.fZValues[j]));
      }
      else
      {
        aVector3.push_back(0.);
      }
    }

    for (j = 0; j < numPlanes - 1; j++)
    {
      totArea += fNumSides * (aVector1[j] + aVector2[j]) + 2.*aVector3[j];
    }

    // Must include top and bottom areas
    //
    if (fOriginalParameters.Rmax[numPlanes - 1] != 0.)
    {
      a = fOriginalParameters.Rmax[numPlanes - 1];
      b = fOriginalParameters.Rmin[numPlanes - 1];
      l2 = UUtils::sqr(a - b);
      aTop = std::sqrt(l2 - UUtils::sqr((a - b) * cosksi)) * (a + b) * cosksi;
    }

    if (fOriginalParameters.Rmax[0] != 0.)
    {
      a = fOriginalParameters.Rmax[0];
      b = fOriginalParameters.Rmin[0];
      l2 = UUtils::sqr(a - b);
      aBottom = std::sqrt(l2 - UUtils::sqr((a - b) * cosksi)) * (a + b) * cosksi;
    }

    Achose1 = 0.;
    Achose2 = fNumSides * (aVector1[0] + aVector2[0]) + 2.*aVector3[0];

    chose = UUtils::Random(0., totArea + aTop + aBottom);
    if ((chose >= 0.) && (chose < aTop + aBottom))
    {
      chose = UUtils::Random(fStartPhi, fStartPhi + totalPhi);
      rang = std::floor((chose - fStartPhi) / ksi - 0.01);
      if (rang < 0)
      {
        rang = 0;
      }
      rang = std::fabs(rang);
      sinphi1 = std::sin(fStartPhi + rang * ksi);
      sinphi2 = std::sin(fStartPhi + (rang + 1) * ksi);
      cosphi1 = std::cos(fStartPhi + rang * ksi);
      cosphi2 = std::cos(fStartPhi + (rang + 1) * ksi);
      chose = UUtils::Random(0., aTop + aBottom);
      if (chose >= 0. && chose < aTop)
      {
        rad1 = fOriginalParameters.Rmin[numPlanes - 1];
        rad2 = fOriginalParameters.Rmax[numPlanes - 1];
        zVal = fOriginalParameters.fZValues[numPlanes - 1];
      }
      else
      {
        rad1 = fOriginalParameters.Rmin[0];
        rad2 = fOriginalParameters.Rmax[0];
        zVal = fOriginalParameters.fZValues[0];
      }
      p0 = UVector3(rad1 * cosphi1, rad1 * sinphi1, zVal);
      p1 = UVector3(rad2 * cosphi1, rad2 * sinphi1, zVal);
      p2 = UVector3(rad2 * cosphi2, rad2 * sinphi2, zVal);
      p3 = UVector3(rad1 * cosphi2, rad1 * sinphi2, zVal);
      return GetPointOnPlane(p0, p1, p2, p3);
    }
    else
    {
      for (j = 0; j < numPlanes - 1; j++)
      {
        if (((chose >= Achose1) && (chose < Achose2)) || (j == numPlanes - 2))
        {
          Flag = j;
          break;
        }
        Achose1 += fNumSides * (aVector1[j] + aVector2[j]) + 2.*aVector3[j];
        Achose2 = Achose1 + fNumSides * (aVector1[j + 1] + aVector2[j + 1])
                  + 2.*aVector3[j + 1];
      }
    }

    // At this point we have chosen a subsection
    // between to adjacent plane cuts...

    j = Flag;

    totArea = fNumSides * (aVector1[j] + aVector2[j]) + 2.*aVector3[j];
    chose = UUtils::Random(0., totArea);

    if ((chose >= 0.) && (chose < fNumSides * aVector1[j]))
    {
      chose = UUtils::Random(fStartPhi, fStartPhi + totalPhi);
      rang = std::floor((chose - fStartPhi) / ksi - 0.01);
      if (rang < 0)
      {
        rang = 0;
      }
      rang = std::fabs(rang);
      rad1 = fOriginalParameters.Rmax[j];
      rad2 = fOriginalParameters.Rmax[j + 1];
      sinphi1 = std::sin(fStartPhi + rang * ksi);
      sinphi2 = std::sin(fStartPhi + (rang + 1) * ksi);
      cosphi1 = std::cos(fStartPhi + rang * ksi);
      cosphi2 = std::cos(fStartPhi + (rang + 1) * ksi);
      zVal = fOriginalParameters.fZValues[j];

      p0 = UVector3(rad1 * cosphi1, rad1 * sinphi1, zVal);
      p1 = UVector3(rad1 * cosphi2, rad1 * sinphi2, zVal);

      zVal = fOriginalParameters.fZValues[j + 1];

      p2 = UVector3(rad2 * cosphi2, rad2 * sinphi2, zVal);
      p3 = UVector3(rad2 * cosphi1, rad2 * sinphi1, zVal);
      return GetPointOnPlane(p0, p1, p2, p3);
    }
    else if ((chose >= fNumSides * aVector1[j])
             && (chose <= fNumSides * (aVector1[j] + aVector2[j])))
    {
      chose = UUtils::Random(fStartPhi, fStartPhi + totalPhi);
      rang = std::floor((chose - fStartPhi) / ksi - 0.01);
      if (rang < 0)
      {
        rang = 0;
      }
      rang = std::fabs(rang);
      rad1 = fOriginalParameters.Rmin[j];
      rad2 = fOriginalParameters.Rmin[j + 1];
      sinphi1 = std::sin(fStartPhi + rang * ksi);
      sinphi2 = std::sin(fStartPhi + (rang + 1) * ksi);
      cosphi1 = std::cos(fStartPhi + rang * ksi);
      cosphi2 = std::cos(fStartPhi + (rang + 1) * ksi);
      zVal = fOriginalParameters.fZValues[j];

      p0 = UVector3(rad1 * cosphi1, rad1 * sinphi1, zVal);
      p1 = UVector3(rad1 * cosphi2, rad1 * sinphi2, zVal);

      zVal = fOriginalParameters.fZValues[j + 1];

      p2 = UVector3(rad2 * cosphi2, rad2 * sinphi2, zVal);
      p3 = UVector3(rad2 * cosphi1, rad2 * sinphi1, zVal);
      return GetPointOnPlane(p0, p1, p2, p3);
    }

    chose = UUtils::Random(0., 2.2);
    if ((chose >= 0.) && (chose < 1.))
    {
      rang = fStartPhi;
    }
    else
    {
      rang = fEndPhi;
    }

    cosphi1 = std::cos(rang);
    rad1 = fOriginalParameters.Rmin[j];
    sinphi1 = std::sin(rang);
    rad2 = fOriginalParameters.Rmax[j];

    p0 = UVector3(rad1 * cosphi1, rad1 * sinphi1,
                  fOriginalParameters.fZValues[j]);
    p1 = UVector3(rad2 * cosphi1, rad2 * sinphi1,
                  fOriginalParameters.fZValues[j]);

    rad1 = fOriginalParameters.Rmax[j + 1];
    rad2 = fOriginalParameters.Rmin[j + 1];

    p2 = UVector3(rad1 * cosphi1, rad1 * sinphi1,
                  fOriginalParameters.fZValues[j + 1]);
    p3 = UVector3(rad2 * cosphi1, rad2 * sinphi1,
                  fOriginalParameters.fZValues[j + 1]);
    return GetPointOnPlane(p0, p1, p2, p3);
  }
  else  // Generic polyhedra
  {
    return GetPointOnSurfaceGeneric();
  }
}

//
// UPolyhedraHistorical stuff
//
UPolyhedraHistorical::UPolyhedraHistorical()
  : fStartAngle(0.), fOpeningAngle(0.), fNumSide(0), fNumZPlanes(0),
    fZValues(0), Rmin(0), Rmax(0)
{
}

UPolyhedraHistorical::~UPolyhedraHistorical()
{
}

UPolyhedraHistorical::
UPolyhedraHistorical(const UPolyhedraHistorical& source)
{
  fStartAngle   = source.fStartAngle;
  fOpeningAngle = source.fOpeningAngle;
  fNumSide       = source.fNumSide;
  fNumZPlanes  = source.fNumZPlanes;

  fZValues = source.fZValues;
  Rmin = source.Rmin;
  Rmax = source.Rmax;
}

UPolyhedraHistorical&
UPolyhedraHistorical::operator=(const UPolyhedraHistorical& right)
{
  if (&right == this) return *this;

  fStartAngle   = right.fStartAngle;
  fOpeningAngle = right.fOpeningAngle;
  fNumSide       = right.fNumSide;
  fNumZPlanes  = right.fNumZPlanes;

  fZValues = right.fZValues;
  Rmin = right.Rmin;
  Rmax = right.Rmax;
  return *this;
}

void UPolyhedra::Extent(UVector3& aMin, UVector3& aMax) const
{
  fEnclosingCylinder->Extent(aMin, aMax);
}


//
// DistanceToIn
//
// This is an override of G4VCSGfaceted::Inside, created in order
// to speed things up by first checking with G4EnclosingCylinder.
//
double UPolyhedra::DistanceToIn(const UVector3& p,
                                const UVector3& v, double aPstep) const
{
  //
  // Quick test
  //
  if (fNoVoxels && fEnclosingCylinder->ShouldMiss(p, v))
    return UUtils::kInfinity;

  //
  // Long answer
  //
  return UVCSGfaceted::DistanceToIn(p, v, aPstep);
}
