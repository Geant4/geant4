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
// UGenericPolycone
//
// 19.10.13 Tatiana Nikitina
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>

#include "UGenericPolycone.hh"
#include "UPolyconeSide.hh"
#include "UPolyPhiFace.hh"
#include "UEnclosingCylinder.hh"
#include "UReduciblePolygon.hh"

using namespace std;
//
// Constructor (generic parameters)
//
UGenericPolycone::UGenericPolycone(const std::string& name,
                                   double phiStart,
                                   double phiTotal,
                                   int   numRZ,
                                   const double r[],
                                   const double z[])
  : UVCSGfaceted(name)
{
  UReduciblePolygon* rz = new UReduciblePolygon(r, z, numRZ);

  Create(phiStart, phiTotal, rz);

  delete rz;
}


//
// Create
//
// Generic create routine, called by each constructor after
// conversion of arguments
//
void UGenericPolycone::Create(double phiStart,
                              double phiTotal,
                              UReduciblePolygon* rz)
{
  //
  // Perform checks of rz values
  //
  if (rz->Amin() < 0.0)
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "				All R values must be >= 0 !";
    UUtils::Exception("UGenericPolycone::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  double rzArea = rz->Area();
  if (rzArea < -VUSolid::Tolerance())
    rz->ReverseOrder();

  else if (rzArea < -VUSolid::Tolerance())
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "				R/Z Cross section is zero or near zero: " << rzArea;
    UUtils::Exception("UGenericPolycone::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  if ((!rz->RemoveDuplicateVertices(VUSolid::Tolerance()))
      || (!rz->RemoveRedundantVertices(VUSolid::Tolerance())))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "				Too few unique R/Z values !";
    UUtils::Exception("UGenericPolycone::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  if (rz->CrossesItself(1 / UUtils::kInfinity))
  {
    std::ostringstream message;
    message << "Illegal input parameters - " << GetName() << std::endl
            << "				R/Z segments Cross !";
    UUtils::Exception("UGenericPolycone::Create()", "GeomSolids0002",
                      UFatalErrorInArguments, 1, message.str().c_str());
  }

  numCorner = rz->NumVertices();

  //
  // Phi opening? Account for some possible roundoff, and interpret
  // nonsense value as representing no phi opening
  //
  if (phiTotal <= 0 || phiTotal > 2 * UUtils::kPi - 1E-10)
  {
    phiIsOpen = false;
    startPhi = 0;
    endPhi = 2 * UUtils::kPi;
  }
  else
  {
    phiIsOpen = true;

    //
    // Convert phi into our convention
    //
    startPhi = phiStart;
    while (startPhi < 0) startPhi += 2 * UUtils::kPi;

    endPhi = phiStart + phiTotal;
    while (endPhi < startPhi) endPhi += 2 * UUtils::kPi;
  }

  //
  // Allocate corner array.
  //
  corners = new UPolyconeSideRZ[numCorner];

  //
  // Copy corners
  //
  UReduciblePolygonIterator iterRZ(rz);

  UPolyconeSideRZ* next = corners;
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
  numFace = phiIsOpen ? numCorner + 2 : numCorner;
  faces = new UVCSGface*[numFace];

  //
  // Construct conical faces
  //
  // But! Don't construct a face if both points are at zero radius!
  //
  UPolyconeSideRZ* corner = corners,
                   *prev = corners + numCorner - 1,
                    *nextNext;
  UVCSGface** face = faces;
  do
  {
    next = corner + 1;
    if (next >= corners + numCorner) next = corners;
    nextNext = next + 1;
    if (nextNext >= corners + numCorner) nextNext = corners;

    if (corner->r < 1 / UUtils::kInfinity && next->r < 1 / UUtils::kInfinity) continue;

    //
    // We must decide here if we can dare declare one of our faces
    // as having a "valid" normal (i.e. allBehind = true). This
    // is never possible if the face faces "inward" in r.
    //
    bool allBehind;
    if (corner->z > next->z)
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
      allBehind = !rz->BisectedBy(corner->r, corner->z,
                                  next->r, next->z, VUSolid::Tolerance());
    }

    *face++ = new UPolyconeSide(prev, corner, next, nextNext,
                                startPhi, endPhi - startPhi, phiIsOpen, allBehind);
  }
  while (prev = corner, corner = next, corner > corners);

  if (phiIsOpen)
  {
    //
    // Construct phi open edges
    //
    *face++ = new UPolyPhiFace(rz, startPhi, 0, endPhi);
    *face++ = new UPolyPhiFace(rz, endPhi,   0, startPhi);
  }

  //
  // We might have dropped a face or two: recalculate numFace
  //
  numFace = face - faces;

  //
  // Make enclosingCylinder
  //
  enclosingCylinder =
    new UEnclosingCylinder(rz->Amax(), rz->Bmax(), rz->Bmin(), phiIsOpen, phiStart, phiTotal);

  InitVoxels(*rz, enclosingCylinder->radius);

  fNoVoxels = fMaxSection < 2;

}


//
// Fake default constructor - sets only member data and allocates memory
//                            for usage restricted to object persistency.
//

/*
UGenericPolycone::UGenericPolycone( __void__& a )
  : UVCSGfaceted(a), startPhi(0.),  endPhi(0.), phiIsOpen(false),
    genericPcon(false), numCorner(0), corners(0),
    fOriginalParameters(0), enclosingCylinder(0)
{
}
*/


//
// Destructor
//
UGenericPolycone::~UGenericPolycone()
{
  delete [] corners;
  delete enclosingCylinder;
}


//
// Copy constructor
//
UGenericPolycone::UGenericPolycone(const UGenericPolycone& source)
  : UVCSGfaceted(source)
{
  CopyStuff(source);
}


//
// Assignment operator
//
UGenericPolycone& UGenericPolycone::operator=(const UGenericPolycone& source)
{
  if (this == &source) return *this;

  UVCSGfaceted::operator=(source);

  delete [] corners;

  delete enclosingCylinder;

  CopyStuff(source);

  return *this;
}


//
// CopyStuff
//
void UGenericPolycone::CopyStuff(const UGenericPolycone& source)
{
  //
  // Simple stuff
  //
  startPhi  = source.startPhi;
  endPhi    = source.endPhi;
  phiIsOpen = source.phiIsOpen;
  numCorner = source.numCorner;

  //
  // The corner array
  //
  corners = new UPolyconeSideRZ[numCorner];

  UPolyconeSideRZ* corn = corners,
                   *sourceCorn = source.corners;
  do
  {
    *corn = *sourceCorn;
  }
  while (++sourceCorn, ++corn < corners + numCorner);

  //
  // Enclosing cylinder
  //
  enclosingCylinder = new UEnclosingCylinder(*source.enclosingCylinder);
}


//
// Reset
//
bool UGenericPolycone::Reset()
{

  std::ostringstream message;
  message << "Solid " << GetName() << " built using generic construct."
          << std::endl << "Not applicable to the generic construct !";
  // UException("UGenericPolycone::Reset(,,)", "GeomSolids1001",
  //            JustWarning, message, "Parameters NOT resetted.");
  return 1;

}


//
// Inside
//
// This is an override of UVCSGfaceted::Inside, created in order
// to speed things up by first checking with UEnclosingCylinder.
//
VUSolid::EnumInside UGenericPolycone::Inside(const UVector3& p) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->MustBeOutside(p)) return eOutside;

  //
  // Long answer
  //
  return UVCSGfaceted::Inside(p);
}


//
// DistanceToIn
//
// This is an override of UVCSGfaceted::Inside, created in order
// to speed things up by first checking with UEnclosingCylinder.
//
double UGenericPolycone::DistanceToIn(const UVector3& p,
                                      const UVector3& v, double aPstep) const
{
  //
  // Quick test
  //
  if (enclosingCylinder->ShouldMiss(p, v))
    return UUtils::kInfinity;

  //
  // Long answer
  //
  return UVCSGfaceted::DistanceToIn(p, v, aPstep);
}

//
// GetEntityType
//
UGeometryType UGenericPolycone::GetEntityType() const
{
  return std::string("GenericPolycone");
}

//
// Make a clone of the object
//
VUSolid* UGenericPolycone::Clone() const
{
  return new UGenericPolycone(*this);
}

//
// Stream object contents to an output stream
//
std::ostream& UGenericPolycone::StreamInfo(std::ostream& os) const
{
  int oldprc = os.precision(16);
  os << "-----------------------------------------------------------\n"
     << "		*** Dump for solid - " << GetName() << " ***\n"
     << "		===================================================\n"
     << " Solid type: UGenericPolycone\n"
     << " Parameters: \n"
     << "		starting phi angle : " << startPhi / (UUtils::kPi / 180.0) << " degrees \n"
     << "		ending phi angle	 : " << endPhi / (UUtils::kPi / 180.0) << " degrees \n";
  int i = 0;

  os << "    number of RZ points: " << numCorner << "\n"
     << "              RZ values (corners): \n";
  for (i = 0; i < numCorner; i++)
  {
    os << "                         "
       << corners[i].r << ", " << corners[i].z << "\n";
  }

  os << "-----------------------------------------------------------\n";
  os.precision(oldprc);

  return os;
}

//
// GetPointOnSurface
//
UVector3 UGenericPolycone::GetPointOnSurface() const
{
  return GetPointOnSurfaceGeneric();

}

void UGenericPolycone::Extent(UVector3& aMin, UVector3& aMax) const
{
  enclosingCylinder->Extent(aMin, aMax);
}

