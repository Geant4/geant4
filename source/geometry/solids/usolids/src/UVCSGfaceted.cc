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
// UVCSGfaceted
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#include "UUtils.hh"
#include <string>
#include <cmath>
#include <sstream>
#include "UVCSGfaceted.hh"
#include "UVCSGface.hh"
#include "UVoxelizer.hh"
#include "UReduciblePolygon.hh"

using namespace std;

//
// Constructor
//
UVCSGfaceted::UVCSGfaceted(const std::string& name)
  : VUSolid(name),
    numFace(0), faces(0), fCubicVolume(0.), fSurfaceArea(0.),
    fMaxSection(0),fBoxShift(0.), fNoVoxels(true),fStatistics(1000000), fCubVolEpsilon(0.001), fAreaAccuracy(-1.)
{
}

//
// Destructor
//
UVCSGfaceted::~UVCSGfaceted()
{
  DeleteStuff();
}

//
// Copy constructor
//
UVCSGfaceted::UVCSGfaceted(const UVCSGfaceted& source)
  : VUSolid(source)
{
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  CopyStuff(source);
}


//
// Assignment operator
//
UVCSGfaceted& UVCSGfaceted::operator=(const UVCSGfaceted& source)
{
  if (&source == this)
  {
    return *this;
  }

  // Copy base class data
  //
  VUSolid::operator=(source);

  // Copy data
  //
  fStatistics = source.fStatistics;
  fCubVolEpsilon = source.fCubVolEpsilon;
  fAreaAccuracy = source.fAreaAccuracy;

  CopyStuff(source);

  return *this;
}


//
// CopyStuff (protected)
//
// Copy the contents of source
//
void UVCSGfaceted::CopyStuff(const UVCSGfaceted& source)
{
  numFace = source.numFace;
  if (numFace == 0)
  {
    return;  // odd, but permissable?
  }

  faces = new UVCSGface*[numFace];

  UVCSGface** face = faces,
              **sourceFace = source.faces;
  do
  {
    *face = (*sourceFace)->Clone();
  }
  while (++sourceFace, ++face < faces + numFace);
  fCubicVolume = source.fCubicVolume;
  fSurfaceArea = source.fSurfaceArea;

  fMaxSection = source.fMaxSection;
  fNoVoxels = source.fNoVoxels;
  fZs = source.fZs;
  fBox = source.fBox;
  fBoxShift = source.fBoxShift;
}


//
// DeleteStuff (protected)
//
// Delete all allocated objects
//
void UVCSGfaceted::DeleteStuff()
{
  if (numFace)
  {
    UVCSGface** face = faces;
    do
    {
      delete *face;
    }
    while (++face < faces + numFace);

    delete [] faces;
  }
}

//
// Inside
//
// It could be a good idea to override this virtual
// member to add first a simple test (such as spherical
// test or whatnot) and to call this version only if
// the simplier test fails.
//

/*
VUSolid::EnumInside UVCSGfaceted::Inside( const UVector3 &p ) const
{
  VUSolid::EnumInside answer=eOutside;
  UVCSGface **face = faces;
  double best = UUtils::kInfinity;
  do
  {
    double distance;
    VUSolid::EnumInside result = (*face)->Inside( p, fgTolerance*0.5, &distance );
    if (result == eSurface) { return eSurface; }
    if (distance < best)
    {
      best = distance;
      answer = result;
    }
  } while( ++face < faces + numFace );

  return answer;
}
*/


//
// Inside
//
// It could be a good idea to override this virtual
// member to add first a simple test (such as spherical
// test or whatnot) and to call this version only if
// the simplier test fails.
//
VUSolid::EnumInside UVCSGfaceted::InsideNoVoxels(const UVector3& p) const
{
  VUSolid::EnumInside answer = eOutside;
  UVCSGface** face = faces;
  double best = UUtils::kInfinity;
  do
  {
    double distance;
    VUSolid::EnumInside result = (*face)->Inside(p, fgTolerance * 0.5, &distance);
    if (result == eSurface)
    {
      return eSurface;
    }
    if (distance < best)
    {
      best = distance;
      answer = result;
    }
  }
  while (++face < faces + numFace);

  return answer;
}



//
// SurfaceNormal
//

bool UVCSGfaceted::Normal(const UVector3& p, UVector3& n) const
{
  UVector3 answer;
  double best = UUtils::kInfinity;
  UVector3 normal;

  UBits bits(numFace);

  int index = GetSection(p.z);
  const vector<int>& candidates = fCandidates[index];
  int size = candidates.size();
  for (int i = 0; i < size; ++i)
  {
    int candidate = candidates[i];
    if (!bits[candidate])
    {
      bits.SetBitNumber(candidate);
      UVCSGface& face = *faces[candidate];
      double distance;
      normal = face.Normal(p, &distance);
      if (distance < best)
      {
        best = distance;
        answer = normal;
        if (distance < fgTolerance)
          break;
      }
    }
  }
  n = answer;
  return true;
}

/*

// non voxelized version:

  bool UVCSGfaceted::Normal( const UVector3 &p, UVector3 &n) const
  {
    UVector3 answer;
    double best = UUtils::kInfinity;
    UVector3 normal;

    for (int i = 0; i < numFace; ++i)
    {
      UVCSGface &face = *faces[i];
      double distance;
      normal = face.Normal( p, &distance);
      if (distance < best)
      {
        best = distance;
        answer = normal;
      }
    }
    n = answer;
    return true;
  }
*/

//
// DistanceToIn(p,v)
//

double UVCSGfaceted::DistanceToInNoVoxels(const UVector3& p,
                                          const UVector3& v) const
{
  double distance = UUtils::kInfinity;
  double distFromSurface = UUtils::kInfinity;
  UVCSGface** face = faces;
  UVCSGface* bestFace = *face;
  static double htol = fgTolerance * 0.5;
  UVector3 faceNormal;

  do
  {
    double   faceDistance, faceDistFromSurface;
    bool    faceAllBehind;
    if ((*face)->Distance(p, v, false, htol,
                          faceDistance, faceDistFromSurface,
                          faceNormal, faceAllBehind))
    {
      //
      // Intersecting face
      //
      if (faceDistance < distance)
      {
        distance = faceDistance;
        distFromSurface = faceDistFromSurface;
        bestFace = *face;
        if (distFromSurface <= 0)
        {
          return 0;
        }
      }
    }
  }
  while (++face < faces + numFace);

  if (distance < UUtils::kInfinity && distFromSurface < htol)
  {
    if (bestFace->Safety(p, false) < htol)
    {
      distance = 0;
    }
  }

  return distance;
}




//
// DistanceToOut(p,v)
//

double UVCSGfaceted::DistanceToOutNoVoxels(const UVector3& p, const UVector3&  v, UVector3& n, bool& aConvex) const
{
  bool allBehind = true;
  double distance = UUtils::kInfinity;
  double distFromSurface = UUtils::kInfinity;
  UVector3 normal, faceNormal;

  UVCSGface** face = faces;
  UVCSGface* bestFace = *face;
  do
  {
    double faceDistance, faceDistFromSurface;
    bool faceAllBehind;
    if ((*face)->Distance(p, v, true, fgTolerance / 2,
                          faceDistance, faceDistFromSurface,
                          faceNormal, faceAllBehind))
    {
      // Intersecting face
      if ((distance < UUtils::kInfinity) || (!faceAllBehind))
      {
        allBehind = false;
      }
      if (faceDistance < distance)
      {
        distance = faceDistance;
        distFromSurface = faceDistFromSurface;
        normal = faceNormal;
        bestFace = *face;
        if (distFromSurface <= 0)
        {
          break;
        }
      }
    }
  }
  while (++face < faces + numFace);

  if (distance < UUtils::kInfinity)
  {
    if (distFromSurface <= 0)
    {
      distance = 0;
    }
    else if (distFromSurface < fgTolerance / 2)
    {
      if (bestFace->Safety(p, true) < fgTolerance / 2)
      {
        distance = 0;
      }
    }

    aConvex = allBehind;
    n = normal;
  }
  else
  {
    if (Inside(p) == eSurface)
    {
      distance = 0;
    }
    aConvex = false;
  }

  return distance;
}




//
// DistanceTo
//
// Protected routine called by DistanceToIn and DistanceToOut
//
double UVCSGfaceted::DistanceTo(const UVector3& p,
                                const bool outgoing) const
{
  UVCSGface** face = faces;
  double best = UUtils::kInfinity;
  do
  {
    double distance = (*face)->Safety(p, outgoing);
    if (distance < best)  best = distance;
  }
  while (++face < faces + numFace);

  return (best < 0.5 * fgTolerance) ? 0 : best;
}

//
// GetEntityType
//
UGeometryType UVCSGfaceted::GetEntityType() const
{
  return std::string("UCSGfaceted");
}


//
// Stream object contents to an output stream
//
std::ostream& UVCSGfaceted::StreamInfo(std::ostream& os) const
{
  os << "-----------------------------------------------------------\n"
     << "		*** Dump for solid - " << GetName() << " ***\n"
     << "		===================================================\n"
     << " Solid type: UVCSGfaceted\n"
     << " Parameters: \n"
     << "		number of faces: " << numFace << "\n"
     << "-----------------------------------------------------------\n";

  return os;
}


//
// GetCubVolStatistics
//
int UVCSGfaceted::GetCubVolStatistics() const
{
  return fStatistics;
}


//
// GetCubVolEpsilon
//
double UVCSGfaceted::GetCubVolEpsilon() const
{
  return fCubVolEpsilon;
}


//
// SetCubVolStatistics
//
void UVCSGfaceted::SetCubVolStatistics(int st)
{
  fCubicVolume = 0.;
  fStatistics = st;
}


//
// SetCubVolEpsilon
//
void UVCSGfaceted::SetCubVolEpsilon(double ep)
{
  fCubicVolume = 0.;
  fCubVolEpsilon = ep;
}


//
// GetAreaStatistics
//
int UVCSGfaceted::GetAreaStatistics() const
{
  return fStatistics;
}


//
// GetAreaAccuracy
//
double UVCSGfaceted::GetAreaAccuracy() const
{
  return fAreaAccuracy;
}


//
// SetAreaStatistics
//
void UVCSGfaceted::SetAreaStatistics(int st)
{
  fSurfaceArea = 0.;
  fStatistics = st;
}


//
// SetAreaAccuracy
//
void UVCSGfaceted::SetAreaAccuracy(double ep)
{
  fSurfaceArea = 0.;
  fAreaAccuracy = ep;
}


//
// Capacity
//
double UVCSGfaceted::Capacity()
{
  if (fCubicVolume != 0.)
  {
    ;
  }
  else
  {
    fCubicVolume = EstimateCubicVolume(fStatistics, fCubVolEpsilon);
  }
  return fCubicVolume;
}


//
// SurfaceArea
//
double UVCSGfaceted::SurfaceArea()
{
  if (fSurfaceArea != 0.)
  {
    ;
  }
  else
  {
    fSurfaceArea = EstimateSurfaceArea(fStatistics, fAreaAccuracy);
  }
  return fSurfaceArea;
}

//
// GetPointOnSurfaceGeneric proportional to Areas of faces
// in case of GenericPolycone or GenericPolyhedra
//
UVector3 UVCSGfaceted::GetPointOnSurfaceGeneric() const
{
  // Preparing variables
  //
  UVector3 answer = UVector3(0., 0., 0.);
  UVCSGface** face = faces;
  double area = 0;
  int i;
  std::vector<double> areas;

  // First step: calculate surface areas
  //
  do
  {
    double result = (*face)->SurfaceArea();
    areas.push_back(result);
    area = area + result;
  }
  while (++face < faces + numFace);

  // Second Step: choose randomly one surface
  //
  UVCSGface** face1 = faces;
  double chose = area * UUtils::Random();
  double Achose1, Achose2;
  Achose1 = 0;
  Achose2 = 0.;
  i = 0;

  do
  {
    Achose2 += areas[i];
    if (chose >= Achose1 && chose < Achose2)
    {
      UVector3 point;
      point = (*face1)->GetPointOnFace();
      return point;
    }
    i++;
    Achose1 = Achose2;
  }
  while (++face1 < faces + numFace);

  return answer;
}

double UVCSGfaceted::SafetyFromOutside(const UVector3& p, bool accurate) const
{
  if (!accurate)
  {
    UVector3 pb(p.x, p.y, p.z - fBoxShift);
    return fBox.SafetyFromOutside(pb);
  }
  return DistanceTo(p, false);
}

double UVCSGfaceted::SafetyFromInsideNoVoxels(const UVector3& p, bool) const
{
  return DistanceTo(p, true);
}


void  UVCSGfaceted::InitVoxels(UReduciblePolygon& rz, double radius)
{
  int size = rz.NumVertices() + 1;
  vector<double> r(size), z(size), zs;
  rz.CopyVertices(&r[0], &z[0]);

  fZs.clear();
  for (int i = 0; i < size; ++i)
  {
    double v = z[i];
    if (std::find(fZs.begin(), fZs.end(), v) == fZs.end())
    {
      fZs.push_back(v);
    }
    std::sort(fZs.begin(), fZs.end());
  }

  size = fZs.size();
  fMaxSection = size - 2;

  for (int i = 0; i <= fMaxSection; ++i)
  {
    vector<int> candidates;

    double left = fZs[i], right = fZs[i + 1];
    double middle = (left + right) / 2;
    FindCandidates(middle, candidates);

    FindCandidates(left, candidates, true);
    FindCandidates(right, candidates, true);

    fCandidates.push_back(candidates);
  }

  fBox.Set(radius, radius, (fZs.back() - fZs.front()) / 2);
  fBoxShift = fZs[0] + fBox.GetZHalfLength();
}

void UVCSGfaceted::FindCandidates(double z, vector <int>& candidates, bool sides)
{
  for (int j = 0; j < numFace; j++)
  {
    UVCSGface* face = faces[j];
    double minZ = -face->Extent(UVector3(0, 0, -1)) ;
    double maxZ = face->Extent(UVector3(0, 0, 1));
    if (z >= minZ - fgTolerance * 10 && z <= maxZ  + fgTolerance * 10)
    {
      if (!sides || std::fabs(minZ - maxZ) < fgTolerance * 10)
        if (std::find(candidates.begin(), candidates.end(), j) == candidates.end())
          candidates.push_back(j);
    }
  }
}


double UVCSGfaceted::DistanceToIn(const UVector3& p, const UVector3& v, double /*aPstep*/) const
{
  if (fNoVoxels) return DistanceToInNoVoxels(p, v);

  UVector3 pb(p.x, p.y, p.z - fBoxShift);

  double idistance, shift;
  idistance = fBox.DistanceToIn(pb, v); // using only box, this appears
  // to be faster than: idistance = enclosingCylinder->DistanceTo(pb, v);
  if (idistance >= UUtils::kInfinity) return idistance;

  // this line can be here or not. not a big difference in performance
  // TODO: fix enclosingCylinder for polyhedra!!! - the current radius appears to be too small
  //  if (enclosingCylinder->ShouldMiss(p, v)) return UUtils::kInfinity;

  // this just takes too much time
  //  idistance = enclosingCylinder->DistanceTo(pb, v);
  //  if (idistance == UUtils::kInfinity) return idistance;

  double z = p.z + idistance * v.z;
  int index = GetSection(z);
  int increment = (v.z > 0) ? 1 : -1;
  if (std::fabs(v.z) < fgTolerance) increment = 0;

  double distance = UUtils::kInfinity;
  double distFromSurface = UUtils::kInfinity;
  UVCSGface* bestFace = 0;
  UBits bits(numFace);
  UVector3 faceNormal;

  do
  {
    const vector<int>& candidates = fCandidates[index];
    int size = candidates.size();
    for (int i = 0; i < size; ++i)
    {
      int candidate = candidates[i];
      if (!bits[candidate])
      {
        bits.SetBitNumber(candidate);
        UVCSGface& face = *faces[candidate];

        double   faceDistance,
                 faceDistFromSurface;
        bool    faceAllBehind;
        if (face.Distance(p, v, false, fgTolerance * 0.5,
                          faceDistance, faceDistFromSurface,
                          faceNormal, faceAllBehind))
        {
          // Intersecting face
          if (faceDistance < distance)
          {
            distance = faceDistance;
            distFromSurface = faceDistFromSurface;
            bestFace = &face;
            if (distFromSurface <= 0) return 0;
          }
        }
      }
    }

    if (!increment)
      break;

    index += increment;
    if (index < 0 || index > fMaxSection)
      break;
    int newz = increment > 0 ? index  : index + 1;
    shift = (fZs[newz] - z) / v.z;
  }
  while (idistance + shift < distance);

  if (distance < UUtils::kInfinity && distFromSurface < fgTolerance / 2)
  {
    if (bestFace->Safety(p, false) < fgTolerance / 2)
    {
      distance = 0;
    }
  }

  return distance;
}



double UVCSGfaceted::DistanceToOut(const UVector3& p, const UVector3&  v, UVector3& n, bool& aConvex, double /*aPstep*/) const
{
  if (fNoVoxels) return DistanceToOutNoVoxels(p, v, n, aConvex);

  int index =  GetSection(p.z);
  int increment = (v.z > 0) ? 1 : -1;

  bool allBehind = true;
  double distance = UUtils::kInfinity;
  double distFromSurface = UUtils::kInfinity;
  UVector3 normal, faceNormal;
  double shift;

  UVCSGface* bestFace = 0;
  UBits bits(numFace);

  do
  {
    const vector<int>& candidates = fCandidates[index];
    int size = candidates.size();
    for (int i = 0; i < size; ++i)
    {
      int candidate = candidates[i];
      if (!bits[candidate])
      {
        bits.SetBitNumber(candidate);
        UVCSGface& face = *faces[candidate];

        double faceDistance, faceDistFromSurface;
        bool faceAllBehind;
        if ((face.Distance(p, v, true, fgTolerance * 0.5, faceDistance, faceDistFromSurface,
                           faceNormal, faceAllBehind)))
        {
          // Intersecting face
          if ((distance < UUtils::kInfinity) || (!faceAllBehind))
          {
            allBehind = false;
          }
          if (faceDistance < distance)
          {
            distance = faceDistance;
            distFromSurface = faceDistFromSurface;
            normal = faceNormal;
            bestFace = &face;
            if (distFromSurface <= 0) break;
          }
        }
      }
    }

    if (distFromSurface <= 0) break;
    if (!increment) break;

    index += increment;
    if (index < 0 || index > fMaxSection)
      break;
    int newz = increment > 0 ? index  : index + 1;
    shift = (fZs[newz] - p.z) / v.z;
  }
  while (shift < distance);

  if (distance < UUtils::kInfinity)
  {
    if (distFromSurface <= 0)
    {
      distance = 0;
    }
    else if (distFromSurface < fgTolerance / 2)
    {
      if (bestFace->Safety(p, true) < fgTolerance * 0.5)
      {
        distance = 0;
      }
    }

    aConvex = allBehind;
    n = normal;
  }
  else
  {
    if (Inside(p) == eSurface)
    {
      distance = 0;
    }
    aConvex = false;
  }

  return distance;
}


VUSolid::EnumInside UVCSGfaceted::Inside(const UVector3& p) const
{
  if (fNoVoxels) return InsideNoVoxels(p);

//  if (fEnclosingCylinder->MustBeOutside(p)) return eOutside;

  int index =  GetSection(p.z);
  double shift;

  UBits bits(numFace);
  double best = UUtils::kInfinity;
  VUSolid::EnumInside answer = eOutside;
  int middle = index;

  do
  {
    const vector<int>& candidates = fCandidates[index];
    int size = candidates.size();
    for (int i = 0; i < size; ++i)
    {
      int candidate = candidates[i];
      if (!bits[candidate])
      {
        UVCSGface& face = *faces[candidate];

        double distance;
        VUSolid::EnumInside result = face.Inside(p, fgTolerance * 0.5, &distance);
        if (result == eSurface) return eSurface;
        if (distance < best)
        {
          best = distance;
          answer = result;
        }
        bits.SetBitNumber(candidate);
      }
    }

    if (index <= middle)
    {
      if (--index >= 0)
      {
        shift = fZs[index + 1] - p.z;
        if (shift < best) continue;
      }
      index = middle;
    }
    if (++index > fMaxSection) break;
    shift = p.z - fZs[index];
  }
  while (shift > best);

  return answer;
}

double UVCSGfaceted::SafetyFromInsideSection(int index, const UVector3& p, UBits& bits) const
{
  double best = UUtils::kInfinity;

  const vector<int>& candidates = fCandidates[index];
  int size = candidates.size();
  for (int i = 0; i < size; ++i)
  {
    int candidate = candidates[i];
    if (!bits[candidate])
    {
      bits.SetBitNumber(candidate);
      UVCSGface& face = *faces[candidate];

      double distance = face.Safety(p, true);
      if (distance < best)  best = distance;
    }
  }
  return best;
}

double UVCSGfaceted::SafetyFromInside(const UVector3& p, bool) const
{
  if (fNoVoxels) return SafetyFromInsideNoVoxels(p);

  int index = UVoxelizer::BinarySearch(fZs, p.z);
  if (index < 0 || index > fMaxSection) return 0;

  UBits bits(numFace);
  double minSafety = SafetyFromInsideSection(index, p, bits);

  if (minSafety > UUtils::kInfinity) return 0;
  if (minSafety < 1e-6) return 0;

  double zbase = fZs[index + 1];
  for (int i = index + 1; i <= fMaxSection; ++i)
  {
    double dz = fZs[i] - zbase;
    if (dz >= minSafety) break;
    double safety = SafetyFromInsideSection(i, p, bits);
    if (safety < minSafety) minSafety = safety;
  }

  if (index > 0)
  {
    zbase = fZs[index - 1];
    for (int i = index - 1; i >= 0; --i)
    {
      double dz = zbase - fZs[i];
      if (dz >= minSafety) break;
      double safety = SafetyFromInsideSection(i, p, bits);
      if (safety < minSafety) minSafety = safety;
    }
  }
  return (minSafety < 0.5 * fgTolerance) ? 0 : minSafety;
}
