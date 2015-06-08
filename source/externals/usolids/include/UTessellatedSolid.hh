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
// UTessellatedSolid
//
// Class description:
//
// UTessellatedSolid is a special Geant4 solid defined by a number of
// facets (UVFacet). It is important that the supplied facets shall form a
// fully enclose space which is the solid.
// Only two types of facet can be used for the construction of
// a UTessellatedSolid, i.e. the UTriangularFacet and UQuadrangularFacet.
//
// How to contruct a UTessellatedSolid:
//
// First declare a tessellated solid:
//
//      UTessellatedSolid* solidTarget = new UTessellatedSolid("Solid_name");
//
// Define the facets which form the solid
//
//      double targetSiz = 10*cm ;
//      UTriangularFacet *facet1 = new
//      UTriangularFacet (UVector3(-targetSize,-targetSize,        0.0),
//                         UVector3(+targetSize,-targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet2 = new
//      UTriangularFacet (UVector3(+targetSize,-targetSize,        0.0),
//                         UVector3(+targetSize,+targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet3 = new
//      UTriangularFacet (UVector3(+targetSize,+targetSize,        0.0),
//                         UVector3(-targetSize,+targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UTriangularFacet *facet4 = new
//      UTriangularFacet (UVector3(-targetSize,+targetSize,        0.0),
//                         UVector3(-targetSize,-targetSize,        0.0),
//                         UVector3(        0.0,        0.0,+targetSize),
//                         ABSOLUTE);
//      UQuadrangularFacet *facet5 = new
//      UQuadrangularFacet (UVector3(-targetSize,-targetSize,      0.0),
//                           UVector3(-targetSize,+targetSize,      0.0),
//                           UVector3(+targetSize,+targetSize,      0.0),
//                           UVector3(+targetSize,-targetSize,      0.0),
//                           ABSOLUTE);
//
// Then add the facets to the solid:
//
//      solidTarget->AddFacet((UVFacet*) facet1);
//      solidTarget->AddFacet((UVFacet*) facet2);
//      solidTarget->AddFacet((UVFacet*) facet3);
//      solidTarget->AddFacet((UVFacet*) facet4);
//      solidTarget->AddFacet((UVFacet*) facet5);
//
// Finally declare the solid is complete:
//
//      solidTarget->SetSolidClosed(true);
//
// 11.07.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UTessellatedSolid_hh
#define UTessellatedSolid_hh 1

#include <iostream>
#include <vector>
#include <set>
#include <map>

#include "VUSolid.hh"
#include "VUFacet.hh"
#include "UVoxelizer.hh"

struct UVertexInfo
{
  int id;
  double mag2;
};

class UVertexComparator
{
  public:
    bool operator()(const UVertexInfo& l, const UVertexInfo& r) const
    {
      return l.mag2 == r.mag2 ? l.id < r.id : l.mag2 < r.mag2;
    }
};

class UTessellatedSolid : public VUSolid
{
  public:

    UTessellatedSolid();
    virtual ~UTessellatedSolid();

    UTessellatedSolid(const std::string& name);

    UTessellatedSolid(__void__&);
    // Fake default constructor for usage restricted to direct object
    // persistency for clients requiring preallocation of memory for
    // persistifiable objects.

    UTessellatedSolid(const UTessellatedSolid& s);
    UTessellatedSolid& operator= (const UTessellatedSolid& s);
    UTessellatedSolid& operator+= (const UTessellatedSolid& right);

    bool AddFacet(VUFacet* aFacet);
    inline VUFacet* GetFacet(int i) const  { return fFacets[i]; }
    int GetNumberOfFacets() const;

    virtual double GetSurfaceArea();

    virtual VUSolid::EnumInside Inside(const UVector3& p) const;

    virtual bool Normal(const UVector3& p, UVector3& aNormal) const;

    virtual double SafetyFromOutside(const UVector3& p, bool aAccurate = false) const;

    virtual double SafetyFromInside(const UVector3& p, bool aAccurate = false) const;
    virtual UGeometryType GetEntityType() const;

    void SetSolidClosed(const bool t);

    bool GetSolidClosed() const;

    virtual UVector3 GetPointOnSurface() const;

    virtual std::ostream& StreamInfo(std::ostream& os) const;

    virtual double Capacity()  { return 0; }
    virtual double SurfaceArea()  { return GetSurfaceArea(); }

    inline virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}
    inline virtual void ComputeBBox(UBBox* /*aBox*/, bool /*aStore = false*/) {}

    inline void SetMaxVoxels(int max)  { fVoxels.SetMaxVoxels(max); }

    inline UVoxelizer& GetVoxels()  { return fVoxels; }

    virtual VUSolid* Clone() const;

    double      GetMinXExtent() const;
    double      GetMaxXExtent() const;
    double      GetMinYExtent() const;
    double      GetMaxYExtent() const;
    double      GetMinZExtent() const;
    double      GetMaxZExtent() const;

    virtual double DistanceToIn(const UVector3& p, const UVector3& v,
                                double aPstep = UUtils::kInfinity) const
    {
      return DistanceToInCore(p, v, aPstep);
    }

    virtual double DistanceToOut(const UVector3& p,
                                 const UVector3& v,
                                       UVector3& aNormalVector,
                                       bool&     aConvex,
                                       double aPstep = UUtils::kInfinity
                                ) const
    {
      return DistanceToOutCore(p, v, aNormalVector, aConvex, aPstep);
    }

    void Extent(UVector3& aMin, UVector3& aMax) const;

    int AllocatedMemoryWithoutVoxels();
    int AllocatedMemory();
    void DisplayAllocatedMemory();

  private:

    double DistanceToOutNoVoxels(const UVector3& p,
                                 const UVector3& v,
                                       UVector3& aNormalVector,
                                       bool&     aConvex,
                                 double aPstep = UUtils::kInfinity
                                ) const;

    double DistanceToInCandidates(const std::vector<int>& candidates, const UVector3& aPoint, const UVector3& aDirection /*, double aPstep, const UBits &bits*/) const;
    void DistanceToOutCandidates(const std::vector<int >& candidates, const UVector3& aPoint, const UVector3& direction, double& minDist, UVector3& minNormal, int& minCandidate/*, double aPstep*/ /*,  UBits &bits*/) const;
    double DistanceToInNoVoxels(const UVector3& p, const UVector3& v, double aPstep = UUtils::kInfinity) const;

    void SetExtremeFacets();

    VUSolid::EnumInside InsideNoVoxels(const UVector3& p) const;
    VUSolid::EnumInside InsideVoxels(const UVector3& aPoint) const;

    void Voxelize();

    void CreateVertexList();

    void PrecalculateInsides();

    void SetRandomVectors();

    double DistanceToInCore(const UVector3& p,
                            const UVector3& v,
                                  double aPstep = UUtils::kInfinity) const;
    double DistanceToOutCore(const UVector3& p,
                             const UVector3& v,
                                   UVector3& aNormalVector,
                                   bool&     aConvex,
                                   double aPstep = UUtils::kInfinity) const;

    int SetAllUsingStack(const std::vector<int>& voxel,
                         const std::vector<int>& max,
                               bool status, UBits& checked);

    void DeleteObjects();
    void CopyObjects(const UTessellatedSolid& s);

    static bool CompareSortedVoxel(const std::pair<int, double>& l,
                                   const std::pair<int, double>& r);

    double MinDistanceFacet(const UVector3& p, bool simple, VUFacet*& facet) const;

    inline bool OutsideOfExtent(const UVector3& p, double tolerance = 0) const
    {
      return (p.x() < fMinExtent.x() - tolerance || p.x() > fMaxExtent.x() + tolerance ||
              p.y() < fMinExtent.y() - tolerance || p.y() > fMaxExtent.y() + tolerance ||
              p.z() < fMinExtent.z() - tolerance || p.z() > fMaxExtent.z() + tolerance);
    }

    void Initialize();

  private:

    std::vector<VUFacet*>  fFacets;
    std::set<VUFacet*> fExtremeFacets;  // Does all other facets lie on or behind this surface?

    UGeometryType          fGeometryType;
    double                 fCubicVolume;
    double                 fSurfaceArea;

    std::vector<UVector3>  fVertexList;

    std::set<UVertexInfo, UVertexComparator> fFacetList;

    UVector3 fMinExtent, fMaxExtent;
    bool fSolidClosed;

    static const double dirTolerance;
    std::vector<UVector3> fRandir;

    double fgToleranceHalf;

    int fMaxTries;

    UVoxelizer fVoxels;  // voxelized solid

    UBits fInsides;
};

#endif
