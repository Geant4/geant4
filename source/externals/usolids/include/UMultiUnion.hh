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
// UMultiUnion
//
// Class description:
//
//   An instance of "UMultiUnion" constitutes a grouping of several solids
//   deriving from the "VUSolid" mother class. The subsolids are stored with
//   their respective location in an instance of "UNode". An instance of
//   "UMultiUnion" is subsequently composed of one or several nodes.
//
// 19.10.12 Marek Gayer
// --------------------------------------------------------------------

#ifndef USOLIDS_UMultiUnion
#define USOLIDS_UMultiUnion

#include <vector>

#include "VUSolid.hh"
#include "UUtils.hh"
#include "UTransform3D.hh"
#include "UBits.hh"
#include "UVoxelizer.hh"

class UMultiUnion : public VUSolid
{
    friend class UVoxelizer;

  public:
    UMultiUnion() : VUSolid() {}
    UMultiUnion(const std::string& name);
    ~UMultiUnion();

    // Build the multiple union by adding nodes
    void AddNode(VUSolid& solid, UTransform3D& trans);

    UMultiUnion(const UMultiUnion& rhs);
    UMultiUnion& operator=(const UMultiUnion& rhs);

    // Accessors
    inline const UTransform3D& GetTransformation(int index) const;
    inline VUSolid* GetSolid(int index) const;
    inline int GetNumberOfSolids()const;

    // Navigation methods
    EnumInside Inside(const UVector3& aPoint) const;

    EnumInside InsideIterator(const UVector3& aPoint) const;

    double SafetyFromInside(const UVector3& aPoint,
                            bool aAccurate = false) const;

    double SafetyFromOutside(const UVector3& aPoint,
                             bool aAccurate = false) const;

    double DistanceToInNoVoxels(const UVector3& aPoint,
                                const UVector3& aDirection,
                                double aPstep = UUtils::kInfinity) const;

    double DistanceToIn(const UVector3& aPoint,
                        const UVector3& aDirection,
                        double aPstep) const;

    double DistanceToOut(const UVector3& aPoint,
                         const UVector3& aDirection,
                         UVector3&       aNormalVector,
                         bool&           aConvex,
                         double aPstep = UUtils::kInfinity) const;

    double DistanceToOutVoxels(const UVector3& aPoint,
                               const UVector3& aDirection,
                               UVector3&       aNormalVector,
                               bool&           aConvex,
                               double aPstep = UUtils::kInfinity) const;

    double DistanceToOutVoxelsCore(const UVector3& aPoint,
                                   const UVector3& aDirection,
                                   UVector3&       aNormalVector,
                                   bool&           aConvex,
                                   std::vector<int>& candidates) const;

    double DistanceToOutNoVoxels(const UVector3& aPoint,
                                 const UVector3& aDirection,
                                 UVector3&       aNormalVector,
                                 bool&           aConvex,
                                 double aPstep = UUtils::kInfinity) const;

    bool Normal(const UVector3& aPoint, UVector3& aNormal) const;

    void Extent(EAxisType aAxis, double& aMin, double& aMax) const;
    void Extent(UVector3& aMin, UVector3& aMax) const;

    double Capacity();
    double SurfaceArea();

    VUSolid* Clone() const ;

    UGeometryType GetEntityType() const { return "MultiUnion"; }
    void ComputeBBox(UBBox* aBox, bool aStore = false);

    virtual void GetParametersList(int /*aNumber*/, double* /*aArray*/) const {}

    // Finalize and prepare for use. User MUST call it once before
    // navigation use.
    void Voxelize();
    EnumInside InsideNoVoxels(const UVector3& aPoint) const;

    inline UVoxelizer& GetVoxels() const;


    std::ostream& StreamInfo(std::ostream& os) const;

    UVector3 GetPointOnSurface() const;

  private:

    void SetVoxelFinder(const UVoxelizer& finder);
    EnumInside InsideWithExclusion(const UVector3& aPoint, UBits* bits = NULL) const;
    int SafetyFromOutsideNumberNode(const UVector3& aPoint, bool aAccurate, double& safety) const;
    double DistanceToInCandidates(const UVector3& aPoint, const UVector3& aDirection, double aPstep, std::vector<int>& candidates, UBits& bits) const;

    std::vector<VUSolid*> fSolids;
    std::vector<UTransform3D> fTransformObjs;
    UVoxelizer fVoxels;  // Pointer to the vozelized solid
    double       fCubicVolume;   // Cubic Volume
    double       fSurfaceArea;   // Surface Area
};

inline UVoxelizer&  UMultiUnion:: GetVoxels() const
{
  return (UVoxelizer&)fVoxels;
}
inline const UTransform3D& UMultiUnion::GetTransformation(int index) const
{
  return fTransformObjs[index];
}
inline VUSolid* UMultiUnion::GetSolid(int index) const
{
  return fSolids[index];
}
inline int UMultiUnion::GetNumberOfSolids()const
{
  return fSolids.size();
}

#endif
