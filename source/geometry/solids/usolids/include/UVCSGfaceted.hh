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
// Class description:
//
//   Virtual class defining CSG-like type shape that is built entire
//   of UCSGface faces.
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UVCSGfaceted_hh
#define UVCSGfaceted_hh

#include "VUSolid.hh"
#include "UVoxelizer.hh"
#include "UBox.hh"
#include "UReduciblePolygon.hh"

class UVCSGface;
class UVisExtent;

class UVCSGfaceted : public VUSolid
{
  public: // with description

    UVCSGfaceted(const std::string& name);
    virtual ~UVCSGfaceted();

    UVCSGfaceted(const UVCSGfaceted& source);
    UVCSGfaceted& operator=(const UVCSGfaceted& source);


    VUSolid::EnumInside InsideNoVoxels(const UVector3& p) const;

    virtual VUSolid::EnumInside Inside(const UVector3& p) const;

    virtual bool Normal(const UVector3& p, UVector3& n) const;

    double DistanceToInNoVoxels(const UVector3& p,
                                const UVector3& v) const;

    virtual double DistanceToIn(const UVector3& p,
                                const UVector3& v, double aPstep = UUtils::kInfinity) const;


    virtual double SafetyFromOutside(const UVector3& aPoint, bool aAccurate = false) const;


    double DistanceTo(const UVector3& p, const bool outgoing) const;

    double DistanceToOutNoVoxels(const UVector3& p,
                                 const UVector3&  v,
                                 UVector3&       n,
                                 bool&           aConvex) const;

    virtual double DistanceToOut(const UVector3& p,
                                 const UVector3&  v,
                                 UVector3&       n,
                                 bool&           aConvex,
                                 double aPstep = UUtils::kInfinity) const;


    virtual double SafetyFromInside(const UVector3& aPoint, bool aAccurate = false) const;

    virtual double SafetyFromInsideNoVoxels(const UVector3& aPoint, bool aAccurate = false) const;

    virtual UGeometryType GetEntityType() const;

    virtual std::ostream& StreamInfo(std::ostream& os) const;

    int GetCubVolStatistics() const;
    double GetCubVolEpsilon() const;
    void SetCubVolStatistics(int st);
    void SetCubVolEpsilon(double ep);
    int GetAreaStatistics() const;
    double GetAreaAccuracy() const;
    void SetAreaStatistics(int st);
    void SetAreaAccuracy(double ep);

    virtual double Capacity();
    // Returns an estimation of the geometrical cubic volume of the
    // solid. Caches the computed value once computed the first time.
    virtual double SurfaceArea();
    // Returns an estimation of the geometrical surface area of the
    // solid. Caches the computed value once computed the first time.

  public: // without description

  protected:  // without description

    double SafetyFromInsideSection(int index, const UVector3& p, UBits& bits) const;

    inline int GetSection(double z) const
    {
      int section = UVoxelizer::BinarySearch(fZs, z);
      if (section < 0) section = 0;
      else if (section > fMaxSection) section = fMaxSection;
      return section;
    }

    int   numFace;
    UVCSGface** faces;
    double fCubicVolume;
    double fSurfaceArea;


    std::vector<double> fZs; // z coordinates of given sections
    std::vector<std::vector<int> > fCandidates; // precalculated candidates for each of the section
    int fMaxSection; // maximum index number of sections of the solid (i.e. their number - 1). regular polyhedra with z = 1,2,3 section has 2 sections numbered 0 and 1, therefore the fMaxSection will be 1 (that is 2 - 1 = 1)
    mutable UBox fBox; // bounding box of the polyhedra, used in some methods
    double fBoxShift; // z-shift which is added during evaluation, because bounding box center does not have to be at (0,0,0)
    bool fNoVoxels; // if set to true, no voxelized algorithms will be used

    UVector3 GetPointOnSurfaceGeneric()const;
    // Returns a random point located on the surface of the solid
    // in case of generic Polycone or generic Polyhedra.

    void CopyStuff(const UVCSGfaceted& source);
    void DeleteStuff();

    void FindCandidates(double z, std::vector <int>& candidates, bool sides = false);

    void InitVoxels(UReduciblePolygon& z, double radius);

  private:

    int   fStatistics;
    double fCubVolEpsilon;
    double fAreaAccuracy;
    // Statistics, error accuracy for volume estimation.

};

#endif
