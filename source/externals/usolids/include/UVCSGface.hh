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
// UVCSGface
//
// Class description:
//
//   Definition of the virtual base class UVCSGface, one side (or face)
//   of a CSG-like solid. It should be possible to build a CSG entirely out of
//   connecting CSG faces.
//   Each face has an inside and outside surface, the former represents
//   the inside of the volume, the latter, the outside.
//
// 19.09.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UVCSGface_hh
#define UVCSGface_hh

#include "UTypes.hh"
#include "VUSolid.hh"

class UVoxelLimits;
class UAffineTransform;
class USolidExtentList;

class UVCSGface
{
  public: // with description

    UVCSGface() {}
    virtual ~UVCSGface() {}

    virtual bool Distance(const UVector3& p, const UVector3& v,
                          bool outgoing, double surfTolerance,
                          double& distance, double& distFromSurface,
                          UVector3& normal, bool& allBehind) = 0;

    virtual double Safety(const UVector3& p, bool outgoing) = 0;

    virtual VUSolid::EnumInside Inside(const UVector3& p, double tolerance,
                                       double* bestDistance) = 0;

    virtual UVector3 Normal(const UVector3& p,
                            double* bestDistance) = 0;

    virtual double Extent(const UVector3 axis) = 0;

    /*  virtual void CalculateExtent( const EAxisType axis,
                                    const UVoxelLimits &voxelLimit,
                                    const UAffineTransform &tranform,
                                    USolidExtentList &extentList       ) = 0;*/

    virtual UVCSGface* Clone() = 0;

    virtual double SurfaceArea() = 0;
    virtual UVector3 GetPointOnFace() = 0;
};

#endif
