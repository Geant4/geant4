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
// UQuadrangularFacet
//
// Class description:
//
// The UQuadrangularFacet class is used for the contruction of
// UTessellatedSolid.
// It is defined by four fVertices, which shall be in the same plane and be
// supplied in anti-clockwise order looking from the outsider of the solid
// where it belongs. Its constructor
//
//     UQuadrangularFacet (const UVector3 Pt0, const UVector3 vt1,
//                         const UVector3 vt2, const UVector3 vt3,
//                         UFacetVertexType);
//
// takes 5 parameters to define the four fVertices:
//   1) UFacetvertexType = "ABSOLUTE": in this case Pt0, vt1, vt2 and vt3
//      are the four fVertices required in anti-clockwise order when looking
//      from the outsider.
//   2) UFacetvertexType = "RELATIVE": in this case the first vertex is Pt0,
//      the second vertex is Pt0+vt, the third vertex is Pt0+vt2 and
//      the fourth vertex is Pt0+vt3, in anti-clockwise order when looking
//      from the outsider.
//
// 17.10.12 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UQuadrangularFacet_HH
#define UQuadrangularFacet_HH 1

#include "VUFacet.hh"
#include "UTriangularFacet.hh"
#include "UVector3.hh"

class UQuadrangularFacet : public VUFacet
{
  public:  // with description

    UQuadrangularFacet(const UVector3& Pt0, const UVector3& vt1,
                       const UVector3& vt2, const UVector3& vt3,
                       UFacetVertexType);
    virtual ~UQuadrangularFacet();

    UQuadrangularFacet(const UQuadrangularFacet& right);
    UQuadrangularFacet& operator=(const UQuadrangularFacet& right);

    VUFacet* GetClone();

    UVector3 Distance(const UVector3& p);
    double Distance(const UVector3& p, const double minDist);
    double Distance(const UVector3& p, const double minDist,
                    const bool outgoing);
    double Extent(const UVector3 axis);
    bool Intersect(const UVector3& p, const UVector3& v,
                   const bool outgoing, double& distance,
                   double& distFromSurface, UVector3& normal);

    double GetArea();
    UVector3 GetPointOnFace() const;

    virtual UGeometryType GetEntityType() const;

    inline int GetNumberOfVertices() const
    {
      return 4;
    }

    UVector3 GetVertex(int i) const
    {
      return i == 3 ? fFacet2.GetVertex(2) : fFacet1.GetVertex(i);
    }

    UVector3 GetSurfaceNormal() const;

    inline double GetRadius() const
    {
      return fRadius;
    }

    inline UVector3 GetCircumcentre() const
    {
      return fCircumcentre;
    }

    inline void SetVertex(int i, const UVector3& val)
    {
      switch (i)
      {
        case 0:
          fFacet1.SetVertex(0, val);
          fFacet2.SetVertex(0, val);
          break;
        case 1:
          fFacet1.SetVertex(1, val);
          break;
        case 2:
          fFacet1.SetVertex(2, val);
          fFacet2.SetVertex(1, val);
          break;
        case 3:
          fFacet2.SetVertex(2, val);
          break;
      }
    }

    inline void SetVertices(std::vector<UVector3>* v)
    {
      fFacet1.SetVertices(v);
      fFacet2.SetVertices(v);
    }

    inline bool IsDefined() const
    {
      return fFacet1.IsDefined();
    }

  protected:
  private:

    inline int GetVertexIndex(int i) const
    {
      return i == 3 ? fFacet2.GetVertexIndex(2) : fFacet1.GetVertexIndex(i);
    }

    inline void SetVertexIndex(int i, int val)
    {
      switch (i)
      {
        case 0:
          fFacet1.SetVertexIndex(0, val);
          fFacet2.SetVertexIndex(0, val);
          break;
        case 1:
          fFacet1.SetVertexIndex(1, val);
          break;
        case 2:
          fFacet1.SetVertexIndex(2, val);
          fFacet2.SetVertexIndex(1, val);
          break;
        case 3:
          fFacet2.SetVertexIndex(2, val);
          break;
      }
    }

    double fRadius;

    UVector3 fCircumcentre;

    int AllocatedMemory()
    {
      return sizeof(*this) + fFacet1.AllocatedMemory() + fFacet2.AllocatedMemory();
    }

    UTriangularFacet fFacet1, fFacet2;
};

#endif
