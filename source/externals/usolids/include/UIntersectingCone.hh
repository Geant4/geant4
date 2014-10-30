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
// UIntersectingCone
//
// Class description:
//
//   Utility class which calculates the intersection
//   of an arbitrary line with a fixed cone
//
// 19.02.13 Marek Gayer
//          Created from original implementation in Geant4
// --------------------------------------------------------------------

#ifndef UIntersectingCone_hh
#define UIntersectingCone_hh

#include "UTypes.hh"

class UIntersectingCone
{
  public:

    UIntersectingCone(const double r[2], const double z[2]);
    virtual ~UIntersectingCone();

    int LineHitsCone(const UVector3& p, const UVector3& v, double& s1, double& s2);

    bool HitOn(const double r, const double z);

    inline double RLo() const
    {
      return rLo;
    }
    inline double RHi() const
    {
      return rHi;
    }
    inline double ZLo() const
    {
      return zLo;
    }
    inline double ZHi() const
    {
      return zHi;
    }

  public: // without description

    /*
    UIntersectingCone(__void__&);
      // Fake default constructor for usage restricted to direct object
      // persistency for clients requiring preallocation of memory for
      // persistifiable objects.
    */


  protected:

    double zLo, zHi,  // Z bounds of side
           rLo, rHi;  // R bounds of side

    bool   type1;   // True if cone is type 1
    // (std::fabs(z1-z2)>std::fabs(r1-r2))
    double A, B;     // Cone radius parameter:
    // type 1: r = A + B*z
    // type 2: z = A + B*r

    // int Solution (const UVector3 &p, const UVector3 &v, double a, double b, double c, double &s1, double &s2);

    int LineHitsCone1(const UVector3& p, const UVector3& v,
                      double& s1, double& s2);

    int LineHitsCone1Optimized(const UVector3& p, const UVector3& v,
                               double& s1, double& s2);

    int LineHitsCone2(const UVector3& p, const UVector3& v,
                      double& s1, double& s2);

//    const double kInfinity;
    const static double EpsilonQuad;
};

#endif
