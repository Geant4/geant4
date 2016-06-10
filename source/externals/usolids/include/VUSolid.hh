#ifndef USOLIDS_VUSolid
#define USOLIDS_VUSolid
////////////////////////////////////////////////////////////////////////////////
//  "Universal" Solid Interface
//  Authors: J. Apostolakis, G. Cosmo, M. Gayer, A. Gheata, A. Munnich, T. Nikitina (CERN)
//
//  Created: 25 May 2011
//
////////////////////////////////////////////////////////////////////////////////

#include "UTypes.hh"
#include "UVector3.hh"

#include "UUtils.hh"

#define USOLIDS
#define USOLIDSONLY

class VUSolid
{
public:

enum EnumInside { eInside=0, eSurface=1, eOutside=2 }; 
   // Use eInside < eSurface < eOutside: allows "max(,)" to combine Inside of surfaces
   // Potentially replace eSurface with eInSurface, eOutSurface 

enum EAxisType { eXaxis=0, eYaxis=1, eZaxis=2};

protected:
static double     fgTolerance;
static double     frTolerance;
static double     faTolerance;

// =>10 degrees/wedge for complete tube

public:
  VUSolid();
  VUSolid(const std::string &name);
  virtual ~VUSolid();

  // Accessors and modifiers for Tolerance
  inline double GetCarTolerance() const;
  inline double GetRadTolerance() const;
  inline double GetAngTolerance() const;
  void SetCarTolerance(double eps);
  void SetRadTolerance(double eps);
  void SetAngTolerance(double eps);
  
  // Navigation methods
  virtual EnumInside Inside (const UVector3 &aPoint) const = 0;
  //
  // Evaluate if point is inside, outside or on the surface within the tolerance
  virtual double  SafetyFromInside ( const UVector3 &aPoint, 
                                     bool aAccurate=false) const = 0;
  virtual double  SafetyFromOutside( const UVector3 &aPoint, 
                                     bool aAccurate=false) const = 0;
  // 
  // Estimates isotropic distance to the surface of the solid. This must
  // be either accurate or an underestimate. 
  //  Two modes: - default/fast mode, sacrificing accuracy for speed
  //             - "precise" mode,  requests accurate value if available.
  //  For both modes, if at a large distance from solid ( > ? ) 
  //    it is expected that a simplified calculation will be made if available.

  virtual double DistanceToIn(  const UVector3 &aPoint, 
				                    const UVector3 &aDirection,
                                double aPstep = UUtils::kInfinity) const = 0;
  virtual double DistanceToOut( const UVector3 &aPoint,
                                const UVector3  &aDirection,
                                UVector3       &aNormalVector,
                                bool           &aConvex,
                                double aPstep = UUtils::kInfinity) const = 0;
  // 
  // o return the exact distance (double) from a surface, given a direction
  // o compute the normal on the surface, returned as argument, calculated
  //      within the method to verify if it is close to the surface or not
  // o for DistanceToOut(), normal-vector and convexity flag could be optional (to decide).
  //    If normal cannot be computed (or shape is not convex), set 'convex' to 'false'.
  // o for DistanceToIn(), the normal-vector could be added as optional

  virtual bool Normal( const UVector3& aPoint, UVector3 &aNormal ) const = 0; 
  // Computes the normal on a surface and returns it as a unit vector
  //   In case a point is further than tolerance_normal from a surface, set validNormal=false
  //   Must return a valid vector. (even if the point is not on the surface.)
  //
  //   On an edge or corner, provide an average normal of all facets within tolerance

  // Decision: provide or not the Boolean 'validNormal' argument for returning validity 

  virtual void ExtentAxis(EAxisType aAxis, double &aMin, double &aMax) const;
  
  virtual void Extent( UVector3 &aMin, UVector3 &aMax ) const = 0;
  // Return the minimum and maximum extent along all Cartesian axes
  // For both the Extent methods
  // o Expect mostly to use a GetBBox()/CalculateBBox() method internally to compute the extent
  // o Decision: whether to store the computed BBox (containing or representing 6 double values), 
  //     and whether to compute it at construction time.
  // Methods are *not* const to allow caching of the Bounding Box
 virtual UGeometryType  GetEntityType() const = 0;
       // Provide identification of the class of an object.
       // (required for persistency and STEP interface)
  
  const std::string &GetName() const {return fName;}
  void            SetName(const std::string &aName) {fName = aName;}

  // Auxiliary methods
  virtual double Capacity()    = 0 ;   //  like  CubicVolume() 
  virtual double SurfaceArea() = 0 ;
  //  Expect the solids to cache the values of Capacity and Surface Area 
  
  // Sampling
  virtual void    SamplePointsInside(int /*aNpoints*/, UVector3 * /*aArray*/) const {}
  virtual void    SamplePointsOnSurface(int /*aNpoints*/, UVector3 * /*aArray*/) const {}
  virtual void    SamplePointsOnEdge(int /*aNpoints*/, UVector3 * /*aArray*/) const {}
  // o generates points on the edges of a solid - primarily for testing purposes
  // o for solids composed only of curved surfaces(like full spheres or toruses) or 
  //        where an implementation is not available, it defaults to PointOnSurface.

  // Visualisation
  virtual void GetParametersList(int aNumber,double *aArray) const =0;
 
  virtual VUSolid* Clone() const =0; 
  // o provide a new object which is a clone of the solid
  
  // Visualization
   
  static double   Tolerance() {return fgTolerance;}

  virtual std::ostream& StreamInfo( std::ostream& os ) const = 0;

  virtual UVector3 GetPointOnSurface() const = 0;

  double EstimateCubicVolume(int nStat, double epsilon) const;
  // Calculate cubic volume based on Inside() method.
  // Accuracy is limited by the second argument or the statistics
  // expressed by the first argument.

  double EstimateSurfaceArea(int nStat, double ell) const;
  // Calculate surface area only based on Inside() method.
  // Accuracy is limited by the second argument or the statistics
  // expressed by the first argument.

protected:
  virtual void    ComputeBBox(UBBox *aBox, bool aStore = false) = 0;
  // o Compute the bounding box for the solid. Called automatically and stored ? 
   // o Can throw an exception if the solid is invalid
private:
  std::string fName;  // Name of the solid
  //UBBox          *fBBox;   // Bounding box
};
inline double VUSolid::GetCarTolerance() const { return fgTolerance;}
inline double VUSolid::GetRadTolerance() const { return frTolerance;}
inline double VUSolid::GetAngTolerance() const { return faTolerance;}

#endif

