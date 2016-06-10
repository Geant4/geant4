
#include "VUSolid.hh"

////////////////////////////////////////////////////////////////////////////////
//  "Universal" Solid Interface
//  Authors: J. Apostolakis, G. Cosmo, M. Gayer, A. Gheata, A. Munnich, T. Nikitina (CERN)
//
//  Created: 25 May 2011
//
////////////////////////////////////////////////////////////////////////////////

double VUSolid::fgTolerance = 1.0E-9;  // cartesian tolerance; to be changed (for U was 1e-8, but we keep Geant4)
double VUSolid::frTolerance = 1.0E-9;  // radial tolerance; to be changed

double VUSolid::faTolerance = 1.0E-9;  // angular tolerance; to be changed

//______________________________________________________________________________
VUSolid::VUSolid() : fName()
{

}

//______________________________________________________________________________
VUSolid::VUSolid(const std::string &name) :
fName(name)
{
	// Named constructor
	SetName(name);
}
 
//______________________________________________________________________________
VUSolid::~VUSolid()
{

}

/*
int UIntersectingCone::LineHitsCone2( const UVector3 &p,
	const UVector3 &v,
	double &s1, double &s2 )
{
	double x0 = p.x, y0 = p.y, z0 = p.z;
	double tx = v.x, ty = v.y, tz = v.z;

	// Special case which might not be so rare: B = 0 (precisely)
	//
	if (B==0)
	{
		if (std::fabs(tz) < 1/UUtils::kInfinity)	{ return 0; }

		s1 = (A-z0)/tz;
		return 1;
	}

	double B2 = B*B;

	double a = tz*tz - B2*(tx*tx + ty*ty);
	double b = 2*( (z0-A)*tz - B2*(x0*tx + y0*ty) );
	double c = UUtils::sqr(z0-A) - B2*( x0*x0 + y0*y0 );

	double radical = b*b - 4*a*c;

	if (radical < -1E-6*std::fabs(b)) { return 0; }	 // No solution

	if (radical < 1E-6*std::fabs(b))
	{
		//
		// The radical is roughly zero: check for special, very rare, cases
		//
		if (std::fabs(a) > 1/UUtils::kInfinity)
		{
			if ( std::fabs(x0*ty - y0*tx) < std::fabs(1E-6/B) )
			{
				s1 = -0.5*b/a;
				return 1;
			}
			return 0;
		}
	}
	else
	{
		radical = std::sqrt(radical);
	}

	if (a < -1/UUtils::kInfinity)
	{
		double sa, sb, q = -0.5*( b + (b < 0 ? -radical : +radical) );
		sa = q/a;
		sb = c/q;
		if (sa < sb) { s1 = sa; s2 = sb; } else { s1 = sb; s2 = sa; }
		if ((z0 + (s1)*tz - A)/B < 0)	{ return 0; }
		return 2;
	}
	else if (a > 1/UUtils::kInfinity)
	{
		double sa, sb, q = -0.5*( b + (b < 0 ? -radical : +radical) );
		sa = q/a;
		sb = c/q;
		s1 = (tz*B > 0)^(sa > sb) ? sb : sa;
		return 1;
	}
	else if (std::fabs(b) < 1/UUtils::kInfinity)
	{
		return 0;
	}
	else
	{
		s1 = -c/b;
		if ((z0 + (s1)*tz - A)/B < 0)	{ return 0; }
		return 1;
	}
}
int UIntersectingCone::LineHitsCone2( const UVector3 &p,
	const UVector3 &v,
	double &s1, double &s2 )
{
	double x0 = p.x, y0 = p.y, z0 = p.z;
	double tx = v.x, ty = v.y, tz = v.z;

	// Special case which might not be so rare: B = 0 (precisely)
	//
	if (B==0)
	{
		if (std::fabs(tz) < 1/UUtils::kInfinity)	{ return 0; }

		s1 = (A-z0)/tz;
		return 1;
	}

	double B2 = B*B;

	double a = tz*tz - B2*(tx*tx + ty*ty);
	double b = 2*( (z0-A)*tz - B2*(x0*tx + y0*ty) );
	double c = UUtils::sqr(z0-A) - B2*( x0*x0 + y0*y0 );

	double radical = b*b - 4*a*c;

	if (radical < -1E-6*std::fabs(b)) { return 0; }	 // No solution

	if (radical < 1E-6*std::fabs(b))
	{
		//
		// The radical is roughly zero: check for special, very rare, cases
		//
		if (std::fabs(a) > 1/UUtils::kInfinity)
		{
			if ( std::fabs(x0*ty - y0*tx) < std::fabs(1E-6/B) )
			{
				s1 = -0.5*b/a;
				return 1;
			}
			return 0;
		}
	}
	else
	{
		radical = std::sqrt(radical);
	}

	if (a < -1/UUtils::kInfinity)
	{
		double sa, sb, q = -0.5*( b + (b < 0 ? -radical : +radical) );
		sa = q/a;
		sb = c/q;
		if (sa < sb) { s1 = sa; s2 = sb; } else { s1 = sb; s2 = sa; }
		if ((z0 + (s1)*tz - A)/B < 0)	{ return 0; }
		return 2;
	}
	else if (a > 1/UUtils::kInfinity)
	{
		double sa, sb, q = -0.5*( b + (b < 0 ? -radical : +radical) );
		sa = q/a;
		sb = c/q;
		s1 = (tz*B > 0)^(sa > sb) ? sb : sa;
		return 1;
	}
	else if (std::fabs(b) < 1/UUtils::kInfinity)
	{
		return 0;
	}
	else
	{
		s1 = -c/b;
		if ((z0 + (s1)*tz - A)/B < 0)	{ return 0; }
		return 1;
	}
}
*/

////////////////////////////////////////////////////////////////
//
// Returns an estimation of the solid volume in internal units.
// The number of statistics and error accuracy is fixed.
// This method may be overloaded by derived classes to compute the
// exact geometrical quantity for solids where this is possible.
// or anyway to cache the computed value.
// This implementation does NOT cache the computed value.

double VUSolid::Capacity()
{
  int cubVolStatistics = 1000000;
  double cubVolEpsilon = 0.001;
  return EstimateCubicVolume(cubVolStatistics, cubVolEpsilon);
}

////////////////////////////////////////////////////////////////
//
// Calculate cubic volume based on Inside() method.
// Accuracy is limited by the second argument or the statistics
// expressed by the first argument.
// Implementation is courtesy of Vasiliki Despoina Mitsou,
// University of Athens.

double VUSolid::EstimateCubicVolume(int nStat, double epsilon) const
{
  int iInside=0;
  double px,py,pz,volume;
  UVector3 min,max;
  UVector3 p;
  VUSolid::EnumInside in;

  // values needed for CalculateExtent signature

  // min max extents of pSolid along X,Y,Z

  this->Extent(min,max);

  // limits

  if(nStat < 100)		nStat	 = 100;
  if(epsilon > 0.01) epsilon = 0.01;

  for(int i = 0; i < nStat; i++ )
  {
    px = min.x+(max.x-min.x)*UUtils::Random();
    py = min.y+(max.y-min.y)*UUtils::Random();
    pz = min.z+(max.z-min.z)*UUtils::Random();
    p	= UVector3(px,py,pz);
    in = this->Inside(p);
    if(in != eOutside) iInside++;		
  }
  volume = (max.x-min.x)*(max.y-min.y)*(max.z-min.z)*iInside/nStat;
  return volume;
}

////////////////////////////////////////////////////////////////
//
// Returns an estimation of the solid surface area in internal units.
// The number of statistics and error accuracy is fixed.
// This method may be overloaded by derived classes to compute the
// exact geometrical quantity for solids where this is possible.
// or anyway to cache the computed value.
// This implementation does NOT cache the computed value.

double VUSolid::SurfaceArea()
{
  int stat = 1000000;
  double ell = -1.;
  return EstimateSurfaceArea(stat,ell);
}

////////////////////////////////////////////////////////////////
//
// Estimate surface area based on Inside(), DistanceToIn(), and
// DistanceToOut() methods. Accuracy is limited by the statistics
// defined by the first argument. Implemented by Mikhail Kosov.

double VUSolid::EstimateSurfaceArea(int nStat, double ell) const
{
  int inside=0;
  double px,py,pz,surf;
  UVector3 min,max;
  UVector3 p;
  VUSolid::EnumInside in;

  // values needed for CalculateExtent signature

  // min max extents of pSolid along X,Y,Z

  this->Extent(min,max);

  // limits

  if(nStat < 100) { nStat = 100; }

  double dX=max.x-min.x;
  double dY=max.y-min.y;
  double dZ=max.z-min.z;
  if(ell<=0.)					// Automatic definition of skin thickness
  {
    double minval=dX;
    if(dY<dX) { minval=dY; }
    if(dZ<minval) { minval=dZ; }
    ell=.01*minval;
  }

  double dd=2*ell;
  min.x-=ell; min.y-=ell; min.z-=ell; dX+=dd; dY+=dd; dZ+=dd;

  for(int i = 0; i < nStat; i++ )
  {
    px = min.x+dX*UUtils::Random();
    py = min.y+dY*UUtils::Random();
    pz = min.z+dZ*UUtils::Random();
    p	= UVector3(px,py,pz);
    in = this->Inside(p);
    if(in != eOutside)
    {
      if	(SafetyFromInside(p)<ell) { inside++; }
    }
    else if(SafetyFromOutside(p)<ell) { inside++; }
  }
  // @@ The conformal correction can be upgraded
  surf = dX*dY*dZ*inside/dd/nStat;
  return surf;
} 

void VUSolid::ExtentAxis(EAxisType aAxis, double &aMin, double &aMax) const
// Returns the minimum and maximum extent along the specified Cartesian axis
{
  // Returns extent of the solid along a given cartesian axis
  if (aAxis >= 0 && aAxis <= 2)
  {
    UVector3 min, max;
    Extent(min,max);
    aMin = min[aAxis]; aMax = max[aAxis];
  }
#ifdef USPECSDEBUG
  else
    cout << "Extent: unknown axis" << aAxis << std::endl;
#endif
}

void VUSolid::SetCarTolerance(double eps)
{
  fgTolerance=eps;
}
void VUSolid::SetRadTolerance(double eps)
{
  frTolerance=eps;
}
void VUSolid::SetAngTolerance(double eps)
{
 faTolerance=eps;
}


