// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Hype.cc,v 1.1 1999-01-07 16:07:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4Hype: this class implements in G4 the volume equivalent 
//               to the HYPE volume in Geant 3, i.e. a tube with 
//               hyperbolic profile.
//
//  Authors: 
//      Ernesto Lamanna (Ernesto.Lamanna@roma1.infn.it) &
//      Francesco Safai Tehrani (Francesco.SafaiTehrani@roma1.infn.it)
//      Rome, INFN & University of Rome "La Sapienza",  9 June 1998.
//
// $ Original: G4Hype.cc,v 1.0 1998/06/09 16:57:50 safai Exp $
//
// class G4Hype: this class implements in G4 the volume equivalent 
//               to the HYPE volume in Geant 3, i.e. a tube with 
//               hyperbolic profile.
//
// For further informations, please read G4Hype.history and G4Hype.doc,
// and the G4Hype.hh header.

#include "G4Hype.hh"

#include "G4VoxelLimits.hh"
#include "G4AffineTransform.hh"

#include "G4VPVParameterisation.hh"

#include "meshdefs.hh"

#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBS.hh"
#include "G4NURBStube.hh"
#include "G4NURBScylinder.hh"
#include "G4NURBStubesector.hh"
#include "G4VisExtent.hh"

#define  SURFACE_PRECISION  (0.5*kCarTolerance)

// Constructor - check parameters, and fills protected data members
G4Hype::G4Hype(const G4String& pName,
	   const G4double newInnerRadius,
	   const G4double newOuterRadius,
	   const G4double newInnerStereo,
	   const G4double newOuterStereo,
	   const G4double newHalfLenZ) : G4CSGSolid(pName)
{
// Check z-len
    if (newHalfLenZ>0)
	{ 
	  halfLenZ=newHalfLenZ;   
	}
    else
	{
	    G4Exception("Error in G4Hype::G4Hype - invalid z half-length");
	}

// Check radii
    if (newInnerRadius>=0 && newOuterRadius>=0)
      if (newInnerRadius < newOuterRadius) { 
	innerRadius=newInnerRadius;
	outerRadius=newOuterRadius;
      }
      else { // swapping radii  (:-)
	innerRadius=newOuterRadius;
	outerRadius=newInnerRadius;
      }
    else
	{
	    G4Exception("Error in G4Hype::G4Hype - invalid radii");
	}
    
    innerStereo=newInnerStereo;
    outerStereo=newOuterStereo;

    // init of precalculated quantities
    tanInnerStereo2=tan(innerStereo)*tan(innerStereo); 
    tanOuterStereo2=tan(outerStereo)*tan(outerStereo);
    innerRadius2=innerRadius*innerRadius;
    outerRadius2=outerRadius*outerRadius;
    endInnerRadius2=HypeInnerRadius2(halfLenZ);
    endOuterRadius2=HypeOuterRadius2(halfLenZ);
    endInnerRadius=sqrt(endInnerRadius2);
    endOuterRadius=sqrt(endOuterRadius2);
     
}

// Destructor
G4Hype::~G4Hype()
{;}

// Dispatch to parameterisation for replication mechanism dimension
// computation & modification.
 void G4Hype::ComputeDimensions(G4VPVParameterisation* p,
                              const G4int n,
                              const G4VPhysicalVolume* pRep)
{
    p->ComputeDimensions(*this,n,pRep);
}

// Calculate extent under transform and specified limit
G4bool G4Hype::CalculateExtent(const EAxis pAxis,
			      const G4VoxelLimits& pVoxelLimit,
			      const G4AffineTransform& pTransform,
			      G4double& pMin, G4double& pMax) const
{


  //G4cout << "Calculate extent !"<< endl;
  //G4cout << "Surface data: " << outerRadius << endl;

  if (!pTransform.IsRotated())
	{
	  // Special case handling for unrotated solid tubes
	  // Compute x/y/z mins and maxs fro bounding box respecting limits,
	  // with early returns if outside limits. Then switch() on pAxis,
	  // and compute exact x and y limit for x/y case
	    
	  G4double xoffset,xMin,xMax;
	  G4double yoffset,yMin,yMax;
	  G4double zoffset,zMin,zMax;

	  G4double diff1,diff2,maxDiff,newMin,newMax;
	  G4double xoff1,xoff2,yoff1,yoff2;

	  xoffset=pTransform.NetTranslation().x();
	  xMin=xoffset-endOuterRadius;
	  xMax=xoffset+endOuterRadius;
	  //G4cout << "xMin, xMax : " << xMin << " " << xMax << endl;
	  if (pVoxelLimit.IsXLimited())
	    {
	      if (xMin>pVoxelLimit.GetMaxXExtent() ||
		  xMax<pVoxelLimit.GetMinXExtent())
		{
		  return false;
		}
	      else
		{
		  if (xMin<pVoxelLimit.GetMinXExtent())
		    {
		      xMin=pVoxelLimit.GetMinXExtent();
		      //G4cout << xMin << endl;
		    }
		  if (xMax>pVoxelLimit.GetMaxXExtent())
		    {
		      xMax=pVoxelLimit.GetMaxXExtent();
		      //G4cout << xMax << endl;
		    }
		}
	    }
	  
	  yoffset=pTransform.NetTranslation().y();
	  yMin=yoffset-endOuterRadius;
	  yMax=yoffset+endOuterRadius;
	  if (pVoxelLimit.IsYLimited())
	    {
	      if (yMin>pVoxelLimit.GetMaxYExtent()
		  ||yMax<pVoxelLimit.GetMinYExtent())
		{
		  return false;
		}
	      else
		{
		  if (yMin<pVoxelLimit.GetMinYExtent())
		    {
		      yMin=pVoxelLimit.GetMinYExtent();
		    }
		  if (yMax>pVoxelLimit.GetMaxYExtent())
		    {
		      yMax=pVoxelLimit.GetMaxYExtent();
		    }
		}
	    }
	  
	  
	  zoffset=pTransform.NetTranslation().z();
	  zMin=zoffset-halfLenZ;
	  zMax=zoffset+halfLenZ;
	  if (pVoxelLimit.IsZLimited())
	    {
	      if (zMin>pVoxelLimit.GetMaxZExtent()
		  ||zMax<pVoxelLimit.GetMinZExtent())
		{
		  return false;
		}
	      else
		{
		  if (zMin<pVoxelLimit.GetMinZExtent())
		    {
		      zMin=pVoxelLimit.GetMinZExtent();
		    }
		  if (zMax>pVoxelLimit.GetMaxZExtent())
		    {
		      zMax=pVoxelLimit.GetMaxZExtent();
		    }
		}
	    }

// Known to cut cylinder
	  switch (pAxis)
		{
		case kXAxis:
		  //G4cout << "kXAxis" << endl;
		  yoff1=yoffset-yMin;
		  yoff2=yMax-yoffset;
		  if (yoff1>=0&&yoff2>=0)
		    {
		      // Y limits cross max/min x => no change
		      pMin=xMin;
		      pMax=xMax;
		    }
		  else
		    {
		      // Y limits don't cross max/min x => compute max delta x, hence new mins/maxs
		      diff1=sqrt(abs(endOuterRadius2-yoff1*yoff1));
		      diff2=sqrt(abs(endOuterRadius2-yoff2*yoff2));
		      maxDiff=(diff1>diff2) ? diff1:diff2;
		      newMin=xoffset-maxDiff;
		      newMax=xoffset+maxDiff;
		      pMin=(newMin<xMin) ? xMin : newMin;
		      pMax=(newMax>xMax) ? xMax : newMax;
		    }
		  
		  break;
		case kYAxis:
		  //G4cout << "kYAxis" << endl;
		  xoff1=xoffset-xMin;
		  xoff2=xMax-xoffset;
		  if (xoff1>=0&&xoff2>=0)
		   {
		      // X limits cross max/min y => no change
		     //G4cout << "assegnazione diretta, xoff1,2 > 0" << yMin << " " << yMax << endl;
		     pMin=yMin;
		     pMax=yMax;
		   }
		  else
		    {
		      // X limits don't cross max/min y => compute max delta y, hence new mins/maxs
		      diff1=sqrt(abs(endOuterRadius2-xoff1*xoff1));
		      diff2=sqrt(abs(endOuterRadius2-xoff2*xoff2));
		      //G4cout << diff1 << " " << diff2 << endl;
		      maxDiff=(diff1>diff2) ? diff1:diff2;
		      //G4cout << "maxDiff " << maxDiff << endl;
		      newMin=yoffset-maxDiff;
		      newMax=yoffset+maxDiff;
		      //G4cout << "assegnazione indiretta" << yMin << " " << yMax << endl;
		      //G4cout << "newMin, newMax : " << newMin << " " << newMax << endl;
		      
		      pMin=(newMin<yMin) ? yMin : newMin;
		      pMax=(newMax>yMax) ? yMax : newMax;
		    }
		  break;
		case kZAxis:
		  //G4cout << "kZAxis" << endl;
		  pMin=zMin;
		  pMax=zMax;
		  break;
		}
	  
	  pMin-=kCarTolerance;
	  pMax+=kCarTolerance;

	  //G4cout << xoffset << " " << yoffset << " " << zoffset << endl;
	  //G4cout << "pMin, pMax: " << pMin << "," << pMax << endl << endl;
	  return true;
	  
	}
  else
    {
      G4int i,noEntries,noBetweenSections4;
      G4bool existsAfterClip=false;
      
      // Calculate rotated vertex coordinates
      G4ThreeVectorList *vertices;
      
      vertices=CreateRotatedVertices(pTransform);
      
      pMin=+kInfinity;
      pMax=-kInfinity;
      
      noEntries=vertices->entries();
      noBetweenSections4=noEntries-4;
      
      for (i=0;i<noEntries;i+=4)
	{
	  ClipCrossSection(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
	}
      
      for (i=0;i<noBetweenSections4;i+=4)
	{
	  ClipBetweenSections(vertices,i,pVoxelLimit,pAxis,pMin,pMax);
	}
      
      if (pMin!=kInfinity||pMax!=-kInfinity)
	{
	  existsAfterClip=true;
	  
	  // Add 2*tolerance to avoid precision troubles
	  pMin-=kCarTolerance;
	  pMax+=kCarTolerance;
	  
	}
      else
	{
	  // Check for case where completely enveloping clipping volume
	  // If point inside then we are confident that the solid completely
	  // envelopes the clipping volume. Hence set min/max extents according
	  // to clipping volume extents along the specified axis.
	  G4ThreeVector clipCentre(
				   (pVoxelLimit.GetMinXExtent()+pVoxelLimit.GetMaxXExtent())*0.5,
				   (pVoxelLimit.GetMinYExtent()+pVoxelLimit.GetMaxYExtent())*0.5,
				   (pVoxelLimit.GetMinZExtent()+pVoxelLimit.GetMaxZExtent())*0.5);
	  
	  if (Inside(pTransform.Inverse().TransformPoint(clipCentre))!=kOutside)
	    {
	      existsAfterClip=true;
	      pMin=pVoxelLimit.GetMinExtent(pAxis);
	      pMax=pVoxelLimit.GetMaxExtent(pAxis);
	    }
	}
      delete vertices;
      return existsAfterClip;
    }
}

// Decides whether point is inside,outside or on the surface
EInside G4Hype::Inside(const G4ThreeVector& p) const
{
  // Get third component of point "in study"
  double xZ=abs(p.z());
  //G4cout << "Inside::This point is : " << p << endl;
  //G4cout << "Surface data : " << outerRadius << endl;

  // if point's Z component is greater than halfLenZ it is SURELY outside.
  if (xZ<=(1+SURFACE_PRECISION)*halfLenZ)
	{
	  // for performance reasons I work on squared quantities
	  double xR2=p.x()*p.x()+p.y()*p.y();

	  // outer radius value at xZ
	  double oRad2=HypeOuterRadius2(xZ);
	  // inner radius value at xZ
	  double iRad2=HypeInnerRadius2(xZ);

	  // two cases:
	  // 1. w/o inner surface 
	  // 2. w/  inner surface

	  //G4cout << "iRad, xR, oRad : " << sqrt(iRad2) << "  " << sqrt(xR2) << " " << sqrt(oRad2) << endl;

	  if ((innerRadius==0.) && (innerStereo==0.))
	    {
	      // Case 1. w/o inner surface
	      if (xR2<oRad2)
		{ // on endcaps ?
		  if (abs(xZ/halfLenZ-1)<SURFACE_PRECISION) 
		    {
		      //G4cout << "kSurf, no Inner, endcaps" << endl;
		      return kSurface; // yes. On endcap!
		    }
		  else 
		    {
		      //G4cout << "kInside, no Inner, endcaps" << endl;		      
		      return kInside;  // no. It's inside!
		    }
		}
	      // is it on the hyperbolical surfaces ?
	      else
		if (abs(sqrt(xR2/oRad2)-1) < SURFACE_PRECISION)
		  {
		    //G4cout << "kSurf, no Inner, hype surface" << endl;
		    return kSurface;  // yes...it's on the surface!
		  }
		else 
		  {
		    //G4cout << "kOut, no Inner, hype surface" << endl;
		    return kOutside;  // no...outside the hype.
		  }
	    }
	  else
	    {
	      // Case 2. W/ inner surface
	      
	      // is it between the hyperbolical surfaces ? 
	      if (xR2>=iRad2 && xR2<=oRad2) 
		{
		  // on endcaps ??
		  if (abs(xZ/halfLenZ-1)<SURFACE_PRECISION) 
		    {
		      //G4cout << "kSurf, Inner, endcaps" << endl;
		      return kSurface; // yes. On endcap!
		    }
		  else 
		    {
		      //G4cout << "kInside, inner, endcaps" << endl;
		      return kInside;  // no. It's inside!
		    }
		}
	      // is it on the hyperbolical surfaces ?
	      else
		if (abs(sqrt(xR2/oRad2)-1) < SURFACE_PRECISION
		    || 
		    abs(sqrt(xR2/iRad2)-1) < SURFACE_PRECISION) 
		  {
		    //G4cout << "kSurf, inner, hype surface" << endl;
		    return kSurface;  // yes...it's on the surface!
		  }
		else 
		  {
		    //G4cout << "kOutside, Inner, hype surface" << endl;
		    return kOutside;  // no...outside the hype.
		  }
	    }
	}
  //G4cout << "qui ci devo arrivare sse abs(z)>halfLenZ!!!! " << endl;
  return kOutside;
}

// return the normal unit vector to the Hyperbolical Surface at a point 
// p on (or nearly on) the surface
G4ThreeVector G4Hype::SurfaceNormal( const G4ThreeVector& p) const
{
  G4ThreeVector norm(0.,0.,0.);
  double pZ=p.z();
  double zRadius2=p.x()*p.x()+p.y()*p.y();
  double param=0;

  //G4cout << "Surface Normal " << p << endl;
  //G4cout << "Surface data : " << outerRadius << endl;


  // check if the point is on the endcaps
  // within SURFACE_PRECISION
  if (abs((pZ/halfLenZ)-1)<SURFACE_PRECISION)
    {
      if (abs(sqrt(zRadius2/endInnerRadius2)-1)<SURFACE_PRECISION 
	  && abs(sqrt(zRadius2/endOuterRadius2)-1)<SURFACE_PRECISION)
	{
	  // this way the versor has the right direction to exit the surface
	  norm=G4ThreeVector(0.,0.,(abs(p.z())/p.z()));
	}
    }
  else
    {
      // normal to hyperbolical surfaces
      // first I have to choose which surface to use
      // I use SURFACE_PRECISION * halfLenZ as a scale factor for precision 

      // Inner Surface ?
      if (abs(sqrt(zRadius2/HypeInnerRadius2(pZ))-1) < SURFACE_PRECISION) 
	{
	  norm= G4ThreeVector(p.x(),p.y(),-p.z()*tanInnerStereo2);
	}
      
      // Outer Surface ?
      if (abs(sqrt(zRadius2/HypeOuterRadius2(p.z()))-1) < SURFACE_PRECISION) 
	{
	  norm= G4ThreeVector(p.x(),p.y(),-p.z()*tanOuterStereo2);
	}
    }
  // return unit vector parallel to the calculated one
  return norm.unit();
}


// Calculate distance to shape from outside, along normalised vector
// - return kInfinity if no intersection, or intersection distance <= tolerance

G4double G4Hype::DistanceToIn(const G4ThreeVector& p,
			      const G4ThreeVector& v    ) const
{
    G4double snxt=kInfinity;			// snxt = default return value
    EInside Pp=Inside(p);
    double zero_tolerance=1e-6;

    //G4cout << "DistToIn(pt,dir)" << p << " " << v << endl;
    //G4cout << "Surface data : " << outerRadius << endl;

    if (Pp==kSurface || Pp==kInside)
      {
	return 0.; 
      }
    // if the point is inside the Distance to enter the surface is zero
    else
      {

	// Returns distance or kInfinity if there exists or not  an intersection!
	// Set intersection point to an arbitrarily far point

        // Calculate the distance from a point (p)
        // to enter or exit the Hype along a certain direction (v)
        G4ThreeVector direction=v;
        G4ThreeVector workPoint=p;

	// init various distance parameters
	double distanceLeft =kInfinity;
	double distanceRight=kInfinity;
	double distanceInner=kInfinity;
	double distanceOuter=kInfinity;

        // make direction a versor
	direction=direction.unit();

        // line parameters - u vector
        double ux=direction.x();
        double uy=direction.y();
        double uz=direction.z();

        // point parameters - v point
        double vx=workPoint.x();
        double vy=workPoint.y();
        double vz=workPoint.z();

	// ==============================================================================//
	// Calculate intersection point with endcaps

	// uz must be, within tolerance, NOT zero
	if (abs(uz)>=zero_tolerance) 
	  {
	    // *left* (see docs) intersection z=-halfLenZ
	    double tLeft=(-halfLenZ-vz)/uz;
	    if (tLeft>0 || abs(tLeft)<zero_tolerance) 
	      {
		double xIL=ux*tLeft+vx;
		double yIL=uy*tLeft+vy;
		double radiusL=sqrt(xIL*xIL+yIL*yIL);
		// is it a real intersection ? (within tolerance)
		// is the radius compatible with endcaps radii ? 
		if ((radiusL>=(1-SURFACE_PRECISION)*endInnerRadius)
		    && radiusL<=(1+SURFACE_PRECISION)*endOuterRadius) 
		  {
		    distanceLeft=abs(tLeft);
		  }	    
	      }

	    // *right* (see docs) intersection z=halfLenZ
	    double tRight=(halfLenZ-vz)/uz;
	    if (tRight>0 || abs(tRight)<zero_tolerance)
	      {
		double xIR=ux*tRight+vx;
		double yIR=uy*tRight+vy;
		double radiusR=sqrt(xIR*xIR+yIR*yIR);
		// is it a real intersection ? 
		// is the radius compatible with endcaps radii ?   
		if (radiusR>=(1-SURFACE_PRECISION)*endInnerRadius
		    && radiusR<=(1+SURFACE_PRECISION)*endOuterRadius) 
		  {
		    distanceRight=abs(tRight);
		  }
	      }
	  }
	else
	  {
	    if (abs(abs(vz/halfLenZ)-1)<SURFACE_PRECISION) // ON endcaps
	      {
		double xR=sqrt(vx*vx+vy*vy);
		if (xR<=(1+SURFACE_PRECISION)*innerRadius) // inside
		  {
		    if (vz<0) // left EndCap
		      { 
			distanceLeft=abs(xR-innerRadius);
		      }
		    else  // right EndCap
		      {
			distanceRight=abs(xR-innerRadius);
		      }
		  }

		if (xR>=(1-SURFACE_PRECISION)*outerRadius) // outside
		  {
		    if (vz<0) // left EndCap
		      { 
			distanceLeft=abs(xR-outerRadius);
		      }
		    else // right EndCap
		      {
			distanceRight=abs(xR-outerRadius);
		      }
		  }
	      }
	  }      

	// ==============================================================================//

	// Calculate intersections with hyperbolical surfaces

	// Inner first! ====================================
        // equation parameters

	double dist1=kInfinity;
	double dist2=kInfinity;
	double param1, param2;

	if (innerRadius!=0. || innerStereo!=0.) // is there an inner surface ???
	  {
	    // yes...inner & outer are both present!
	    // degenerations are automagically cured (i.e. hyperbolical -> cylindrical)
	    param1=tanInnerStereo2;
	    param2=innerRadius2;

	    double paramA=(ux*ux+uy*uy)-param1*uz*uz;
	    double paramB=(ux*vx+uy*vy)-param1*uz*vz;
	    double paramC=(vx*vx+vy*vy)-param1*vz*vz-param2;
 
	    // intersection equation discriminant
	    double delta=paramB*paramB-paramA*paramC;		  

	    if (delta>0. || abs(delta)<zero_tolerance)
	      { // solution exists!!
		delta=abs(delta); 
		double solution1=(-paramB+sqrt(delta))/paramA;
		double solution2=(-paramB-sqrt(delta))/paramA;
		
		if (solution1>0. || abs(solution1)<zero_tolerance) 
		  {
		    dist1=abs(solution1);
		  }

		if (solution2>0. || abs(solution2)<zero_tolerance)
		  {
		    dist2=abs(solution2);
		  }
		
		// check if they fall in the correct length range

		// if BOTH intersection are real (geometric cut)
		if ((abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		    && (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)) 
		  {
		    distanceInner=((dist1 < dist2) ? dist1 : dist2);
		  }
		else
		  {
		    if (abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		      {
			distanceInner=dist1;
		      }
		    if (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		      {
			distanceInner=dist2;
		      }		     
		  }
	      }
	  }

      	// Inner Done! ==================	

	// and now Outer! ===============================================
	dist1=kInfinity;
	dist2=kInfinity;

	// look for degenerations (i.e. hyperbolical -> cylindrical)
	param1=tanOuterStereo2;
	param2=outerRadius2;

	// equation parameters
	double paramA=(ux*ux+uy*uy)-param1*uz*uz;
	double paramB=(ux*vx+uy*vy)-param1*uz*vz;
	double paramC=(vx*vx+vy*vy)-param1*vz*vz-param2;

        // intersection equation discriminant
        double delta=paramB*paramB-paramA*paramC;
	
        if (delta > 0 || abs(delta)<zero_tolerance)
	  { // solution exists!!
	    delta=abs(delta);
	    double solution1=(-paramB+sqrt(delta))/paramA;
	    double solution2=(-paramB-sqrt(delta))/paramA;
	    
	    if (solution1>0 || abs(solution1)<zero_tolerance) 
	      {
		dist1=abs(solution1);
	      }
	    
	    if (solution2>0 || abs(solution2)<zero_tolerance)
	      {
		dist2=abs(solution2);
	      }
	    
	    // check if they fall in the correct length range		    
	    // if BOTH intersection are real (geometric cut)
	    if ((abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		&& (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ))
	      {
		    distanceOuter=((dist1 < dist2) ? dist1 : dist2);
	      }
	    else
	      {
		if (abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		  distanceOuter=dist1;
		if (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		  distanceOuter=dist2;
	      }		  
	  }
	// Outer Done! ===============================================

	// and now choose che correct snxt value
	if (snxt > distanceInner)
	  {
	    snxt=distanceInner;
	  }
	if (snxt > distanceOuter) 
	  {
	    snxt=distanceOuter;
	  }
	if (snxt > distanceRight) 
	  {
	    snxt=distanceRight;
	  }
	if (snxt > distanceLeft ) 
	  {
	    snxt=distanceLeft;
	  }
      }
    // Done!
    //G4cout << "DistToIn(pt,dir) = " << snxt << endl;
    return snxt;
}
 

// Calculate distance to shape from outside, along perpendicular direction 
// (if one exists)

G4double G4Hype::DistanceToIn(const G4ThreeVector& p) const
{
  double zero_tolerance=1e-6;
  // I use the cylindrical simmetry of the system to optimize calculation
  double pZ=abs(p.z());
  double pR2=p.x()*p.x()+p.y()*p.y();
  double pR=sqrt(pR2);
  double ApointR=pR;
  double ApointZ=pZ;
  double BpointR=0.;
  double BpointZ=pZ;
  
  //G4cout << "DistToIn(pt)" << p << endl;
  //G4cout << "Surface data : " << outerRadius << endl;

  EInside Pp=Inside(p);
  if (Pp==kSurface || Pp==kInside) 
    {
      //G4cout << "Point on surface or inside, returning 0" << endl;
      return 0.;
    }

  // is point in zone 1 ?  
  if (pZ>=(1-SURFACE_PRECISION)*halfLenZ && 
      pR2<=(1+SURFACE_PRECISION)*(1+SURFACE_PRECISION)*endOuterRadius2 && 
      pR2>=(1-SURFACE_PRECISION)*(1-SURFACE_PRECISION)*endInnerRadius2)
    { 
      //G4cout << "zone1 - returning : " <<  abs(pZ)-halfLenZ << endl;
      return abs(pZ)-halfLenZ; 
    }

  // special degenerated case: point with z=0
  if (pZ<zero_tolerance) 
    { 
      if (pR<=(1+SURFACE_PRECISION)*innerRadius) 
	{
	  //G4cout << "z=0, inner nearer,  returning: " << (innerRadius-pR) << endl;
	  return abs(innerRadius-pR);
	}
      else
	{
	  //G4cout << "z=0, outer nearer, returning: " << pR2-outerRadius << endl;
	  return abs(pR-outerRadius);
	}
    }

  // check for degenerations (hyperbolical -> cylindrical)
  if (outerStereo==0.)
    {
      // cylindrical
      if (pR>=(1-SURFACE_PRECISION)*outerRadius &&
	  pZ<=(1+SURFACE_PRECISION)*halfLenZ)
	{
	  //G4cout << "zone2 - Cyl - returning : " <<  pR-outerRadius << endl;
	  return abs(pR-outerRadius);
	}
    }
  else // hyperbolical!
    {
      // some parameters I need for zone 2 and 3
      // zone2: normal has to be calculated in (halfLenZ, endOuterRadius)
      // zone3: normal has to be calculater in (halfLenZ, endInnerRadius)

      // Parameters for Zone2
      // since I use just normR/normZ and normZ/normR I simplify the 2 factor
      // maxNormR=extremeR;
      // maxNormZ=-tanThetaStereo*tanThetaStereo*halfLenZ;
  
      // calculate parameters:
      // paramNorm=maxNormR/maxNormZ
      double paramNormZ2=endOuterRadius/(-tanOuterStereo2*halfLenZ);

      // maxRz2 identify the zone on the theta plane where 
      // I have to use the iterative method to calculate the distance
      // this is an "hyperbolical" triangle with vertex:
      // 1 -> (0,outerRadius)
      // 2 -> (0,maxR);
      // 3 -> (halfLenZ,0);
      // Ok - Verified on 20.02.98 
      
      // maxZ_Z2=-(1/paramNorm)*extremeR+halfLenZ;
      // maxRz2=-paramNormZ2*halfLenZ+endOuterRadius;
      
      // is point in zone 2 ?
    
      if (pR>=(1-SURFACE_PRECISION)*sqrt(HypeOuterRadius2(pZ)) // above the outer hype profile
	  && ApointR<=(paramNormZ2*(ApointZ-halfLenZ)+endOuterRadius)
	  // but under the normal in the extreme point
	  && ApointZ<=(1+SURFACE_PRECISION)*halfLenZ)
	{
	  //Prepare the parameter for the iteration
	  double mydist=0.;       // distance between C and A
	  double backupdist=kInfinity;  // old distance...
	  // needed to check if I'm moving in the right direction
	  double tanThetaStereoSquared=tanOuterStereo2;
	  double param2=1+2*tanThetaStereoSquared;
	  
	  do {
	    // move along the curve
	    // at the first Step distance=0., so there's no problem
	    // then mydist is a value *with* Sign
	    BpointZ=BpointZ+mydist;
	    
	    // calculate the normal to the Hype in
	    // B=(BpointZ, BpointR)
	    // (this is the intersection point between the Hype and the 
	    // line normal to the Hype axis passing in A)
	    BpointR=sqrt(HypeOuterRadius2(BpointZ));
	    
	    // //G4cout << "BpointZ : " << BpointZ << "  BpointR: " << BpointR << endl; 
	    
	    // the components of the normal vector in B are
	    // / normR = 2*BpointR
	    // \ normZ =-2*tanThetaStereoSquared*BpointZ
	    // but, since I'm interested only in the normZ/normR value,
	    // I can simplify the 2 factor
	    // with a bit of calculation (see docs) I get:
	    
	    mydist=BpointZ*(param2-tanThetaStereoSquared*ApointR/BpointR)-ApointZ;
	    
	    // if mydist grows during the iteration, it means I'm 
	    // moving in the wrong direction. So I have to 
	    // 1. revert back to the place where I started from
	    //    i.e. stepping back for the  backupdist
	    // 2. move of backupdist in the opposite direction
	    // so mydist is the double of backupdist with opposite Sign
	    // //G4cout << " mydist : " << mydist << endl;
	    
	    if (abs(mydist)<abs(backupdist)) { backupdist=mydist; }
	    else { mydist=-2*backupdist; }
	    
	  } while (abs(mydist)>zero_tolerance);
	  
	  // Now I have the correct point BpointZ, BpointR stored...
	  //G4cout << "zone2 - hyp - returning : " << 
	  //sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
	  // (ApointR-BpointR)*(ApointR-BpointR)) << endl;

	  return sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
		      (ApointR-BpointR)*(ApointR-BpointR));
	}
    }

  if (innerStereo==0.)
    {
      // cylindrical
      if (pR<=(1+SURFACE_PRECISION)*innerRadius &&
	  pZ<=(1+SURFACE_PRECISION)*halfLenZ)
	{
	  //G4cout << "zone3 - Cyl - returning : " << innerRadius-pR << endl;
	  return abs(innerRadius-pR);
	}
    }
  else // hyperbolical!  
    {
      // Parameters for Zone 3
      // it innerRadius=0 and innerStereo=0 then NO zone 3
      if (innerRadius!=0.)
	{
	  // calculate parameters:
	  // paramNorm=maxNormR/maxNormZ
	  double coParamNormZ3=(-tanInnerStereo2*halfLenZ)/endInnerRadius;
	  
	  // maxRz2 identify the zone on the theta plane where 
	  // I have to use the iterative method to calculate the distance
	  // this is an "hyperbolical" triangle with vertex:
	  // 1 -> (0,innerRadius)
	  // 2 -> (maxZ,0);
	  // 3 -> (halfLenZ,endInnerRadius);
	  
	  double maxZz3=-(coParamNormZ3)*endInnerRadius+halfLenZ;
	  // maxRz3=-paramNormZ3*halfLenZ+endInnerRadius;
	  
	  // is point in zone 3 ?
	  
	  if (pR<=(1+SURFACE_PRECISION)*sqrt(HypeInnerRadius2(pZ)) // under the inner hype profile
	      && ApointZ<=(1+SURFACE_PRECISION)*((coParamNormZ3*(ApointR-endOuterRadius)+ halfLenZ)))
	    {
	      //Prepare the parameter for the iteration
	      double mydist=0.;       // distance between C and A
	      double backupdist=kInfinity;  // old distance...
	      // needed to check if I'm moving in the right direction
	      double tanThetaStereoSquared=tanInnerStereo2;
	      double param2=1+2*tanThetaStereoSquared;
	      
	      do {
		// move along the curve
		// at the first Step distance=0., so there's no problem
		// then mydist is a value *with* Sign
		BpointZ=BpointZ+mydist;
		
		// calculate the normal to the Hype in
		// B=(BpointZ, BpointR)
		// (this is the intersection point between the Hype and the 
		// line normal to the Hype axis passing in A)
		BpointR=sqrt(HypeInnerRadius2(BpointZ));
		
		// //G4cout << "BpointZ : " << BpointZ << "  BpointR: " << BpointR << endl; 
		
		// the components of the normal vector in B are
		// / normR = 2*BpointR
		// \ normZ =-2*tanThetaStereoSquared*BpointZ
		// but, since I'm interested only in the normZ/normR value,
		// I can simplify the 2 factor
		// with a bit of calculation (see docs) I get:
		
		mydist=BpointZ*(param2-tanThetaStereoSquared*ApointR/BpointR)-ApointZ;
		
		// if mydist grows during the iteration, it means I'm 
		// moving in the wrong direction. So I have to 
		// 1. revert back to the place where I started from
		//    i.e. stepping back for the  backupdist
		// 2. move of backupdist in the opposite direction
		// so mydist is the double of backupdist with opposite Sign
		// //G4cout << " mydist : " << mydist << endl;
		
		if (abs(mydist)<abs(backupdist)) { backupdist=mydist; }
		else { mydist=-2*backupdist; }
		
	      } while (abs(mydist) > zero_tolerance);
	      
	      // Now I have the correct point BpointZ, BpointR stored...
	      //G4cout << " zone3 - hyp - returning: " <<
	      //sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
	      //     (ApointR-BpointR)*(ApointR-BpointR)) << endl;
	      return sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
			  (ApointR-BpointR)*(ApointR-BpointR));
	    }
	}
    }
  // zone 4a or 4b 
  if (pR<=(1+SURFACE_PRECISION)*endInnerRadius) 
    { // 4b
      //G4cout << " zone 4b - returning: " <<
      // sqrt((ApointR-endInnerRadius)*(ApointR-endInnerRadius)+
      //      (ApointZ-halfLenZ)*(ApointZ-halfLenZ)) << endl;

      return sqrt( 
		  (ApointR-endInnerRadius)*(ApointR-endInnerRadius)+
		  (ApointZ-halfLenZ)*(ApointZ-halfLenZ));
    }
  else
    { // 4a
      //G4cout << " zone 4a - returning: " <<
      //sqrt((ApointR-endOuterRadius)*(ApointR-endOuterRadius)+
      //     (ApointZ-halfLenZ)*(ApointZ-halfLenZ)) << endl;

      return sqrt( 
		  (ApointR-endOuterRadius)*(ApointR-endOuterRadius)+
		  (ApointZ-halfLenZ)*(ApointZ-halfLenZ));
    }
}



// Calculate distance to surface of shape from `inside', allowing for tolerance
G4double G4Hype::DistanceToOut(const G4ThreeVector& p,const G4ThreeVector& v,
			       const G4bool calcNorm, G4bool *validNorm,G4ThreeVector *n) const
{
    G4double snxt=kInfinity;			// snxt = default return value
    ESide side;
    double xi,yi,zi;  // service parameters
    double zero_tolerance=1E-6;

    if (calcNorm) *validNorm=true;  // All normals are valid
    *n=G4ThreeVector(0.,0.,0.);
    //G4cout << "DistanceToOut(pt,dir): " << p << " " << v << endl;
    //G4cout << "Surface Data " << outerRadius << endl;

    EInside Pp=Inside(p);
    if (Pp==kOutside)
      {
	snxt=0.;
      }
    else
      {
	// Returns distance or kInfinity if there exists or not  an intersection!

	// Calculate the distance from a point (p)
        // to exit the Hype along a certain direction (v)
        G4ThreeVector direction=v;
        G4ThreeVector workPoint=p;

	// init various distance parameters
	double distanceLeft= kInfinity;
	double distanceRight=kInfinity;
	double distanceInner=kInfinity;
	double distanceOuter=kInfinity;

        // make direction a versor
	direction=direction.unit();
	//G4cout << "direzione : " << direction << endl;
	//G4cout << "workPoint : " << workPoint << endl;

        // line parameters - u vector
        double ux=direction.x();
        double uy=direction.y();
        double uz=direction.z();

        // point parameters - v point
        double vx=workPoint.x();
        double vy=workPoint.y();
        double vz=workPoint.z();

	// ==============================================================================//
	// Calculate intersection point with endcaps
	// First I have to check that uz is not within tolerance with zero!

	if (abs(uz)>=zero_tolerance) 
	  {
	    // *left* (see docs) intersection z=-halfLenZ
	    double tLeft=(-halfLenZ-vz)/uz;
	    //   //G4cout << tLeft << endl;
	    if (tLeft>0 || abs(tLeft)<zero_tolerance) 
	      {
		double xIL=ux*tLeft+vx;
		double yIL=uy*tLeft+vy;
		double radiusL=sqrt(xIL*xIL+yIL*yIL);
		
		// is it a real intersection ? 
		if (radiusL>=(1-SURFACE_PRECISION)*endInnerRadius &&
		    radiusL<=(1+SURFACE_PRECISION)*endOuterRadius)
		  {
		    distanceLeft=tLeft;
		  }
	      }
	
	    // *right* (see docs) intersection z=halfLenZ
	    double tRight=(halfLenZ-vz)/uz;
	    // //G4cout << tRight << endl;
	    if (tRight>0 || abs(tRight)<zero_tolerance)
	      {
		double xIR=ux*tRight+vx;
		double yIR=uy*tRight+vy;
		double radiusR=sqrt(xIR*xIR+yIR*yIR);

		// is it a real intersection ? 
		if (radiusR>=(1-SURFACE_PRECISION)*endInnerRadius && 
		    radiusR<=(1+SURFACE_PRECISION)*endOuterRadius) 
		  {
		    distanceRight=abs(tRight);
		  }
	      }
	    	    
	    //G4cout << "tLeft, tRight : " << tLeft << "  " << tRight << endl;
	  }
	  
	// ==============================================================================//

	// Calculate intersections with hyperbolical surfaces

	double param1,param2;

	// Inner first! ====================================
	if (innerRadius!=0. || innerStereo !=0.)
	  {
	    param1=tanInnerStereo2;
	    param2=innerRadius2;
	    
	    double paramA=(ux*ux+uy*uy)-param1*uz*uz;
	    double paramB=(ux*vx+uy*vy)-param1*uz*vz;
	    double paramC=(vx*vx+vy*vy)-param1*vz*vz-param2;
	    
	    // intersection equation discriminant
	    double delta=paramB*paramB-paramA*paramC;

	    if (delta>0 || abs(delta)<zero_tolerance)
	      { 
		// solution exists!!
		delta=abs(delta);
		double solution1=(-paramB+sqrt(delta))/paramA;
		double solution2=(-paramB-sqrt(delta))/paramA;
		double dist1=kInfinity;
		double dist2=kInfinity;
		//G4cout << "Sol 1,2 " << solution1 << "  " << solution2 << endl;

		if (solution1>0 || abs(solution1)<zero_tolerance)
		  dist1=abs(solution1);

		if (solution2>0 || abs(solution2)<zero_tolerance)
		  dist2=abs(solution2);

		//if (abs(solution1)<=zero_tolerance)
		//  dist1=0.


		// check if they fall in the correct length range
		if ((abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ) && 
		    (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)) 
		  {
		    distanceInner=((dist1 < dist2) ? dist1 : dist2);
		  }
		else
		  {
		    if (abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		      distanceInner=dist1;
		    if (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		      distanceInner=dist2;
		  }
	      }
	  }
      	// Inner Done! ==================	

	// and now Outer! ===============================================

	// equation parameters
	double dist1=kInfinity;
	double dist2=kInfinity;

	// look for degenerations (i.e. hyperbolical -> cylindrical)
	param1=tanOuterStereo2;
	param2=outerRadius2;
	
	double paramA=(ux*ux+uy*uy)-param1*uz*uz;
	double paramB=(ux*vx+uy*vy)-param1*uz*vz;
	double paramC=(vx*vx+vy*vy)-param1*vz*vz-param2;

	// intersection equation discriminant
        double delta=paramB*paramB-paramA*paramC;
	//G4cout << "delta : " << delta << endl;
	//G4cout << "A, B, C : " << paramA << "  " << paramB << "   " << paramC << endl;

        if (delta > 0 || abs(delta)<zero_tolerance) 
        { 
	  // solution exists!!
	  delta=abs(delta);
	  double solution1=(-paramB+sqrt(delta))/paramA;
	  double solution2=(-paramB-sqrt(delta))/paramA;

	  //G4cout << "Solution 1, 2  : " << solution1*1E+14 << "   " << solution2*1E+14 << endl;
	  
	  if (solution1>0 || abs(solution1)<zero_tolerance) 
	    {
	      dist1=abs(solution1);
	    }

	  if (solution2>0 || abs(solution2)<zero_tolerance)
	    {
	      dist2=abs(solution2);
	    }

	  // check if they fall in the correct length range
	  if ((abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ) &&
	      (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)) 
	    {
	      distanceOuter=((dist1 < dist2) ? dist1 : dist2);
	    }
	  else
	    {
	      if (abs(dist1*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		distanceOuter=dist1;
	      if (abs(dist2*uz+vz)<=(1+SURFACE_PRECISION)*halfLenZ)
		distanceOuter=dist2;
	    }
        }
	// Outer Done! ===============================================

	// and now choose che correct snxt value
	
	if (snxt > distanceInner) 
	  { 
	    snxt=distanceInner;
	    side=innerFace;
	  }
	if (snxt > distanceOuter) 
	  {
	    snxt=distanceOuter;
	    side=outerFace;
	  }

	if (abs(distanceInner)<zero_tolerance && distanceOuter>0 && distanceOuter<kInfinity)
	  {
	    snxt=distanceOuter;
	    side=outerFace;
	  }

	if (abs(distanceOuter)<zero_tolerance && distanceInner>0 && distanceInner<kInfinity)
	  {
	    snxt=distanceInner;
	    side=innerFace;
	  }

	if (snxt > distanceRight) 
	  {
	    snxt=distanceRight;
	    side=rightCap;
	  }
	if (snxt > distanceLeft)
	  {
	    snxt=distanceLeft;
	    side=leftCap;
	  }

	// Done!
	//G4cout << "distanze left, right, inner, outer :" << distanceLeft << "   " 
	//     << distanceRight << "   " << distanceInner << "    " << distanceOuter << endl;
	
	// normal calculation, where required.
	// normal to EXIT the volume...
	if (calcNorm)
	  {
	    switch(side)
	      {
	      case outerFace:
		xi=vx+snxt*ux;
		yi=vy+snxt*uy;
		zi=vz+snxt*uz;
		*n=G4ThreeVector(xi,yi,-zi*tanOuterStereo2).unit();
		//G4cout << "outerFace : " << *n << endl;
		break;
	  
	      case innerFace:
		xi=vx+snxt*ux;
		yi=vy+snxt*uy;
		zi=vx+snxt*uz;
		*n=-G4ThreeVector(xi,yi,-zi*tanInnerStereo2).unit();
		//G4cout << "innerFace: " << *n << endl;
		break;
		
	      case leftCap:
		*n=G4ThreeVector(0,0,-1);
		//G4cout << "leftCap: " << *n << endl;
		break;
		
	      case rightCap:
		*n=G4ThreeVector(0,0,1);
		//G4cout << "rightCap: " << *n << endl;
		break;
		
	      default: *validNorm=false;
	      }
	  }
      }
    snxt+=zero_tolerance*halfLenZ; 
    //G4cout << "valore restituito: " << snxt << endl;
    //G4cout << "*n = " << *n << endl;
    return snxt;
}

// Calculate distance (<=actual) to closest surface of shape from inside
G4double G4Hype::DistanceToOut(const G4ThreeVector& p) const
{

  // I use the cylindrical simmetry of the system to optimize calculation
  double zero_tolerance=1e-6;
  double pZ=abs(p.z());
  double pR2=p.x()*p.x()+p.y()*p.y();
  double pR=sqrt(pR2);

  //G4cout << "pZ : " << pZ << " param: " << pZ/halfLenZ << " SP: "
  //      << SURFACE_PRECISION <<endl;

  // distance parameters
  double distEndCap=kInfinity;
  double distSurfac=kInfinity;

  //G4cout << "DistToOut(pt)" << p << endl;
  //G4cout << "Surface data: " << outerRadius << endl;

  EInside Pp=Inside(p);
  if (Pp==kOutside)
    {
      //G4cout << "Point outside, returning 0." << endl;
      return 0;
    }
      
  if (Pp==kSurface) 
    {
      //G4cout << "Point on surface, returning Step" << endl;
      return halfLenZ*SURFACE_PRECISION;
    }

  // point is surely inside from now on!
  // calculate distEndCap <------------------
  if (pR<=(1+SURFACE_PRECISION)*endOuterRadius && 
      pR>=(1-SURFACE_PRECISION)*endInnerRadius)
    { 
      distEndCap=abs(halfLenZ-pZ);
    }
  else
    { 
      distEndCap=sqrt((halfLenZ-pZ)*(halfLenZ-pZ)+
		      (pR-endInnerRadius)*(pR-endInnerRadius));
    }
  
  // special case: point with z=0
  if (pZ<zero_tolerance)
    { 
      double a1=kInfinity;
      if (innerRadius!=0) 
	{
	  a1= abs(pR-innerRadius);
	}

      double a2=abs(outerRadius-pR);
      distSurfac= (a1<a2?a1:a2);
      //G4cout << "a1: " << a1 << "  a2: " << a2 << "  distEndCap: " << 
      //	distEndCap << endl; 
      //G4cout << "Point with z=0, returning: " <<(distSurfac<distEndCap?distSurfac:distEndCap) << endl;
      return (distSurfac<distEndCap?distSurfac:distEndCap); // Done
    }
  
  // case w/o inner surface
  if (innerRadius==0. && innerStereo==0.) 
    {
      // no inner surface 
      if (outerStereo==0.)
	{
	  // cylindrical
	  if (pR<=(1+SURFACE_PRECISION)*outerRadius &&
	      pZ<=(1+SURFACE_PRECISION)*halfLenZ)
	    {
	      distSurfac=abs(outerRadius-pR);
	      //G4cout << "Cyl, no inner surf, returning: " << 
	      //(distSurfac<distEndCap?distSurfac:distEndCap) << endl;
	      return (distSurfac<distEndCap?distSurfac:distEndCap);
	    }
	}
      else // hyperbolical!  
	{
	  // only outer surface!! No need to check the nearer surface to the point
	  // since I use just normR/normZ and normZ/normR I simplify the 2 factor
	  // maxNormR=extremeR;
	  // maxNormZ=-tanThetaStereo*tanThetaStereo*halfLenZ;
	  
	  double ApointR=pR;
	  double ApointZ=pZ;
	  double BpointR=0.;
	  double BpointZ=pZ;
	  //Prepare the parameter for the iteration
	  double mydist=0.;             // distance between C and A
	  double backupdist=kInfinity;  // old distance...
	  
	  // needed to check if I'm moving in the right direction
	  double tanThetaStereoSquared=tanOuterStereo2;
	  double param2=1+2*tanThetaStereoSquared;
	  
	  do {
	    // move along the curve
	    // at the first Step distance=0., so there's no problem
	    // then distance is a value *with* Sign
	    BpointZ=BpointZ+mydist;
	    
	    // calculate the normal to the Hype in
	    // B=(BpointZ, BpointR)
	    // (this is the intersection point between the Hype and the 
	    // line normal to the Hype axis passing in A)
	    BpointR=sqrt(HypeOuterRadius2(BpointZ));
	    
	    // //G4cout << "BpointZ : " << BpointZ << "  BpointR: " << BpointR << endl; 
	    
	    // the components of the normal vector in B are
	    // / normR = 2*BpointR
	    // \ normZ =-2*tanThetaStereoSquared*BpointZ
	    // but, since I'm interested only in the normZ/normR value,
	    // I can simplify the 2 factor
	    // with a bit of calculation (see docs) I get:
	    
	    mydist=BpointZ*(param2-tanThetaStereoSquared*ApointR/BpointR)-ApointZ;
	    
	    // if mydist grows during the iteration, it means I'm 
	    // moving in the wrong direction. So I have to 
	    // 1. revert back to the place where I started from
	    //    i.e. stepping back for the  backupdist
	    // 2. move of backupdist in the opposite direction
	    // so mydist is the double of backupdist with opposite Sign
	    // //G4cout << " mydist : " << mydist << endl;
	    
	    if (abs(mydist)<abs(backupdist)) { backupdist=mydist; }
	    else { mydist=-2*backupdist; }
	    
	  } while (abs(mydist)>zero_tolerance);
	  
	  // Now I have the correct point BpointZ, BpointR stored...
	  //distSurfac=sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
	  //                (ApointR-BpointR)*(ApointR-BpointR));
	  //G4cout << "Hyperb, no inner surface, returning : " << 
	  //  (distSurfac<distEndCap?distSurfac:distEndCap) << endl;
	  return (distSurfac<distEndCap?distSurfac:distEndCap);
	}
    }
  // case w/ inner surface
  else // if(innerRadius==0 && innerStereo==0)
    {
      if (abs(pR-sqrt(HypeOuterRadius2(pZ))) < abs(pR-sqrt(HypeInnerRadius2(pZ))))
	{
	  // outer surface is nearer...calculate on that one
	  if (outerStereo==0.)
	    {
	      // cylindrical
	      if (pR<=(1+SURFACE_PRECISION)*outerRadius &&
		  pZ<=(1+SURFACE_PRECISION)*halfLenZ)
		{
		  distSurfac=abs(outerRadius-pR);
		  //G4cout << "Cyl, with inner surface, outer is nearer,  returning : " << 
		  //  (distSurfac<distEndCap?distSurfac:distEndCap) << endl;
		  return (distSurfac<distEndCap?distSurfac:distEndCap);
		}
	    }
	  else // hyperbolical!  
	    {
	      // since I use just normR/normZ and normZ/normR I simplify the 2 factor
	      // maxNormR=extremeR;
	      // maxNormZ=-tanThetaStereo*tanThetaStereo*halfLenZ;
	      
	      double ApointR=pR;
	      double ApointZ=pZ;
	      double BpointR=0.;
	      double BpointZ=pZ;
	      //Prepare the parameter for the iteration
	      double mydist=0.;       // distance between C and A
	      double backupdist=kInfinity;  // old distance...
	      
	      // needed to check if I'm moving in the right direction
	      double tanThetaStereoSquared=tanOuterStereo2;
	      double param2=1+2*tanThetaStereoSquared;
	      
	      do {
		// move along the curve
		// at the first Step distance=0., so there's no problem
		// then distance is a value *with* Sign
		BpointZ=BpointZ+mydist;
		
		// calculate the normal to the Hype in
		// B=(BpointZ, BpointR)
		// (this is the intersection point between the Hype and the 
		// line normal to the Hype axis passing in A)
		BpointR=sqrt(HypeOuterRadius2(BpointZ));
		
		// //G4cout << "BpointZ : " << BpointZ << "  BpointR: " << BpointR << endl; 
		
		// the components of the normal vector in B are
		// / normR = 2*BpointR
		// \ normZ =-2*tanThetaStereoSquared*BpointZ
		// but, since I'm interested only in the normZ/normR value,
		// I can simplify the 2 factor
		// with a bit of calculation (see docs) I get:
	    
		mydist=BpointZ*(param2-tanThetaStereoSquared*ApointR/BpointR)-ApointZ;
	    
		// if mydist grows during the iteration, it means I'm 
		// moving in the wrong direction. So I have to 
		// 1. revert back to the place where I started from
		//    i.e. stepping back for the  backupdist
		// 2. move of backupdist in the opposite direction
		// so mydist is the double of backupdist with opposite Sign
		// //G4cout << " mydist : " << mydist << endl;
		
		if (abs(mydist)<abs(backupdist)) { backupdist=mydist; }
		else { mydist=-2*backupdist; }
	
	      } while (abs(mydist)>zero_tolerance);
	      
	      // Now I have the correct point BpointZ, BpointR stored...
	      distSurfac=sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
			      (ApointR-BpointR)*(ApointR-BpointR));
	      //G4cout << "Hyp, with inner surface, outer is nearer,  returning : " << 
	      //(distSurfac<distEndCap?distSurfac:distEndCap) << endl;

	      return (distSurfac<distEndCap?distSurfac:distEndCap);
	    }
	}
      else //((sqrt(pR2)-sqrt(HypeOuterRadius2(pZ))) < (sqrt(pR2)-sqrt(HypeInnerRadius2(pZ))))
	{
	  // working on inner surface
	  if (innerStereo==0.)
	    {
	      // cylindrical
	      if (pR>=(1-SURFACE_PRECISION)*innerRadius &&
		  pZ<=(1+SURFACE_PRECISION)*halfLenZ)
		{
		  distSurfac=abs(pR-innerRadius);
		  //G4cout << "Cyl, with inner surface, inner is nearer,  returning : " << 
		  //  (distSurfac<distEndCap?distSurfac:distEndCap) << endl;

		  return (distSurfac<distEndCap?distSurfac:distEndCap);
		}
	    }
	  else // hyperbolical!  
	    {
	      // inner surface is nearer, but there are blank zones
	      // i.e. zones where it is not possible to calculate distance along normal
	      
	      double ApointR=pR;
	      double ApointZ=pZ;
	      double BpointR=0.;
	      double BpointZ=pZ;
	      
	      // calculate parameters:
	      double paramNorm=endInnerRadius/(-tanInnerStereo2*halfLenZ);
	      if (ApointR<=(1+SURFACE_PRECISION)*(paramNorm*(ApointZ-halfLenZ)+endInnerRadius))
		{
		  // //G4cout << "HERE!" << endl;
		  // not in a blank zone!!!
		  //Prepare the parameter for the iteration
		  double mydist=0.;       // distance between C and A     
		  double backupdist=kInfinity;  // old distance...
		  
		  // needed to check if I'm moving in the right direction
		  double tanThetaStereoSquared=tanInnerStereo2;
		  double param2=1+2*tanThetaStereoSquared;
		  
		  do {
		    // move along the curve
		    // at the first Step distance=0., so there's no problem
		    // then distance is a value *with* Sign
		    BpointZ=BpointZ+mydist;
		    
		    // calculate the normal to the Hype in
		    // B=(BpointZ, BpointR)
		    // (this is the intersection point between the Hype and the 
		    // line normal to the Hype axis passing in A)
		    BpointR=sqrt(HypeInnerRadius2(BpointZ));
		    
		    // //G4cout << "BpointZ : " << BpointZ << "  BpointR: " << BpointR << endl; 
		    
		    // the components of the normal vector in B are
		    // / normR = 2*BpointR
		    // \ normZ =-2*tanThetaStereoSquared*BpointZ
		    // but, since I'm interested only in the normZ/normR value,
		    // I can simplify the 2 factor
		    // with a bit of calculation (see docs) I get:
		    
		    mydist=BpointZ*(param2-tanThetaStereoSquared*ApointR/BpointR)-ApointZ;
		    
		    // if mydist grows during the iteration, it means I'm 
		    // moving in the wrong direction. So I have to 
		    // 1. revert back to the place where I started from
		    //    i.e. stepping back for the  backupdist
		    // 2. move of backupdist in the opposite direction
		    // so mydist is the double of backupdist with opposite Sign
		    // //G4cout << " mydist : " << mydist << endl;
		    
		    if (abs(mydist)<abs(backupdist)) { backupdist=mydist; }
		    else { mydist=-2*backupdist; }
		    
		  } while (abs(mydist)<zero_tolerance);
		  
		  // Now I have the correct point BpointZ, BpointR stored...
		  distSurfac=sqrt((ApointZ-BpointZ)*(ApointZ-BpointZ)+
				  (ApointR-BpointR)*(ApointR-BpointR));
		  //G4cout << "Hyp, with inner surface, inner is nearer,  returning : " << 
		  //  (distSurfac<distEndCap?distSurfac:distEndCap) << endl;

		  return (distSurfac<distEndCap?distSurfac:distEndCap);
		}
	      else // in a blank zone!!!
		{ 
		  double dist1=sqrt((halfLenZ-pZ)*(halfLenZ-pZ)+
				    (pR-endInnerRadius)*(pR-endInnerRadius));
		  double dist2=sqrt((halfLenZ-pZ)*(halfLenZ-pZ)+
				    (pR-endOuterRadius)*(pR-endOuterRadius));
		  distSurfac=(dist1<dist2?dist1:dist2);
		  //G4cout << "Blank zone,  returning : " << 
		  //  (distSurfac<distEndCap?distSurfac:distEndCap) << endl;

		  return (distSurfac<distEndCap?distSurfac:distEndCap); // Done
		}
	    } // outer hyperbolical
	}
    }
  return 0.;
}

// Create a List containing the transformed vertices

G4ThreeVectorList* G4Hype::CreateRotatedVertices(const G4AffineTransform& pTransform) const
{
  G4ThreeVectorList *vertices;
  G4ThreeVector vertex0,vertex1,vertex2,vertex3;
  G4double meshAngle,meshRMax,crossAngle,cosCrossAngle,sinCrossAngle,sAngle;
  G4double rMaxX,rMaxY,rMinX,rMinY;
  G4int crossSection,noCrossSections;

  // Compute # of cross-sections necessary to mesh tube
  // not too few, not too many
  noCrossSections=(G4int)((kMinMeshSections+kMaxMeshSections)/2)+1;
  meshAngle=M_PI*2.0/(noCrossSections-1);
  meshRMax=endOuterRadius/cos(meshAngle*0.5);

  sAngle=-meshAngle*0.5;
 
  vertices=new G4ThreeVectorList(noCrossSections*4);
  if (vertices)
    {
      for (crossSection=0;crossSection<noCrossSections;crossSection++)
	{
	  // Compute coordinates of cross section at section crossSection
	  crossAngle=sAngle+crossSection*meshAngle;
	  cosCrossAngle=cos(crossAngle);
	  sinCrossAngle=sin(crossAngle);
	  
	  rMaxX=meshRMax*cosCrossAngle;
	  rMaxY=meshRMax*sinCrossAngle;
	  rMinX=endInnerRadius*cosCrossAngle;
	  rMinY=endInnerRadius*sinCrossAngle;
	  vertex0=G4ThreeVector(rMinX,rMinY,-halfLenZ);
	  vertex1=G4ThreeVector(rMaxX,rMaxY,-halfLenZ);
	  vertex2=G4ThreeVector(rMaxX,rMaxY,+halfLenZ);
	  vertex3=G4ThreeVector(rMinX,rMinY,+halfLenZ);

	  vertices->insert(pTransform.TransformPoint(vertex0));
	  vertices->insert(pTransform.TransformPoint(vertex1));
	  vertices->insert(pTransform.TransformPoint(vertex2));
	  vertices->insert(pTransform.TransformPoint(vertex3));
	}
    }
  else
    {
      G4Exception("G4Hype::CreateRotatedVertices Out of memory - Cannot alloc vertices");
    }
  return vertices;
}

void G4Hype::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}

G4VisExtent G4Hype::GetExtent() const 
{
  // Define the sides of the box into which the G4Hype instance would fit.
  return G4VisExtent (-endOuterRadius, endOuterRadius,
		      -endOuterRadius, endOuterRadius,
		      -halfLenZ, halfLenZ);
}

G4Polyhedron* G4Hype::CreatePolyhedron () const 
{
  //  return new G4PolyhedronHype (fRMin, fRMax, fDz, fSPhi, fDPhi);
  return 0;
}

G4NURBS* G4Hype::CreateNURBS () const 
{
  return new G4NURBStube(endInnerRadius, endOuterRadius, halfLenZ); // Tube for now!!!
}
