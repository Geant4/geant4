// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidPolyhedra.cc,v 1.8 1999-05-25 17:51:43 sgiani Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Corrections by S.Giani:
// - Xaxis now corresponds to phi=0
// - partial angle = phiTotal / Nsides
// - end planes exact boundary calculation for phiTotal < 2pi 
//   (also including case with RMIN=RMAX) 
// - Xaxis now properly rotated to compute correct scope of vertixes
// - corrected surface orientation for outer faces parallel to Z
// - completed explicit setting of the orientation for all faces
// - some comparison between doubles avoided by using tolerances. 
//
// 
// The polygonal solid G4BREPSolidPolyhedra is a shape defined by an inner 
// and outer polygonal surface and two planes perpendicular to the Z axis. 
// Each polygonal surface is created by linking a series of polygons created 
// at different planes perpendicular to the Z-axis. All these polygons all 
// have the same number of sides (sides) and are defined at the same Z planes 
// for both inner and outer polygonal surfaces. 
//

#include "G4BREPSolidPolyhedra.hh"
#include "G4FPlane.hh"


G4BREPSolidPolyhedra::G4BREPSolidPolyhedra(G4String name,
					   const G4double phi1,
					   const G4double dphi,
					   const int sides,
					   const int  num_z_planes,      
					   const G4double z_start,
					   const G4double z_values[],
					   const G4double RMIN[],
					   const G4double RMAX[]     
					   ) : G4BREPSolid(name)
{
  const int sections= num_z_planes - 1;
  
  if(dphi >= 2*pi-perMillion)
    nb_of_surfaces = 2*(sections * sides) + 2;
  else
    nb_of_surfaces = 2*(sections * sides) + 4;


  SurfaceVec = new G4Surface*[nb_of_surfaces];
  
  G4Vector3D Axis(0,0,1);
  G4Vector3D XAxis(1,0,0);
  G4Vector3D TmpAxis;
  G4Point3D  Origin(0,0,z_start);    
  G4Point3D  LocalOrigin(0,0,z_start);    
  G4double   Length;
  G4int      Count     = 0 ;
  G4double   PartAngle = (dphi)/sides;


  ///////////////////////////////////////////////////
   

  for(G4int a=0;a<sections;a++)
  {
    TmpAxis= XAxis;
    TmpAxis.rotateZ(phi1);
    Length = z_values[a+1] - z_values[a];

    // L. Broglia
    // Be careful in the construction of the planes
    // See G4FPlane 
    
    // Create sides
    for(int b=0;b<sides;b++)
    {
      G4Point3DVector PointList(4);
      // Create inner side
      // Calc points for the planar surface boundary
      // The order of the point give the sense 
      PointList[0] = LocalOrigin + (RMIN[a] * TmpAxis);
      PointList[3] = LocalOrigin + (Length*Axis) + (RMIN[a+1] * TmpAxis);
      TmpAxis.rotateZ(PartAngle);
      PointList[2] = LocalOrigin + (Length*Axis) + (RMIN[a+1] * TmpAxis);
      PointList[1] = LocalOrigin + (RMIN[a] * TmpAxis);   

      // Add to surface list and reverse sense	  
      SurfaceVec[Count] = new G4FPlane( &PointList, 0, 0);

      Count++;
      
      // Rotate axis back for the other surface point calculation
      TmpAxis.rotateZ(-PartAngle);

      // Create outer side
      // Calc points for the planar surface boundary
      // The order of the point give the sense 	  
      G4Point3DVector PointList2(4);
      PointList2[0] = LocalOrigin + (RMAX[a] * TmpAxis);
      PointList2[3] = LocalOrigin + (Length*Axis) + (RMAX[a+1] * TmpAxis);
      TmpAxis.rotateZ(PartAngle);
      PointList2[2] = LocalOrigin + (Length*Axis) + (RMAX[a+1] * TmpAxis);
      PointList2[1] = LocalOrigin + (RMAX[a] * TmpAxis);	  
           
      // Add to surface list and set sense	   
      SurfaceVec[Count] = new G4FPlane(&PointList2);

      Count++;
    }
    
    LocalOrigin = LocalOrigin + (Length*Axis);
  }
  
  // Create end planes

  if(dphi >= 2*pi-perMillion)
  {
    // Create only end planes
    G4Point3DVector EndPointList(sides);
    G4Point3DVector InnerPointList(sides);
    G4Point3DVector EndPointList2(sides);
    G4Point3DVector InnerPointList2(sides);	
    TmpAxis = XAxis;
    TmpAxis.rotateZ(phi1);	
    TmpAxis.rotateZ(dphi);
    
    for(int c=0;c<sides;c++)
    {
      // outer polyline for origin end
      EndPointList[c] = Origin + (RMAX[0] * TmpAxis);
      InnerPointList[c] = Origin + (RMIN[0] * TmpAxis);
      EndPointList2[c] = LocalOrigin + (RMAX[sections] * TmpAxis);
      InnerPointList2[c] = LocalOrigin + (RMIN[sections] * TmpAxis);
      TmpAxis.rotateZ(-PartAngle);
    }
    
    // Add to surface list and set sense
    SurfaceVec[nb_of_surfaces-2] =
      new G4FPlane(&EndPointList, &InnerPointList);
    
    // Add to surface list and reverse sense
    SurfaceVec[nb_of_surfaces-1] =
      new G4FPlane(&EndPointList2, &InnerPointList2, 0);
  }
  else
  {
    // If phi section, create a single boundary (case with RMIN=0 included)
    TmpAxis = XAxis;
    TmpAxis.rotateZ(phi1);	
    TmpAxis.rotateZ(dphi);
    
    // Create end planes 
    G4Point3DVector EndPointList ((sides+1)*2);
    G4Point3DVector EndPointList2((sides+1)*2);	      
    
    for(int c=0;c<sides+1;c++)
    {
      // outer polylines for origin end and opposite side
      EndPointList[c]  = Origin + (RMAX[0] * TmpAxis);
      EndPointList[(sides+1)*2-1-c]  = Origin + (RMIN[0] * TmpAxis);
      EndPointList2[c] = LocalOrigin + (RMAX[sections] * TmpAxis);
      EndPointList2[(sides+1)*2-1-c] = LocalOrigin + (RMIN[sections] * TmpAxis);
      TmpAxis.rotateZ(-PartAngle);
    }
    
    // Create the lateral planars
    TmpAxis = XAxis;
    G4ThreeVector TmpAxis2 = XAxis;
    TmpAxis.rotateZ(phi1);
    TmpAxis2.rotateZ(phi1);
    TmpAxis2.rotateZ(dphi);	
    
    LocalOrigin=Origin;
    int points = sections*2+2;
    G4Point3DVector GapPointList(points);
    G4Point3DVector GapPointList2(points);
    Count=0;
    
    for(int d=0;d<sections+1;d++)
    {
      GapPointList[Count] = LocalOrigin + (RMAX[d]*TmpAxis);
      GapPointList[points-1-Count] = LocalOrigin + (RMIN[d]*TmpAxis);	    
      
      GapPointList2[Count] = LocalOrigin + (RMAX[d]*TmpAxis2);
      GapPointList2[points-1-Count] = LocalOrigin + (RMIN[d]*TmpAxis2);	 
   	         
      Count++;

      Length = z_values[d+1] - z_values[d];
      LocalOrigin = LocalOrigin+(Length*Axis);
    }
    
    // Add the lateral planars to the surfaces list and set/reverse sense
    
    SurfaceVec[nb_of_surfaces-4] = new G4FPlane(&GapPointList);
    SurfaceVec[nb_of_surfaces-3] = new G4FPlane(&GapPointList2, 0, 0);
    
    //Add the end planes to the surfaces list and set/reverse sense
    
    if(RMAX[0]-RMIN[0] >= perMillion){
      SurfaceVec[nb_of_surfaces-2] = new G4FPlane(&EndPointList);
    }
    else{
      nb_of_surfaces -= 1;
    };
    
    if(RMAX[sections]-RMIN[sections] >= perMillion){
      SurfaceVec[nb_of_surfaces-1] = new G4FPlane(&EndPointList2, 0, 0);
    }
    else{
      nb_of_surfaces -= 1;
    };    	
    
  }
  
  // Store the original parameters, to be used in visualisation
  // Note radii are scaled because this BREP uses the radius of the
  // inscribed circle but graphics_reps/G4Polyhedron uses the radius of
  // the circumscribed circle.
  original_parameters.Start_angle= phi1;
  original_parameters.Opening_angle= dphi;		   
  original_parameters.Sides= sides;
  original_parameters.Num_z_planes= num_z_planes; 
  original_parameters.Z_values= new G4double[num_z_planes];
  original_parameters.Rmin= new G4double[num_z_planes];
  original_parameters.Rmax= new G4double[num_z_planes];
  G4double rFactor = cos(dphi/(2*sides));

  for(int is=0;is<num_z_planes;is++)
  {
    original_parameters.Z_values[is]= z_values[is]; 
    original_parameters.Rmin[is]= RMIN[is]/rFactor;
    original_parameters.Rmax[is]= RMAX[is]/rFactor;
  }

  ///////////////////////////////////////////////////
  // Print for debugging

#ifdef G4VERBOSE
  static G4int print_pgone_parameters = 1;
  
  if(print_pgone_parameters)
  {
    G4cout << "Parameters of the G4PGone " << name << endl;
    G4cout << "  starting angle =" << original_parameters.Start_angle << endl;
    G4cout << "  opening angle =" << original_parameters.Opening_angle << endl;
    G4cout << "  sides =" << original_parameters.Sides << endl;
    G4cout << "  nb of z planes=" << original_parameters.Num_z_planes << endl;

    for (G4int nb = 0; nb <= sections; nb++)
    {
      G4cout << "   Z[" << nb << "] = " << original_parameters.Z_values[nb];
      G4cout << "   Rmin[" << nb << "] = " << original_parameters.Rmin[nb];
      G4cout << "   Rmax[" << nb << "] = " << original_parameters.Rmax[nb] 
	     << endl;
    }   
  }
#endif


 
  // z_values[0]  should be equal to z_start, for consistency 
  //   with what the constructor does.
  // Otherwise the z_values that are shifted by (z_values[0] - z_start) , 
  //   because z_values are only used in the form 
  //      length = z_values[d+1] - z_values[d];         // JA Apr 2, 97
  
  if( z_values[0] != z_start )
  {
    G4cerr << "ERROR in creating G4BREPSolidPolyhedra: "  << 
      " z_values[0]= " << z_values[0] << " is not equal to " << 
      " z_start= " , z_start;
    // G4Exception(" Error in creating G4BREPSolidPolyhedra: z_values[0] must be equal to z_start" );
    original_parameters.Z_values[0]= z_start; 
  }

  active=1;
  Initialize(); 
}


G4BREPSolidPolyhedra::~G4BREPSolidPolyhedra()
{
  delete[] original_parameters.Z_values;
  delete[] original_parameters.Rmin;
  delete[] original_parameters.Rmax;
}


void G4BREPSolidPolyhedra::Initialize()
{
  // Calc bounding box for solids and surfaces
  // Convert concave planes to convex     
  ShortestDistance=1000000;
  CheckSurfaceNormals();
  if(!Box || !AxisBox)
    IsConvex();
  
  CalcBBoxes();
}


EInside G4BREPSolidPolyhedra::Inside(register const G4ThreeVector& Pt) const
{
  // This function find if the point Pt is inside, 
  // outside or on the surface of the solid


  G4double halfTolerance = kCarTolerance*0.5;

  G4Vector3D v(1, 0, 0.01);
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(v);
  G4Ray r(Pttmp, Vtmp);
  
  // Check if point is inside the Polyhedra bounding box
  if( !GetBBox()->Inside(Pttmp) )
    return kOutside;

  // Set the surfaces to active again
  Reset();
  
  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  G4int hits=0, samehit=0;

  for(G4int a=0; a < nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->Active())
    {
      // count the number of intersections.
      // if this number is odd, the start of the ray is
      // inside the volume bounded by the surfaces, so 
      // increment the number of intersection by 1 if the 
      // point is not on the surface and if this intersection 
      // was not founded before
      if( (SurfaceVec[a]->Intersect(r)) & 1 )
      {
	// test if the point is on the surface
	if(SurfaceVec[a]->Distance() <= kCarTolerance*kCarTolerance)
	  return kSurface;
	
	// test if this intersection was founded before
	for(G4int i=0; i<a; i++)
	  if(SurfaceVec[a]->Distance() == SurfaceVec[i]->Distance())
	  {
	    samehit++;
	    break;
	  }
	
	// count the number of surfaces intersected by the ray
	if(!samehit)
	  hits++;
      }
    }
  }
   
  // if the number of surfaces intersected is odd,
  // the point is inside the solid
  if(hits&1)
    return kInside;
  else
    return kOutside;
}


G4ThreeVector G4BREPSolidPolyhedra::SurfaceNormal
                                         (const G4ThreeVector& Pt) const
{
  // This function calculates the normal of the surface
  // at a point on the surface
  // Note : the sense of the normal depends on the sense of the surface 

  G4Vector3D   n(0,0,0);
  G4int        iplane;
    
  G4Vector3D norm;
  G4Ray r( Pt, G4Vector3D(1, 0, 0) );

  // Find on which surface the point is
  for(iplane = 0; iplane < nb_of_surfaces; iplane++)
  {
    if(SurfaceVec[iplane]->HowNear(Pt) < kCarTolerance)
      // the point is on this surface
      break;
  }
  
  // calcul of the normal at this point
  norm = SurfaceVec[iplane]->SurfaceNormal(Pt);

  n = G4ThreeVector ( norm.x(), norm.y(), norm.z() );
  n = n.unit();

  return n;
}


G4double G4BREPSolidPolyhedra::DistanceToIn(const G4ThreeVector& Pt) const
{
  // Calculates the shortest distance ("safety") from a point
  // outside the solid to any boundary of this solid.
  // Return 0 if the point is already inside.


  G4double *dists = new G4double[nb_of_surfaces];
  G4double halfTolerance = kCarTolerance*0.5;  
  G4int a;

  // Set the surfaces to active again
  Reset();
  
  // calcul of the shortest distance of the point to each surfaces
  // Be carreful : it's a signed value
  for(a=0; a< nb_of_surfaces; a++)  
    dists[a] = SurfaceVec[a]->HowNear(Pt);
     
  G4double Dist = kInfinity;
  
  // if dists[] is positive, the point is outside
  // so take the shortest of the shortest positive distances
  // dists[] can be equal to 0 : point on a surface
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is inside to return 0, utilize the Inside
  //   function. But I don`t know if it is really needed because dToIn is 
  //   called only if the point is outside )
  for(a = 0; a < nb_of_surfaces; a++)
    if( fabs(Dist) > fabs(dists[a]) ) 
      //if( dists[a] >= 0)
	Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    // the point is inside the solid or on a surface
    return 0;
  else 
    //return Dist;
    return fabs(Dist);
}


G4double G4BREPSolidPolyhedra::DistanceToIn(register const G4ThreeVector& Pt, 
					register const G4ThreeVector& V) const
{
  // Calculates the distance from a point outside the solid
  // to the solid`s boundary along a specified direction vector.
  // 	
  // Note : Intersections with boundaries less than the 
  //	    tolerance must be ignored if the direction 
  // 	    is away from the boundary
  
  G4int a;
  
  // Set the surfaces to active again
  Reset();
  
  G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
  
  for(a=0; a< nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->Active())
    {
      // test if the ray intersect the surface
      if( (SurfaceVec[a]->Intersect(r)) )
      {
	// if more than 1 surface is intersected,
	// take the nearest one
	if( SurfaceVec[a]->Distance() < ShortestDistance )
	  if( SurfaceVec[a]->Distance() > halfTolerance )
	  {
	    ShortestDistance = SurfaceVec[a]->Distance();
	  }
	  else
	  {
	    // the point is within the boundary
	    // ignored it if the direction is away from the boundary	     
	    G4Vector3D Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);

	    if( (Norm * Vtmp) < 0 )
	      ShortestDistance = SurfaceVec[a]->Distance();
	  }
      }
    }
  }
  
  // Be carreful !
  // SurfaceVec->Distance is in fact the squared distance
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  else
    // no intersection, return kInfinity
    return kInfinity; 
}


G4double G4BREPSolidPolyhedra::DistanceToOut(register const G4ThreeVector& Pt, 
					 register const G4ThreeVector& V, 
					 const G4bool calcNorm, 
					 G4bool *validNorm, 
					 G4ThreeVector *n            ) const
{
  // Calculates the distance from a point inside the solid
  // to the solid`s boundary along a specified direction vector.
  // Return 0 if the point is already outside.
  //
  // Note : If the shortest distance to a boundary is less 
  // 	    than the tolerance, it is ignored. This allows
  // 	    for a point within a tolerant boundary to leave
  //	    immediately

  // Set the surfaces to active again
  Reset();

  const G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Ptv = Pt;
  G4int a;

  // I don`t understand this line
  if(validNorm)
    *validNorm=false;

  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
 
  for(a=0; a< nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->Active())
    {
      // test if the ray intersect the surface
      if( (SurfaceVec[a]->Intersect(r)) )
      {
	// if more than 1 surface is intersected,
	// take the nearest one
	if( SurfaceVec[a]->Distance() < ShortestDistance )
	  if( SurfaceVec[a]->Distance() > halfTolerance )
	  {
	    ShortestDistance = SurfaceVec[a]->Distance();
	  }
	  else
	  {
	    // the point is within the boundary: ignored it
	  }
      }
    }
  }

  // Be carreful !
  // SurfaceVec->Distance is in fact the squared distance
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  else
    // if no intersection is founded, the point is outside
    // so return 0
    return 0; 
}


G4double G4BREPSolidPolyhedra::DistanceToOut(const G4ThreeVector& Pt) const
{
  // Calculates the shortest distance ("safety") from a point
  // inside the solid to any boundary of this solid.
  // Return 0 if the point is already outside.	

  G4double *dists = new G4double[nb_of_surfaces];
  G4double halfTolerance = kCarTolerance*0.5;  
  G4int a;

  // Set the surfaces to active again
  Reset();
  
  // calcul of the shortest distance of the point to each surfaces
  // Be carreful : it's a signed value
  for(a=0; a< nb_of_surfaces; a++)
    dists[a] = SurfaceVec[a]->HowNear(Pt);  

  G4double Dist = kInfinity;
  
  // if dists[] is negative, the point is inside
  // so take the shortest of the shortest negative distances
  // dists[] can be equal to 0 : point on a surface
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is outside to return 0, utilize the Inside
  //   function. But I don`t know if it is really needed because dToOut is 
  //   called only if the point is inside )

  for(a = 0; a < nb_of_surfaces; a++)
    if( fabs(Dist) > fabs(dists[a]) ) 
      //if( dists[a] <= 0)
	Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    // the point is ouside the solid or on a surface
    return 0;
  else
    // return Dist;
    return fabs(Dist);
}


//  In graphics_reps:   
#include "G4Polyhedron.hh"   

G4Polyhedron* G4BREPSolidPolyhedra::CreatePolyhedron() const
{
  return new G4PolyhedronPgon( original_parameters.Start_angle, 
			       original_parameters.Opening_angle, 
			       original_parameters.Sides, 
			       original_parameters.Num_z_planes, 
			       original_parameters.Z_values,
			       original_parameters.Rmin,
			       original_parameters.Rmax);
}
