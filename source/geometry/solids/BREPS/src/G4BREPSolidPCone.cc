// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidPCone.cc,v 1.1 1999-01-07 16:07:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#include "G4BREPSolidPCone.hh"
#include "G4FCylindricalSurface.hh"
#include "G4FConicalSurface.hh"
#include "G4CircularCurve.hh"
#include "G4FPlane.hh"


G4BREPSolidPCone::G4BREPSolidPCone(G4String name,
				   const G4double start_angle,
				   const G4double opening_angle,
				   const int      num_z_planes, // sections,
				   const G4double z_start,		   
				   const G4double z_values[],
				   const G4double RMIN[],
				   const G4double RMAX[]
				   ): G4BREPSolid(name)

{
  const int  sections= num_z_planes-1;
  nb_of_surfaces = 2*sections+2;
  SurfaceVec = new G4Surface*[nb_of_surfaces];
  G4ThreeVector Axis(0,0,1);
  G4ThreeVector Origin(0,0,z_start);    
  G4double Length;
  G4ThreeVector LocalOrigin(0,0,z_start);    
  G4int a, b = 0;

  G4ThreeVector PlaneAxis(0, 0, 1);
  G4ThreeVector PlaneDir (0, 1, 0);   

  ///////////////////////////////////////////////////
  // Temporary
  for (G4int x = 0; x <= sections; x++)
  {
    G4cout<<"Z"<<x<<"="<<z_values[x];
    G4cout<<" Rmin"<<x<<"="<<RMIN[x];
    G4cout<<" Rmax"<<x<<"="<<RMAX[x]<<endl;
  }

  G4cout<<"start angle ="<<start_angle<<endl;
  G4cout<<"open angle ="<<opening_angle<<endl;
  G4cout<<"zstart ="<<z_start<<endl;

  ///////////////////////////////////////////////////
  // Test the validity of the R values
  
  // RMIN[0] and RMIN[num_z_planes-1] cannot be = 0
  // when RMIN[0] or RMIN[num_z_planes-1] are = 0
  if(  ((RMIN[0] == 0)              && (RMAX[0] == 0))             || 
       ((RMIN[num_z_planes-1] == 0) && (RMAX[num_z_planes-1] == 0))   )
    G4Exception("RMIN at the extremities can not be nul when RMAX = 0");  
  
  // only RMAX[0] and RMAX[num_z_planes-1] can be = 0
  for(a = 1; a < num_z_planes-1; a++)
    if (RMAX[a] == 0)
      G4Exception("RMAX inside the solid  can not be nul");
  
  // RMAX[a] must be greater than RMIN[a] 
  for(a = 0; a < num_z_planes; a++)
    if (RMIN[a] >= RMAX[a])
      G4Exception("RMAX must be greater than RMIN");
  

  ///////////////////////////////////////////////////
  // Create cylindrical et conical surfaces
  
  for(a=0; a<sections; a++)
  {
    // Surface length
    Length = z_values[a+1] - z_values[a];
    
    if (Length == 0)
    {
      // The surface to create is planar

      G4double R1, R2;
      // test where is the plane surface
      if(RMAX[a] != RMAX[a+1])
      {
	R1 = RMAX[a];
	R2 = RMAX[a+1];
      }
      else if(RMIN[a] != RMIN[a+1])
      {
	R1 = RMIN[a];
	R2 = RMIN[a+1];
      }
      else
      {
        G4cerr << "Error in construction of G4BREPSolidPCone: "
               <<  "Exactly the same z, rmin and rmax given for "
               <<  "consecutive indices, " << a << " and " << a+1 << endl;
        // G4Exception("G4BREPSolidPCone constructor: Error in parameter values");
        continue; 
      }

      // Create plane surface
      G4Point3D ArcStart1 = LocalOrigin + (R1*PlaneDir);
      G4Point3D ArcStart2 = LocalOrigin + (R2*PlaneDir);
      
      G4CurveVector cv1;
      G4CircularCurve *tmp1, *tmp2;
      
      if(R1 != 0)  
      {
	tmp1 = new G4CircularCurve;
	tmp1->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), R1);
	tmp1->SetBounds(ArcStart1, ArcStart1);
	if(R1>R2)
	  tmp1->SetSameSense(1);
	else
	  tmp1->SetSameSense(0);

	cv1.append(tmp1);
      }

      if(R2 != 0)  
      {
	tmp2 = new G4CircularCurve;
	tmp2->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), R2);
	tmp2->SetBounds(ArcStart2, ArcStart2);
	if(R1>R2)
	  tmp2->SetSameSense(0);
	else
	  tmp2->SetSameSense(1);
	cv1.append(tmp2);
      }
	
      SurfaceVec[b]   = new G4FPlane(PlaneDir, PlaneAxis, LocalOrigin);
      SurfaceVec[b]->SetBoundaries(&cv1);

      nb_of_surfaces--;
      b++;

    }
    else
    {
      // The surface to create is conical or cylindrical

      // Inner PCone
      if(RMIN[a] != RMIN[a+1])
      {
	// Create cone
	if(RMIN[a] > RMIN[a+1])
	{
	  G4Vector3D ConeOrigin = LocalOrigin ; 
	  
	  SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis, Length, 
						RMIN[a+1], RMIN[a]);
	}
	else
	{      
	  G4Vector3D Axis2 = (-1*Axis);
	  G4Vector3D LocalOrigin2 = LocalOrigin + (Length*Axis);
	  G4Vector3D ConeOrigin = LocalOrigin2 ; 
	  
	  SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis2, 
						Length, RMIN[a], RMIN[a+1]); 
	}
      
	b++;  
      }
      else
      {
	if (RMIN[a] == 0)
	{
	  // Do not create any surface
	  // and decrease nb_of_surfaces
	  nb_of_surfaces--;
	}
	else
	{
	  // Create cylinder
	  G4Vector3D CylOrigin = LocalOrigin ; 
	  
	  SurfaceVec[b] = new G4FCylindricalSurface(CylOrigin, Axis, 
						    RMIN[a], Length );
	  b++;
	}	    
      }
    
      // Outer PCone
      if(RMAX[a] != RMAX[a+1])
      {
	// Create cone
	if(RMAX[a] > RMAX[a+1])
	{
	  G4Vector3D ConeOrigin = LocalOrigin ;
	  
	  SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis, 
						Length, RMAX[a+1], RMAX[a]);
	}
	else
	{
	  G4Vector3D Axis2 = (-1*Axis);
	  G4Vector3D LocalOrigin2 = LocalOrigin + (Length*Axis);
	  G4Vector3D ConeOrigin = LocalOrigin2 ;

	  SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis2, 
						Length, RMAX[a], RMAX[a+1]);
	}
	
	b++;
      }
      else
      {
	// Create cylinder
	G4Vector3D CylOrigin = LocalOrigin ;
	
	if (RMAX[a] == 0)
	{
	  // Do not create any surface
	  // and decrease nb_of_surfaces
	  nb_of_surfaces--;
	}
	else
	{
	  // Create cylinder
	  G4Vector3D CylOrigin = LocalOrigin ; 
	  
	  SurfaceVec[b] = new G4FCylindricalSurface(CylOrigin, Axis, 
						    RMAX[a], Length );
	  b++;
	}	  	
      }
    }

    // Move surface origin to next section
    LocalOrigin = LocalOrigin + (Length*Axis);
  }
  
  ///////////////////////////////////////////////////
  // Create two end planes  
   
  // Create start G4Plane & boundaries    
  G4Point3D ArcStart1a = Origin + (RMIN[0]*PlaneDir);
  G4Point3D ArcStart1b = Origin + (RMAX[0]*PlaneDir);
 
  G4CurveVector cv;
  G4CircularCurve* tmp;

  tmp = new G4CircularCurve;
  tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, Origin), RMIN[0]);
  tmp->SetBounds(ArcStart1a, ArcStart1a);
  tmp->SetSameSense(0);
  cv.append(tmp);
  
  tmp = new G4CircularCurve;
  tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, Origin), RMAX[0]);
  tmp->SetBounds(ArcStart1b, ArcStart1b);
  tmp->SetSameSense(1);
  cv.append(tmp);

  SurfaceVec[nb_of_surfaces-2]   = new G4FPlane(PlaneDir, PlaneAxis, Origin);
  SurfaceVec[nb_of_surfaces-2]->SetBoundaries(&cv);
 

  // Create end G4Plane & boundaries    
  G4Point3D ArcStart2a = LocalOrigin + (RMIN[sections]*PlaneDir);  
  G4Point3D ArcStart2b = LocalOrigin + (RMAX[sections]*PlaneDir);    
  
  cv.clear();

  tmp = new G4CircularCurve;
  tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), 
	    RMIN[sections]);
  tmp->SetBounds(ArcStart2a, ArcStart2a);
  tmp->SetSameSense(0);
  cv.append(tmp);

  tmp = new G4CircularCurve;
  tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), 
	    RMAX[sections]);
  tmp->SetBounds(ArcStart2b, ArcStart2b);
  tmp->SetSameSense(1);
  cv.append(tmp);
  
  SurfaceVec[nb_of_surfaces-1]= new G4FPlane(PlaneDir, PlaneAxis, LocalOrigin);
  SurfaceVec[nb_of_surfaces-1]->SetBoundaries(&cv);
  

  active=1;
  Initialize();

  // Store the original parameters, to be used in visualisation
  original_parameters.Start_angle= start_angle;
  original_parameters.Opening_angle= opening_angle;		   
  original_parameters.Num_z_planes= num_z_planes; 
  // original_parameters.z_start= z_start;		   
  original_parameters.Z_values= new G4double[num_z_planes];
  original_parameters.Rmin= new G4double[nb_of_surfaces];
  original_parameters.Rmax= new G4double[nb_of_surfaces];
  
  for(int is=0;is<num_z_planes;is++)
  {
    original_parameters.Z_values[is]= z_values[is]; 
    original_parameters.Rmin[is]= RMIN[is];
    original_parameters.Rmax[is]= RMAX[is];
  }

  // z_values[0]  should be equal to z_start, for consistency 
  //   with what the constructor does.
  // Otherwise the z_values that are given are used 
  //   shifted by   z_values[0] - z_start: 
  //  (because z_values are only used in 
  //    line 26:      Length = z_values[a+1] - z_values[a]; 
  //  )                                                      // JA Apr 2, 97
  
  /*
  if( z_values[0] != z_start )
  {
  G4cerr << "ERROR in creating G4BREPSolidPCone: "  
  << " z_values[0]= " << z_values[0] << " is not equal to " 
  << " z_start= " , z_start; 
  // G4Exception(" Error in creating G4BREPSolidPCone: z_values[0] must be equal to z_start" );
  original_parameters.Z_values[0]= z_start;

  }
  */ 
  
}


G4BREPSolidPCone::~G4BREPSolidPCone()
{
  delete[] original_parameters.Z_values;
  delete[] original_parameters.Rmin;
  delete[] original_parameters.Rmax;
}


void G4BREPSolidPCone::Initialize()
{ 
  // Calc bounding box for solids and surfaces
  // Convert concave planes to convex     
  ShortestDistance=1000000;
  CheckSurfaceNormals();
  
  if(!Box || !AxisBox)
    IsConvex();
  
  CalcBBoxes();
}


EInside G4BREPSolidPCone::Inside(register const G4ThreeVector& Pt) const
{
  // Check if point lies between end planes of PCone
  G4double dist1 = SurfaceVec[nb_of_surfaces-1]->ClosestDistanceToPoint(Pt);
  G4double dist2 = SurfaceVec[nb_of_surfaces-2]->ClosestDistanceToPoint(Pt);

  if((dist1 < -kCarTolerance && dist2 <-kCarTolerance)||
     (dist1 > kCarTolerance && dist2 >kCarTolerance)    )
    return kOutside;

  G4Vector3D v(1,0,0);
  G4double Dist;
  G4double halfTolerance = kCarTolerance*0.5;
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(v);
  G4Ray r(Pttmp, Vtmp);
  TestSurfaceBBoxes(r);
  G4int hits=0;

  for(G4int a=0; a < nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->Active())
      if(SurfaceVec[a]->Intersect(r))
      {
	if(SurfaceVec[a]->Distance() < kCarTolerance)
	  return kSurface;
	
	hits++;
      }
  }

  // Set the surfaces to active again
  for(G4int b=0; b < nb_of_surfaces; b++)
    SurfaceVec[b]->Reset();


  if(hits&1)
    return kInside;
  
  return kOutside;
}


G4ThreeVector G4BREPSolidPCone::SurfaceNormal(const G4ThreeVector& Pt) const
{
  G4cout<<" SurfaceNormal() of G4BREPSolidPCone modified by L. Broglia";
  
  G4Vector3D   Ptv          = Pt;
  G4Vector3D   n(0,0,0);
  G4double     zCoord       = Pt.z();
  const G4int  num_z_planes = original_parameters.Num_z_planes;
  G4int        iplane;
 
  // Find the appropriate z "slice" 
  //
  for(iplane=0; iplane< num_z_planes; iplane++) 
    if (  (zCoord < original_parameters.Z_values[iplane+1]) &&
	  (zCoord >= original_parameters.Z_values[iplane])      ) 
      break;
    
  G4Vector3D norm;
  G4Ray r( Pt, G4Vector3D(1, 0, 0) );

  // We must find which is the correct surface, the inner or the outer one
  // (if they exist)
  for(iplane = 0; iplane < num_z_planes; iplane++)
  {
    // check if the point is on the surface
    if(SurfaceVec[iplane]->Intersect(r))
      if(SurfaceVec[iplane]->Distance() < kCarTolerance)
	// the point is on the surface
	break;
  }
  
  norm =  SurfaceVec[iplane]->SurfaceNormal(Pt);

  n = G4ThreeVector ( norm.x(), norm.y(), norm.z());
  n = n.unit();

/*
  if ( SurfaceVec[innerSurface]->WithinBoundary(Ptv) == 1 ) 
  {
    norm =  SurfaceVec[ innersurface ]->SurfaceNormal(Pt);
  }
  else  if ( SurfaceVec[outerSurface]->WithinBoundary(Ptv) == 1 ) 
  {
    norm =  SurfaceVec[ outerSurface]->SurfaceNormal(Pt);
  } 

  // Check if it is on one of the top/bottom planes
  //
  if ( fabs(zCoord - original_parameters.Z_values[0]) < kCarTolerance ) 
  {
    //    n = G4ThreeVector (0., 0., sign( original_parameters.Z_values[0]
    // 				    -original_parameters.Z_values[1]) );
    
    n = G4ThreeVector (0., 0., original_parameters.Z_values[0]
		       -original_parameters.Z_values[1] );
    n = n.unit();
  }
  else if (fabs(zCoord - original_parameters.Z_values[num_z_planes-1]) < 
	   kCarTolerance) 
    {
    n = G4ThreeVector(0., 0., original_parameters.Z_values[num_z_planes]
		      -original_parameters.Z_values[num_z_planes-1] ); 
    n = n.unit();
  }
*/

  return n;
}


G4double G4BREPSolidPCone::DistanceToIn(const G4ThreeVector& Pt) const
{
  G4double *dists = new G4double[nb_of_surfaces];
  G4double halfTolerance = kCarTolerance*0.5;  
  G4int a;

  for(a=0; a< nb_of_surfaces;a++)
    dists[a] = fabs(SurfaceVec[a]->HowNear(Pt));
  
  G4double Dist=kInfinity;
  
  for(a=0; a< nb_of_surfaces;a++)
    if(Dist>dists[a]) Dist = dists[a];
  
  delete[] dists;

  // Set the surfaces to active again
  for(G4int b=0; b < nb_of_surfaces; b++)
    SurfaceVec[b]->Reset();
 
  return Dist;
}


G4double G4BREPSolidPCone::DistanceToIn(register const G4ThreeVector& Pt, 
					register const G4ThreeVector& V) const
{
  int a;
  Reset();
  G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  //  G4double kInfinity = ;
  G4Ray r(Pttmp, Vtmp);
  TestSurfaceBBoxes(r);
  ShortestDistance = kInfinity;
  
  for(a=0; a< nb_of_surfaces;a++)
  {
    if(SurfaceVec[a]->Active())
      if(SurfaceVec[a]->Intersect( r ))
      {
	if(ShortestDistance > SurfaceVec[a]->Distance())
	  if(SurfaceVec[a]->Distance()> halfTolerance)
	  {
	    ShortestDistance = SurfaceVec[a]->Distance();
	  }
	  else
	  {
	    G4Vector3D Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);
	    if((Norm * Vtmp)<0)
	      ShortestDistance = SurfaceVec[a]->Distance();
	  }
      }
  }

  // Set the surfaces to active again
  for(G4int b=0; b < nb_of_surfaces; b++)
    SurfaceVec[b]->Reset();

  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  
  return kInfinity; 
}


G4double G4BREPSolidPCone::DistanceToOut(register const G4ThreeVector& Pt, 
					 register const G4ThreeVector& V, 
					 const G4bool calcNorm, 
					 G4bool *validNorm, 
					 G4ThreeVector *n            ) const
{
  const G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Ptv = Pt;
  G4double wb = 0.0;
  G4int a;

  for( a=0; a< nb_of_surfaces-2; a++) 
  {
    wb = fabs( SurfaceVec[a]->HowNear(Ptv) );
  
    //  If we are on a surface and exiting it return Zero
    if ( (wb < halfTolerance) && (V.dot(SurfaceVec[a]->Normal(Ptv))>0) ) 
    {
      return (0.0);
    }
  }
  
  if(validNorm)
    *validNorm=false;
  
  Reset();  

  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  
  // G4double kInfinity = 10e20;
  G4Ray r(Pttmp, Vtmp);
  TestSurfaceBBoxes(r);
  ShortestDistance = kInfinity;
  
  for(a=0; a< nb_of_surfaces;a++)
    if(SurfaceVec[a]->Active()) 
      if(SurfaceVec[a]->Intersect( r ))
	if(ShortestDistance > SurfaceVec[a]->Distance()&&
	   SurfaceVec[a]->Distance()> halfTolerance) 	
	  ShortestDistance = SurfaceVec[a]->Distance();
  
  // Set the surfaces to active again
  for(G4int b=0; b < nb_of_surfaces; b++)
    SurfaceVec[b]->Reset();

  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  
  return kInfinity; 
}


G4double G4BREPSolidPCone::DistanceToOut(const G4ThreeVector& Pt) const
{
  int a;

  G4double *dists = new G4double[nb_of_surfaces];
  G4double halfTolerance = kCarTolerance*0.5;    

  for(a=0; a< nb_of_surfaces; a++)
    dists[a] = fabs(SurfaceVec[a]->HowNear(Pt));
  
  G4double Dist=kInfinity;
  
  for(a=0; a< nb_of_surfaces;a++)
     if( Dist>dists[a] ) Dist = dists[a];

  // Set the surfaces to active again
  for(G4int b=0; b < nb_of_surfaces; b++)
    SurfaceVec[b]->Reset();


  // If we are on a surface, the return value Dist must be zero!
  delete[] dists;	
  return Dist;
}


//  In graphics_reps:   

#include "G4Polyhedron.hh"   

G4Polyhedron* G4BREPSolidPCone::CreatePolyhedron() const
{
  return new G4PolyhedronPcon( original_parameters.Start_angle, 
			       original_parameters.Opening_angle, 
			       original_parameters.Num_z_planes, 
			       original_parameters.Z_values,
			       original_parameters.Rmin,
			       original_parameters.Rmax);
}














