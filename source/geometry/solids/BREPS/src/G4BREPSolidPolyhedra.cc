// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BREPSolidPolyhedra.cc,v 1.1 1999-01-07 16:07:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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
  
  if(dphi == 2*pi)
    nb_of_surfaces = 2*(sections * sides) + 2;
  else
    nb_of_surfaces = 2*(sections * sides) + 4;


  SurfaceVec = new G4Surface*[nb_of_surfaces];
  
  G4Vector3D Axis(0,0,1);
  G4Vector3D YAxis(0,1,0);
  G4Vector3D TmpAxis;
  G4Point3D Origin(0,0,z_start);    
  G4Point3D LocalOrigin(0,0,z_start);    
  G4double Length;
  int Count =0;
  G4double  PartAngle = (dphi - phi1)/sides;

  ///////////////////////////////////////////////////
  // Temporary

  for (G4int x = 0; x <= sections; x++)
  {
    cout<<"Z"<<x<<"="<<z_values[x];
    cout<<" Rmin"<<x<<"="<<RMIN[x];
    cout<<" Rmax"<<x<<"="<<RMAX[x]<<endl;
  }

  cout<<"phi1 ="<<phi1<<endl;
  cout<<"dphi ="<<dphi<<endl;
  cout<<"sides ="<<sides<<endl;
  cout<<"zstart ="<<z_start<<endl;


  ///////////////////////////////////////////////////
   

  for(G4int a=0;a<sections;a++)
  {
    TmpAxis= YAxis;
    TmpAxis.rotateZ(phi1);
    Length = z_values[a+1] - z_values[a];
    
    // Create sides
    for(int b=0;b<sides;b++)
    {
      G4Point3DVector PointList(4);
      // Create inner side
      // Calc points for the planar surface boundary
      PointList[0] = LocalOrigin + (RMIN[a] * TmpAxis);
      PointList[1] = LocalOrigin + (Length*Axis) + (RMIN[a+1] * TmpAxis);
      TmpAxis.rotateZ(PartAngle);
      PointList[2] = LocalOrigin + (Length*Axis) + (RMIN[a+1] * TmpAxis);
      PointList[3] = LocalOrigin + (RMIN[a] * TmpAxis);	  
      SurfaceVec[Count] = new G4FPlane( &PointList);
      Count++;
      
      // Rotate axis back for the other surface point calculation
      TmpAxis.rotateZ(-PartAngle);
      
      // Calc points for the planar surface boundary	  
      G4Point3DVector PointList2(4);
      PointList2[0] = LocalOrigin + (RMAX[a] * TmpAxis);
      PointList2[1] = LocalOrigin + (Length*Axis) + (RMAX[a+1] * TmpAxis);
      TmpAxis.rotateZ(PartAngle);
      PointList2[2] = LocalOrigin + (Length*Axis) + (RMAX[a+1] * TmpAxis);
      PointList2[3] = LocalOrigin + (RMAX[a] * TmpAxis);	  
      SurfaceVec[Count] = new G4FPlane(&PointList2);
      Count++;
    }
    
    LocalOrigin = LocalOrigin + (Length*Axis);
  }
  
  // Create end planes

  if(dphi == 2*pi)
  {
    // Create only end planes
    G4Point3DVector EndPointList(sides);
    G4Point3DVector InnerPointList(sides);
    G4Point3DVector EndPointList2(sides);
    G4Point3DVector InnerPointList2(sides);	
    TmpAxis = YAxis;
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
    
    SurfaceVec[nb_of_surfaces-2] =
      new G4FPlane(&EndPointList, &InnerPointList);
    
    SurfaceVec[nb_of_surfaces-1] =
      new G4FPlane(&EndPointList2, &InnerPointList2);	
  }
  else
  {
    TmpAxis = YAxis;
    TmpAxis.rotateZ(phi1);	
    TmpAxis.rotateZ(dphi);
    
    // Create end planes & two planes for the "missing" part 
    G4Point3DVector   EndPointList(sides+2);
    G4Point3DVector InnerPointList(sides+2);
    G4Point3DVector EndPointList2(sides+2);
    G4Point3DVector InnerPointList2(sides+2);	
    TmpAxis = YAxis;      
    
    for(int c=0;c<sides+1;c++)
    {
      // outer polyline for origin end
      EndPointList[c] = Origin + (RMAX[0] * TmpAxis);
      InnerPointList[c] = Origin + (RMIN[0] * TmpAxis);
      EndPointList2[c] = LocalOrigin + (RMAX[sections] * TmpAxis);
      InnerPointList2[c] = LocalOrigin + (RMIN[sections] * TmpAxis);
      TmpAxis.rotateZ(-PartAngle);
    }

    // Create the extra points on the axis
    TmpAxis = YAxis;
    TmpAxis.rotateZ(phi1);
    EndPointList[sides+1] = Origin;
    InnerPointList[sides+1] = Origin;
    EndPointList2[sides+1] = LocalOrigin;
    InnerPointList2[sides+1] = LocalOrigin;
    int points = sides+2;
    
    SurfaceVec[nb_of_surfaces-4] =
      new G4FPlane(&EndPointList, &InnerPointList);
    
    SurfaceVec[nb_of_surfaces-3] =
      new G4FPlane(&EndPointList2, &InnerPointList2);    
    
    // Create the planars for the "gap"
    TmpAxis = YAxis;
    G4ThreeVector TmpAxis2 = YAxis;
    TmpAxis.rotateZ(phi1);
    TmpAxis2.rotateZ(phi1);
    TmpAxis2.rotateZ(dphi);	
    
    LocalOrigin=Origin;
    points = sections*2+2;
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
    
    SurfaceVec[nb_of_surfaces-2] = new G4FPlane(&GapPointList);
    SurfaceVec[nb_of_surfaces-1] = new G4FPlane(&GapPointList2);	
    
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
  G4double *dists =  new G4double[nb_of_surfaces];
  G4double Dist = kInfinity;
  G4double tmpdist=kInfinity;
  G4double halfTolerance = kCarTolerance*0.5;  
  
  for(int a=0; a< nb_of_surfaces;a++)
  {
    tmpdist = (SurfaceVec[a]->HowNear(Pt));
    
    if(fabs(Dist) > fabs(tmpdist))
      Dist = tmpdist;
  }

  if(Dist > halfTolerance)
    return kOutside;
  
  if(Dist < -halfTolerance)
    return kInside;    
  
  return kSurface;
}


G4ThreeVector G4BREPSolidPolyhedra::SurfaceNormal
                                         (const G4ThreeVector& Pt) const
{
  G4cout<<" SurfaceNormal() of G4BREPSolidPolyhedra modified by L. Broglia";

  /* 
  //G4Exception(" SurfaceNormal() of G4BREPSolidPolyhedra is not yet implemented." );
 
#ifdef WILL_IMPLEMENT
  G4double  zCoord= Pt.z();
  G4int     a,  zSlice, phiSlice;

  // Try to find the appropriate z "slice"
  for(a=0; a< nb_of_surfaces-2;a++)
    if (  (zCoord <  original_parameters.Z_values[a+1]) 
	  &&(zCoord >= original_parameters.Z_values[a]) ) 
      break; 
    
  zSlice= a; 
  
  // Try to find the appropriate phi plane
  phiSlice = 0; 
  
  unsigned int isurface=0;
  
  G4ThreeVec norm =  SurfaceVec[isurface]->SurfaceNormal(Pt);
  G4ThreeVector normalVector = G4ThreeVector ( norm.GetX(), 
					       norm.GetY(), 
					       norm.GetZ());
#else
  G4ThreeVector normalVector = G4ThreeVector ( 0.0, 0.0, 1.0);
#endif
*/

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

  return n;
}


G4double G4BREPSolidPolyhedra::DistanceToIn(const G4ThreeVector& Pt) const
{
  G4double Dist=kInfinity;
  G4double tmpdist=kInfinity;  
  G4double halfTolerance = kCarTolerance*0.5;
  
  for(int a=0; a< nb_of_surfaces;a++)
  {
    tmpdist = fabs(SurfaceVec[a]->HowNear(Pt));
    if(Dist>tmpdist && tmpdist> halfTolerance) Dist = tmpdist;
  }
  
  return Dist;
}


G4double G4BREPSolidPolyhedra::DistanceToIn(register const G4ThreeVector& Pt, 
					    register const G4ThreeVector& V
					    ) const
{
  Reset();  
  G4Point3D Pttmp(Pt);
  G4Vector3D Vtmp(V);
  G4double halfTolerance = kCarTolerance*0.5;  
 
  //  G4double kInfinity = kInfinity;
  G4Ray r(Pttmp, Vtmp);
  TestSurfaceBBoxes(r);
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  ShortestDistance = kInfinity;
  
  for(int a=0; a< nb_of_surfaces;a++)
  {
    if(SurfaceVec[a]->Active())
      // L. Broglia : old
      // if(SurfaceVec[a]->Intersect(r))
      if( (G4FPlane*)(SurfaceVec[a])->Evaluate(r) )
	if(ShortestDistance > SurfaceVec[a]->Distance())
	  if(SurfaceVec[a]->Distance() > halfTolerance)
	    ShortestDistance = SurfaceVec[a]->Distance();
	  else
	  {
	    G4ThreeVector Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);
	    
	    if((Norm * Vtmp)<0)
	      ShortestDistance = SurfaceVec[a]->Distance();
	  }
    
  }
  
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  
  return kInfinity; 
}


G4double G4BREPSolidPolyhedra::DistanceToOut(register const G4ThreeVector& Pt,
					     register const G4ThreeVector& V,
					     const G4bool calcNorm, 
					     G4bool *validNorm, 
					     G4ThreeVector *n) const
{
  if(validNorm)
    *validNorm = false;

  Reset();  

  G4double halfTolerance = kCarTolerance*0.5;  
  G4Point3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  
  //  G4double kInfinity = 10e20;
  G4Ray r(Pttmp, Vtmp);
  TestSurfaceBBoxes(r);
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  ShortestDistance = kInfinity;
  
  for(int a=0; a< nb_of_surfaces;a++)
  {
    if(SurfaceVec[a]->Active())
      if(SurfaceVec[a]->Intersect(r))
	if(ShortestDistance > SurfaceVec[a]->Distance()&&
	   SurfaceVec[a]->Distance() > halfTolerance)
	  ShortestDistance = SurfaceVec[a]->Distance();
  }
  
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  
  return kInfinity; 
}


G4double G4BREPSolidPolyhedra::DistanceToOut(const G4ThreeVector& Pt) const
{
  G4double Dist=kInfinity;
  G4double tmpdist=kInfinity;  
  G4double halfTolerance = kCarTolerance*0.5;    
  
  for(int a=0; a< nb_of_surfaces;a++)
  {
    tmpdist = fabs(SurfaceVec[a]->HowNear(Pt));
  
    if(Dist>tmpdist  && tmpdist> halfTolerance) 
      Dist = tmpdist;
  }
  
  return Dist;
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



