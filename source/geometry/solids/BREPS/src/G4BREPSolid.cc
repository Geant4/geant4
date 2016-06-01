#include "G4BREPSolid.hh"
#include "G4AffineTransform.hh"
#include "G4VGraphicsScene.hh"
#include "G4Polyhedron.hh"
#include "G4NURBSbox.hh"
#include "G4BoundingBox3D.hh"
#include "G4FPlane.hh"
#include "G4BSplineSurface.hh"
#include "G4ToroidalSurface.hh"
#include "G4SphericalSurface.hh"

G4Ray G4BREPSolid::Track;
G4double G4BREPSolid::ShortestDistance=1000000;
int G4BREPSolid::NumberOfSolids=0;
InstMgr G4BREPSolid::InstanceList;


G4BREPSolid::G4BREPSolid(const G4String name) : G4VSolid(name)
{
  place=0;
  Box=0;
  Convex=0;
  AxisBox=0;
  PlaneSolid=0;
  active=1;
  intersectionDistance=kInfinity;
  startInside=0;
  solidname = name;
}

G4BREPSolid::G4BREPSolid( const G4String name        , 
			  G4Surface**    srfVec      , 
			  G4int          numberOfSrfs  ) : G4VSolid(name)
{
  place                = 0;
  Box                  = 0;
  Convex               = 0;
  AxisBox              = 0;
  PlaneSolid           = 0;
  active               = 1;
  intersectionDistance = kInfinity;
  startInside          = 0;
  nb_of_surfaces       = numberOfSrfs;
  SurfaceVec           = srfVec;

  Initialize();
}



G4BREPSolid::~G4BREPSolid()
{
  if(place)
    delete place;
  
  delete bbox;
  
  for(G4int a=0;a<nb_of_surfaces;a++)
    delete SurfaceVec[a];
  
  delete [] SurfaceVec;
}


void G4BREPSolid::Initialize()
{
  if(active)
  {
    // Calc bounding box for solids and surfaces
    // Convert concave planes to convex     
    ShortestDistance=1000000;
    IsBox();	
    CheckSurfaceNormals();
    
    if(!Box || !AxisBox)
      IsConvex();
    
    CalcBBoxes();
  }
}


void G4BREPSolid::CheckSurfaceNormals()
{
  if(!PlaneSolid)
    return; // All faces must be planar
 
  Convex=1;
 
  // Checks that the normals of the surfaces point outwards.
  // If not, turns the Normal to point out.
  
  // Loop through each face and check the G4Vector3D of the Normal
  G4Surface* srf;
  G4Vector3D *Normal1, Normal2;
  G4Point3D V;
  
  G4int PointNum=0;
  G4int SrfNum = 0;
  G4double YValue=0;
  G4Point3D Pt;
  
  //    const int Faces = all_surfaces.number_of_elements;
  //    const int Faces = surfaces.entries();
  G4int a, b;
  for(a=0; a<nb_of_surfaces; a++)
  {
    //Find vertex point containing extreme y value
    //	srf = all_surfaces.GetSurface(a);
    //	srf = surfaces(a);
    srf = SurfaceVec[a];		
    int Points = srf->GetNumberOfPoints();	
    
    for(b =0; b<Points; b++)
    {
      Pt = (G4Point3D)srf->GetPoint(b);	    
      if(YValue < Pt.y())
      {
	YValue = Pt.y();
	PointNum = b;// Save point number
	SrfNum = a;  // Save srf number
      }
    }
  }

  // Move the selected face to the first in the List
  //    srf = all_surfaces.GetSurface(SrfNum);
  //    srf = surfaces(SrfNum);
  srf = SurfaceVec[SrfNum];            
  
  //    all_surfaces.MoveToFirst(srf);

  // Start handling the surfaces in order and compare
  // the neighbouring ones and turn their normals if they
  // point inwards
  G4Point3D Pt1;
  G4Point3D Pt2;
  G4Point3D Pt3;
  G4Point3D Pt4;
  
  int ConnectingPoints=0;
  
  G4Vector3D N1;
  G4Vector3D N2;    
  G4Vector3D N3;    
  G4Vector3D N4;    

  int* ConnectedList = new int[nb_of_surfaces];    

  for(a=0; a<nb_of_surfaces; a++)
    ConnectedList[a]=0;
  
  int Connections=0;
  
  G4Surface* ConnectedSrf;

  for(a=0; a<nb_of_surfaces-1; a++)
  {
    
    if(ConnectedList[a] == 0) 
      break;
    else	
      ConnectedList[a]=1;
    
    srf = SurfaceVec[a];
    int SrfPoints = srf->GetNumberOfPoints();
    N1 = (srf->Norm())->GetDir();

    for(b=a+1; b<nb_of_surfaces; b++)
    {
      if(ConnectedList[b] == 1) 
	break;
      else	
	ConnectedList[b]=1;
     
      // Get next in List
      //	    ConnectedSrf = all_surfaces.GetSurface(b);
      //	    ConnectedSrf = surfaces(b);
      ConnectedSrf = SurfaceVec[b];	    	    
	
      // Check if it is connected to srf by looping through the
      // points.
      int ConnSrfPoints = ConnectedSrf->GetNumberOfPoints();

      for(G4int c=0;c<SrfPoints;c++)
      {
	Pt1 = srf->GetPoint(c);
	
	for(G4int d=0;d<ConnSrfPoints;d++)
	{
	  // Find common points
	  Pt2 = (ConnectedSrf)->GetPoint(d);
	 
	  if( Pt1 == Pt2 ) 
	  {
	    // Common point found.
	    
	    // Compare normals
	    N2 = ((ConnectedSrf)->Norm())->GetDir();
	    
	    // Check cross product.
	    G4Vector3D CP1 = N1.cross(N2);
	    G4double CrossProd1 = CP1.x()+CP1.y()+CP1.z();
	    
	    // Create the other normals
	    if(c==0) 
	      Pt3 = srf->GetPoint(c+1);
	    else
	      Pt3 = srf->GetPoint(0);
	    
	    N3 = (Pt1-Pt3);
	    
	    if(d==0) 
	      Pt4 = (ConnectedSrf)->GetPoint(d+1);
	    else
	      Pt4 = (ConnectedSrf)->GetPoint(0);
	    
	    N4 = (Pt1-Pt4);
		      
	    G4Vector3D CP2 = N3.cross(N4);
	    G4double CrossProd2 = CP2.x()+CP2.y()+CP2.z();
	    
	    G4cout << "\nCroosProd2: " << CrossProd2;
	    
	    if( (CrossProd1 < 0 && CrossProd2 < 0) ||
		(CrossProd1 > 0 && CrossProd2 > 0)    )
	      {
		// Turn Normal
		(ConnectedSrf)->Norm()
		  ->SetDir(-1 * (ConnectedSrf)->Norm()->GetDir());
		
		// Take the CrossProd1 again as the other Normal was turned.
		CP1 = N1.cross(N2);
		CrossProd1 = CP1.x()+CP1.y()+CP1.z();			    
	      }
	    
	    if(CrossProd1 > 0) 
	      Convex=0;
	  }
	}
      }
    }
  }
  
  delete []ConnectedList;
}


int G4BREPSolid::IsBox()
{
  // This is Done by checking that the solid consist of 6 planes.
  // Then the type is checked to be planar face for each face.
  // For each G4Plane the Normal is computed. The dot product
  // of one face Normal and each other face Normal is computed.
  // One result should be 1 and the rest 0 in order to the solid
  // to be a box.

  Box=0;
  G4Surface* srf1, *srf2;
  register G4int a;
  
  // Calc the Normal for the planes
  for(a=0; a < nb_of_surfaces;a++)    
  {
    srf1 = SurfaceVec[a];		
    
    if(srf1->MyType()==1)
      (srf1)->Project();     // Calc the projection
    else
    {
      PlaneSolid=0;
      return 0;
    }
  }

  // Check that all faces are planar
  for(a=0; a < nb_of_surfaces;a++)    
  {
    srf1 = SurfaceVec[a];		
    
    if (srf1->MyType()!=1)
      return 0;
  }

  PlaneSolid = 1;
  
  // Check that the amount of faces is correct
  if(nb_of_surfaces!=6) return 0;
  
  G4Point3D Pt;
  G4int Points;
  G4int Sides=0;
  G4int Opposite=0;

  srf1 = SurfaceVec[0];
  Points = (srf1)->GetNumberOfPoints();
  
  if(Points!=4)
    return 0;
    
  G4Vector3D Normal1 = (srf1->Norm())->GetDir();
  G4double Result;
  
  for(G4int b=1; b < nb_of_surfaces;b++)
  {
    srf2 = SurfaceVec[b];			
    G4Vector3D Normal2 = ((srf2)->Norm())->GetDir();
    Result = fabs(Normal1 * Normal2);
    
    if((Result != 0) && (Result != 1))
      return 0;
    else
    {
      if(!(int)Result) 
	Sides++;
      else
	if(((int)Result) == 1)
	  Opposite++;
    }
  }

  if((Opposite != 1) && (Sides != nb_of_surfaces-2))
    return 0;
  
  G4Vector3D x_axis(1,0,0);
  G4Vector3D y_axis(0,1,0);
  
  if(((fabs(x_axis * Normal1) == 1) && (fabs(y_axis * Normal1) == 0)) ||
     ((fabs(x_axis * Normal1) == 0) && (fabs(y_axis * Normal1) == 1)) ||
     ((fabs(x_axis * Normal1) == 0) && (fabs(y_axis * Normal1) == 0)))
    AxisBox=1;
  else 
    Box=1;
  
  return 1;
}	


G4bool G4BREPSolid::IsConvex()
{
  if(!PlaneSolid)
    return 0; // All faces must be planar

  // This is not robust. There can be concave solids
  // where the concavity comes for example from
  // three triangles.

  // Additional  checking 20.8. For each face the connecting faces are
  // found and the cross product computed between the face and each
  // connecting face. If the result changes value at any point the
  // solid is concave.
  
  G4Surface* Srf;
  G4Surface* ConnectedSrf;
  int Result;
  Convex = 1;
  
  G4int a, b, c, d;
  for(a=0;a<nb_of_surfaces;a++)
  {
    Srf = SurfaceVec[a];		
    
    // Primary test. Test wether any one of the faces
    // is concave -> solid is concave. This is not enough to
    // distinguish all the cases of concavity.
    Result = Srf->IsConvex();
    
    if(Result != -1)
    {
      Convex = 0;
      return 0;
    }	
  }

  Srf = SurfaceVec[0];        
  G4Point3D Pt1;
  G4Point3D Pt2;
  
  int ConnectingPoints=0;
  
  G4Vector3D N1;
  G4Vector3D N2;    

  // L. Broglia
  // The number of connecting points can be 
  // (nb_of_surfaces-1) * nb_of_surfaces  (loop a & loop b)
  // HandledList is not used : why ?
  
  // G4int* ConnectedList = new G4int[nb_of_surfaces];  
  G4int* ConnectedList = new G4int[(nb_of_surfaces-1) * nb_of_surfaces];  
  G4int* HandledList = new G4int[nb_of_surfaces];    

  for(a=0; a<nb_of_surfaces; a++)
  {
    HandledList[a]=0;
    ConnectedList[a]=0;
  }
  
  HandledList[0]=1;
  G4int Connections=0;

  for(a=0; a<nb_of_surfaces-1; a++)
  {
    Srf = SurfaceVec[a];	
    G4int SrfPoints = Srf->GetNumberOfPoints();
    Result=0;
    
    for(b=0; b<nb_of_surfaces; b++)
    {
      if(b==a)
	b++;
      
      if(b==nb_of_surfaces)
	break;
      
      // Get next in List
      ConnectedSrf = SurfaceVec[b];
      
      // Check if it is connected to Srf by looping through the
      // points.
      G4int ConnSrfPoints = ConnectedSrf->GetNumberOfPoints();
      
      for(c=0; c<SrfPoints; c++)
      {
	const G4Point3D& Pts1 =Srf->GetPoint(c);
	
	for(d=0; d<ConnSrfPoints; d++)
	{
	  // Find common points
	  const G4Point3D& Pts2 = ConnectedSrf->GetPoint(d);
	  
	  if(Pts1 == Pts2)
	    ConnectingPoints++;
	}
	
	if(ConnectingPoints > 0) 
	  break;
      }
      
      if( ConnectingPoints > 0 )
      {
	Connections++;
	ConnectedList[Connections]=b;
      }
      
      ConnectingPoints=0;
    }
  }
  
  // If connected, check for concavity.
  // Get surfaces from ConnectedList and compare their normals
  for(c=0; c<Connections; c++)
  {
    G4int Left=0; 
    G4int Right =0;
    G4int tmp = ConnectedList[c];
    
    Srf = SurfaceVec[tmp];
    ConnectedSrf = SurfaceVec[tmp+1];
    
    // Get normals.
    N1 = Srf->Norm()->GetDir();	    
    N2 = ConnectedSrf->Norm()->GetDir();
    // Check cross product.
    G4Vector3D CP = N1.cross(N2); 
    G4double CrossProd = CP.x()+CP.y()+CP.z();
    
    if( CrossProd > 0 )	
      Left++;
    
    if(CrossProd < 0)
      Right++;
    
    if(Left&&Right)
    {
      Convex = 0;		
      return 0;
    }
        
    Connections=0;
  } 
  
  Convex=1;
  // L. Broglia
  // Problems with this delete when there are many solids to create
  // delete []ConnectedList;
  // delete []HandledList;    
  return 1;
}


G4bool G4BREPSolid::CalculateExtent(const     EAxis pAxis,
				    const     G4VoxelLimits& pVoxelLimit,
				    const     G4AffineTransform& pTransform,
				    G4double& pMin, 
				    G4double& pMax           ) const
{
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();
  
  G4ThreeVector VMin(Min.x(),Min.y(),Min.z());
  G4ThreeVector VMax(Max.x(),Max.y(),Max.z());  
  
  switch (pAxis)
  {
  case kXAxis:
    pMin=Min.x();
    pMax=Max.x();
    break;

  case kYAxis:
    pMin=Min.y();
    pMax=Max.y();
    break;

  case kZAxis:
    pMin=Min.z();
    pMax=Max.z();
    break;
  }
  
  pMin-=kCarTolerance;
  pMax+=kCarTolerance;
  
  return true;

  /*
     if (!pTransform.IsRotated())
     {
     //      if(pVoxelLimit.Inside(VMin))
     //	if(pVoxelLimit.Inside(VMax))
     {
     switch (pAxis)
     {
     case kXAxis:
     pMin=VMin.x();
     pMax=VMax.x();
     break;
     case kYAxis:
     pMin=VMin.y();
     pMax=VMax.y();
     break;
     case kZAxis:
     pMin=VMin.z();
     pMax=VMax.z();
     break;
     }

     //	    pMin-=kCarTolerance;
     //	    pMax+=kCarTolerance;
     return true;
     }
     //      else
     {
     //	  return false;
     }
     }
     //  else
     {
     //      return false;
     }

     return false;*/
}


EInside G4BREPSolid::Inside(register const G4ThreeVector& Pt)const
{

  G4cout<<"\n Solid Id="<<GetId();

  // Tolerance has also to be considered i.e. the cases
  // where the Point is very close to the box. -1 is 
  // returned when the Point is within the boundary.
  Reset();  

  // Get the bounding box extents
  const G4Point3D min = bbox->GetBoxMin();
  const G4Point3D max = bbox->GetBoxMax();
  
  // First check if the point is Inside the bbox
  // of the solid.
  if((Pt.x() < min.x() || Pt.x() > max.x())||
     (Pt.y() < min.y() || Pt.y() > max.y())||
     (Pt.z() < min.z() || Pt.z() > max.z()))    
    return kOutside;

  // If the point is Inside the bbox, it is possibly Inside 
  // the solid
  
  // Create the ray.
  G4Vector3D v(0,0,1);
  Track.Init( Pt, v);
  
  // Check the bboxes of the surfaces first to get the favorable 
  // surfaces to process.
  TestSurfaceBBoxes( Track );

  // At this point it might be useful to have a routine that checks
  // if the point is outside all the face bboxes, Inside any of them
  // or Inside all of them. The last case would mean Inside the solid for
  // sure. Also should be considered if various types of solids should be 
  // handled separately i.e. convex polyhedras versus b-spline solids.
  int Hits = 0;
  
  // Continue by shooting an arbitrary ray (G4Vector3D = v) starting from Point
  // Repeat this until no more hits are found. The amount of hits gives the 
  // answer (par : outside; odd : inside)
  // Tolerance is a problem as there might be cases where we don't
  // get the actual intersection on the surface and the Step is not long
  // enough to jump over the real boundary. 
  G4double Tolerance = 0.0001;
  G4double Dist=Tolerance+1;
  
  const G4Vector3D& RayDir = Track.GetDir();    
  
  while(Intersect( Track ))    
    if(FinalEvaluation( Track ))      
    {
      // Step over the tolerance and shoot again in the same
      // direction using the Hit point as a starting point
      // for the new ray 
      Hits++;	
	
      // Tolerance Step
      if(Hits==1) 
	Dist=ShortestDistance;
      
      const G4Point3D& NewStart = intersection_point+(Tolerance*RayDir);
      
      // Set ray starting point
      Track.SetStart(NewStart);	  
      
      Reset();
      
      // To make the routine more robust a check could be made that
      // Hits doesn't grow too much say > 10. If this happens
      // v could be changed to be something else. This may avoid
      // situations where the direction chosen for the test
      // is not plausible since at that direction there are many
      // faces very close to each other or a face is highly complex or 
      // the direction is in the direction of the face.
      //if(Hits > 10)
      // 	{
      //	r.SetStart(Point);Hits = 0;
      //	r.set_G4Vector3D(G4Vector3D(1,0,0));
      //	while(first_intersect( r ))
      //	{
      //		Hits++;	
      //		intersection_point.X() += Tolerance;
      //	      	r.SetStart(intersection_point);
      //	}
      //	}
    }
    // else
    //  break;
  
  // Set the surfaces to active again
  G4Surface* Srf;
  
  for(G4int a=0;a<nb_of_surfaces;a++)
  {
    Srf = SurfaceVec[a];
    Srf->Reset();
  }

  //    G4cout << "\n Inside::Hits = " << Hits << "\n";
  if(Dist < (0.5*kCarTolerance))
    return kSurface;
  
  // Take the mod of hits with 2
  if(Hits&1)
    // Inside
    return kInside;
  else
    // Outside
    return kOutside;
}


G4ThreeVector G4BREPSolid::SurfaceNormal(const G4ThreeVector& GetStart)const
{  
  return GetStart;;
}


G4double G4BREPSolid::DistanceToIn(const G4ThreeVector& Pt)const
{
  // This function calculates the approximative distance 
  // of a point from a surface. For exemple, the DistanceToIn
  // from a point to a G4FPlane is the closest distance

  G4cout<<"\n Solid Id="<<GetId();
  
  Reset();  
  
  G4double halfTolerance = kCarTolerance*0.5;
  G4Surface* srf;
  G4double PointDistance = INFINITY;
  G4double TmpDistance = 0;
  
  for(G4int a=0; a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];
    TmpDistance = fabs(srf->ClosestDistanceToPoint(Pt));

    if( (TmpDistance < PointDistance) && (TmpDistance > halfTolerance) )
      PointDistance = TmpDistance;
    else if(TmpDistance < halfTolerance)
      PointDistance = 0; // the point is on the surface, 
                         // but maybe not into the boundary 
  
    if(PointDistance < ShortestDistance)
      ShortestDistance = PointDistance;
  }
  
  return ShortestDistance;
}


G4double G4BREPSolid::DistanceToIn(register const G4ThreeVector& Pt, 
				   register const G4ThreeVector& V  )const
{
  G4cout<<"\n Solid Id="<<GetId();
  
#ifdef G4VERBOSE
  if(V.mag2() == 0.0)
    G4Exception("Error in G4BREPSolid::DistanceToIn(Pt, Vec) : Vec = 0");
#endif
  
  Reset();
  G4double halfTolerance = kCarTolerance*0.5;  
  Track.Init(G4Point3D(Pt),V);
  
  
  if(!bbox->Test(Track))	
    return kInfinity;

/*  
  if(AxisBox)
  {
    ShortestDistance = bbox->GetDistance();
    
    if(ShortestDistance > halfTolerance)
      return ShortestDistance;

    //  L. Broglia
    // Do not see this surface
    //else
      // ShortestDistance = kInfinity;
      //return 0; // the point is on the surface
  }
*/

  TestSurfaceBBoxes(Track);
  
  if(Intersect( Track ))
    if(FinalEvaluation(Track,1))
      return sqrt(ShortestDistance);
 
  return kInfinity;
  
  /*
     G4Vector3D Pttmp(Pt);
     G4Vector3D Vtmp(V);   
     //  G4double kInfinity = 10e20;
     G4Ray r(Pttmp, Vtmp);
     if(!bbox->Test3dBBox(r))	
     return kInfinity;
     if(AxisBox)
     {
     ShortestDistance = bbox->GetDistance();
     if(ShortestDistance > halfTolerance)
     return ShortestDistance;
     else
     ShortestDistance = kInfinity;	  
     }
     //  if(number_of_Surfaces>75)
     //    RemoveHiddenFaces(r, 1);
     TestSurfaceBBoxes(r);
     if(Intersect( r ))
     if(FinalEvaluation(r,1))
     return sqrt(ShortestDistance);
     return kInfinity;
     */
}


G4double G4BREPSolid::DistanceToOut(const G4ThreeVector& Pt)const
{

  G4cout<<"\n Solid Id="<<GetId();
  
  Reset();
  
  G4double halfTolerance = kCarTolerance*0.5;  
  G4Surface* srf;
  G4double PointDistance = INFINITY;
  G4double TmpDistance = 0;
  
  for(G4int a=0; a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];
    TmpDistance = fabs(srf->ClosestDistanceToPoint(Pt));
    
    if( (TmpDistance < PointDistance) && (TmpDistance > halfTolerance) )
      PointDistance = TmpDistance; 
    else if(TmpDistance < halfTolerance)
      PointDistance = 0 ; // the point is on the surface
  }
  
  if(PointDistance < ShortestDistance) 
    ShortestDistance = PointDistance;
 
  return ShortestDistance;
}


G4double G4BREPSolid::DistanceToOut(register const G4ThreeVector& P, 
				    register const G4ThreeVector& D, 
				    const G4bool calcNorm, 
				    G4bool *validNorm, 
				    G4ThreeVector *n                ) const
{
  G4cout<<"\n Solid Id="<<GetId();
  
  if(validNorm)
    *validNorm = false;

#ifdef G4VERBOSE
  if(D.mag2() == 0.0)
    G4Exception("Error in G4BREPSolid::DistanceToOut(Pt, Vec) : Vec = 0");
#endif
     
  Reset();
  
  Track.Init(G4Point3D(P),D);
 
  TestSurfaceBBoxes(Track);
 
  if(Intersect(Track))
    if(FinalEvaluation(Track))
      if(ShortestDistance > kCarTolerance/2) // if d=0, do not see the surface
	return sqrt(ShortestDistance);
  
  return kInfinity; // This should never happen
}


void G4BREPSolid::DescribeYourselfTo (G4VGraphicsScene& scene) const 
{
  scene.AddThis (*this);
}


G4VisExtent G4BREPSolid::GetExtent() const 
{
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();  
  return G4VisExtent (Min.x(), Max.x(), Min.y(), Max.y(), Min.z(), Max.z());
}


G4Polyhedron* G4BREPSolid::CreatePolyhedron () const
{
  // temporary
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();  

  return new G4PolyhedronBox (Max.x(), Max.y(), Max.z());
}


G4NURBS* G4BREPSolid::CreateNURBS () const 
{
  // temporary
  G4Point3D Min = bbox->GetBoxMin();
  G4Point3D Max = bbox->GetBoxMax();  

  return new G4NURBSbox (Max.x(), Max.y(), Max.z());
}


int G4BREPSolid::CreateSTEPData()
{
  // create the solid entity
  //STEPentity* ent = new STEPentity();
  //stateEnum *sEnu = new stateEnum("newSE");
  //MgrNode *mnode = new MgrNode(ent, sEnu);

  // create the attributelist & attributes for this solid
  //    STEPattributeList *aList = new STEPattributeList();
  //AttrDescriptor *aDesc = new AttrDescriptor();
  //  STEPattribute *name_attr = new STEPattribute(aDesc, name);
  
  // create the mgrnode neede by the instance List

  // append node to instance List
  
  // call entoity creation routines for child entities.
  return 0; // to shut up compilers
}


void G4BREPSolid::CalcBBoxes()
{
  // First initialization
  // Calculates the bounding boxes for the surfaces and
  // for the solid.

  G4Surface* srf;
  register G4Point3D min, max;
  
  if(active)
  {
    min =  PINFINITY;
    max = -PINFINITY;
    
    for(G4int a = 0;a < nb_of_surfaces;a++)
    {
      // Get first in List
      srf = SurfaceVec[a];
      G4int convex=1; 
      G4int concavepoint=-1;
	
      if (srf->MyType() == 1) 
      {
	concavepoint = srf->IsConvex();
	convex = srf->GetConvex();
      }
      // Make bbox for face
      //	    if(convex && Concavepoint==-1)
      {
	srf->CalcBBox();
	G4Point3D box_min = srf->bbox->GetBoxMin();
	G4Point3D box_max = srf->bbox->GetBoxMax();
	// Find max and min of face bboxes to make 
	// solids bbox.
	
	// max < box_max;
        if(max.x() < box_max.x()) max.setX(box_max.x()); 
        if(max.y() < box_max.y()) max.setY(box_max.y()); 
        if(max.z() < box_max.z()) max.setZ(box_max.z()); 
	
        // min > box_min;
        if(min.x() > box_min.x()) min.setX(box_min.x()); 
        if(min.y() > box_min.y()) min.setY(box_min.y()); 
        if(min.z() > box_min.z()) min.setZ(box_min.z());
      }
    }
    
    bbox =  new G4BoundingBox3D(min, max);
    //  G4cout << "\nBox " << min.X() << "  " << min.Y() << " " << min.Z();
    //G4cout << "\n--- " << max.X() << "  " << max.Y() << " " << max.Z();
    
    return;
  }
  
  G4cout << "\n No bbox calculated for solid. Error.";
}


void G4BREPSolid::RemoveHiddenFaces(register const G4Ray& rayref, int In) const
{
  // Deactivates the planar faces that are
  // on the "back" side of a solid.
  // B-splines are not handled by this routine. Also cases
  // where the ray starting point is Inside the bbox of the solid
  // are ignored as we don't know if the starting point is Inside the actual 
  // solid except for axisoriented boxlike solids
  
  register G4Surface* srf;
  register const G4Vector3D& RayDir = rayref.GetDir();
  register G4double Result;
  G4int a;
  //    if(!AxisBox)
  // In all other cases the ray starting point is outside the solid

  if(!In)// In all other cases the ray starting point is outside the solid
    for(a=0; a<nb_of_surfaces; a++)
    {
      // Deactivates the solids faces that are hidden
      srf = SurfaceVec[a];
      
      if(srf->MyType()==1)
      {
	const G4Vector3D& Normal = (srf->Norm())->GetDir();
	Result = (RayDir * Normal);
	
	if( Result >= 0 )
	  srf->Deactivate();
      }
    }
  else
    for(a=0; a<nb_of_surfaces; a++)
    {
      // Deactivates the AxisBox type solids faces whos normals 
      // point in the G4Vector3D opposite to the rays G4Vector3D
      // i.e. are behind the ray starting point as in this case the
      // ray starts from Inside the solid.
      srf = SurfaceVec[a];
      
      if(srf->MyType()==1)
      {
	const G4Vector3D& Normal = (srf->Norm())->GetDir();
	Result = (RayDir * Normal);
	
	if( Result < 0 )
	  srf->Deactivate();
      }
    }
}


void G4BREPSolid::TestSurfaceBBoxes(register const G4Ray& rayref) const 
{
  register G4Surface* srf;
  G4int active_srfs = nb_of_surfaces;
  
  // Do the bbox tests to all surfaces in List
  // for planar faces the intersection is instead evaluated.
  G4int intersection=0;
  
  for(G4int a=0;a<nb_of_surfaces;a++)
  {
    // Get first in List
    srf = SurfaceVec[a];
    
    if(srf->Active())
    {
      // Get type
      if(srf->MyType() != 1) // 1 == planar face
      {
	if(srf->bbox->Test(rayref))
	  srf->Distance(bbox->GetDistance());	
	else
	{
	  // Test failed. Flag as inactive.
	  srf->Deactivate();
	  active_srfs--;
	}
      }
      else
      {
	// Type was convex planar face
	intersection = srf->Intersect(rayref);
	
	if(!intersection) 
	  active_srfs--;
      }
    }
    else
      active_srfs--;		    
  }
  
  if(!active_srfs) Active(0);
}


int G4BREPSolid::Intersect(register const G4Ray& rayref) const
{
  // Gets the roughly calculated closest
  // intersection point for a b_spline & accurate point for others
  register G4Surface* srf;
  G4double HitDistance = -1;
  const G4Point3D& RayStart = rayref.GetStart();
  const G4Point3D& RayDir   = rayref.GetDir();

  G4int result=1;
   
  // Sort List of active surfaces according to
  // bbox distances to ray starting point.
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  G4int Number=0;
   
  // Start handling active surfaces in order
  for(register G4int a=0;a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];	
    int included = 0;
    
    if(srf->Active())
    {
      result = srf->Intersect(rayref);
      if(result)
      {
	register G4Surface* tmp;
	   
	// L. Broglia
	// What is the utility of this test ?
	// If only 1 surface is intersected, it return 0
	// instead the intersection exist
 
	/*
	if((a+1)<nb_of_surfaces)
	  tmp = SurfaceVec[a+1];		
	else
	  tmp = 0;
	
	if((tmp)  && ( tmp->Active()))
	{      
	*/

	// Get the evaluated point on the surface
	G4Point3D& closest_point = srf->closest_hit;

	// Test for DistanceToIn(pt, vec)
	// if d = 0 and vec.norm > 0, do not see the surface
	if( !( (srf->Distance() < kCarTolerance/2)  || 
	     (RayDir.dot(srf->SurfaceNormal(closest_point)) > 0) ) )
	{
  
	  if(srf->MyType()==1)
	    HitDistance = srf->Distance();
	  else
	  {
	    // Check if the evaluated point is in front of the 
	    // bbox of the next surface.
	    // Took sqrt away to gain speed. Squaredistances may
	    // as well be used.	 
	    HitDistance = RayStart.distance2(closest_point);
	  }
	  
	  /*
	    included = 0;
	    
	    
	    if(tmp->MyType()==1)
	    {
	    if(HitDistance >= tmp->Distance())
	    included = 1;
	    }
	    else
	    {
	    G4double Dist = tmp->bbox->GetDistance();
	    Dist = Dist*Dist;
	    
	    if(HitDistance >= Dist)		    
	    included = 1;
	    }
	    
	    if(included) // Check also which other surfaces it is included in
	    {
	    if(Number+2<nb_of_surfaces)
	    for(G4int a=Number+2; a < nb_of_surfaces; a++)
	    {
	    tmp = SurfaceVec[a];
	    
	    if(HitDistance < tmp->Distance())
	    tmp->Deactivate();
	    }
	    }
	    else 
	    {
	    // Mark rest surfaces as inactive
	    for(register int c=a+1;c<nb_of_surfaces;c++)
	    {
	    tmp = SurfaceVec[c];
	    
	    if(tmp->MyType()!=1)
	    tmp->Deactivate();
	    }
	    } 
	    }*/
	}
      }   
      else
      {
	// No Hit.
	included = 1;
	srf->Deactivate();
      }//if(result...
    }//if(srf->act...
    
    Number++;
  } // while...
  
  if(HitDistance < 0)
    return 0;
  
  QuickSort(SurfaceVec, 0, nb_of_surfaces-1);
  
  if(!(SurfaceVec[0]->Active()))
    return 0;     
  
  ((G4BREPSolid*)this)->intersection_point = SurfaceVec[0]->closest_hit;

  
  bbox->SetDistance(HitDistance);

  return 1;
}


int G4BREPSolid::FinalEvaluation(register const G4Ray& rayref, 
				 const int ToIn               ) const
{
  register G4Surface* srf;
  G4double halfTolerance = 0.5*kCarTolerance;
  G4double Dist=0;
  G4int count=0;
  ((G4BREPSolid*)this)->intersectionDistance = kInfinity;
  
  for(register G4int a=0;a<nb_of_surfaces;a++)
  {
    srf = SurfaceVec[a];
    
    if(srf->Active())
    {
      const G4Point3D& srf_intersection = srf->Evaluation(rayref);
      
      // Calc Hit point distance from ray starting point.
      if(srf->MyType() != 1)
      {
	// took sqrt away to gain speed...
	// Square distances are used instead
	
	G4Point3D start = rayref.GetStart();
	Dist = srf_intersection.distance2(start);  
      }
      else 
	Dist = srf->Distance(); 
	
      // Skip point wichare on the surface i.e. within
      // tolerance of the surface
      // Special handling for DistanceToIn & reflections 
      if(sqrt(Dist) < halfTolerance)
      {
	if(ToIn) 
	{
	  const G4Vector3D& Dir = rayref.GetDir();
	  const G4Point3D& Hit = srf->closest_hit;
	  const G4Vector3D& Norm = srf->SurfaceNormal(Hit);
	   
	  if(( Dir * Norm ) >= 0)
	  {
	    Dist = INFINITY;
	    srf->Deactivate();
	  }
	  
	  // else continue with the distance,
	  // even though < tolerance
	}
	else
	{
	  Dist = INFINITY;
	  srf->Deactivate();
	}
      }
      
      // If more than one surfaces are evaluated til the
      // final stage, only the closest point is taken
      if(Dist < intersectionDistance)
      {
	// Check that Hit is in the direction of the ray 
	// from the starting point
	const G4Point3D& Pt = rayref.GetStart();
	const G4Vector3D& Dir = rayref.GetDir();
	
	G4Point3D TestPoint = (0.00001*Dir) + Pt;
	G4double TestDistance = srf_intersection.distance2(TestPoint);
	  
	if(TestDistance > Dist)
	{
	  // Hit behind ray starting point, no intersection.
	  Dist = INFINITY;
	  srf->Deactivate();
	}
	else
	{
	  ((G4BREPSolid*)this)->intersectionDistance = Dist;
	  ((G4BREPSolid*)this)->intersection_point = srf_intersection;
	}
	  
	// Check that the intersection is closer than the
	// next surfaces approx point.
	if(srf->Active())
	{
	  if(a+1<nb_of_surfaces)
	  {
	    const G4Vector3D& Dir = rayref.GetDir();
	    const G4Point3D& Hit = srf->closest_hit;
	    const G4Vector3D& Norm = srf->SurfaceNormal(Hit);
	      
	    // L. Broglia
	    //if(( Dir * Norm ) >= 0)
	    if(( Dir * Norm ) < 0)
	    {
	      Dist = INFINITY;
	      srf->Deactivate();
	    }
	    
	    // else continue with the distance,
	    // even though < tolerance

	    // L. Broglia
            // I think that this line has been forgeted
	    ShortestDistance = Dist;
	  }
	  else
	  {
	    ShortestDistance = Dist;			    
	    return 1;
	  }
	}
    
      }//if(Dist...
      
    }//if srf->Active()
    else
    {
     /* if(intersectionDistance < kInfinity)
	return 1;
      return 0;*/
    }
  }//for...
  
  if(intersectionDistance < kInfinity)    
    return 1;
  
  return 0;
}
 

G4Point3D G4BREPSolid::Scope()
{
  G4Point3D scope;
  G4Point3D Max = bbox->GetBoxMax();
  G4Point3D Min = bbox->GetBoxMin();  
  
  scope.setX(fabs(Max.x()) - fabs(Min.x()));
  scope.setY(fabs(Max.y()) - fabs(Min.y()));
  scope.setZ(fabs(Max.z()) - fabs(Min.z()));	  
  
  return scope;
}

















