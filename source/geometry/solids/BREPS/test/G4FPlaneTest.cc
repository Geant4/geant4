//////////////////////////////////////////////////////////////////////////
// $Id: G4FPlaneTest.cc,v 1.4 2000-08-28 08:58:04 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//////////////////////////////////////////////////////////////////////////
//
// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// Test the G4FPlane class
// Created by L. Broglia, 22 October 1998

#include "G4FPlane.hh"
#include "G4Surface.hh"
#include "G4Axis2Placement3D.hh"



int main()
{

/////////////////////////////////////////////////////////////
//
// I want to compare the 2 different creator.
// I utilized a piece of the G4BREPSolidPolyhedra.cc file

  G4cout<<"\n\n//////////////////////////////////////////////////////";
  
  G4Vector3D      Dir         ;
  G4Vector3D      Ax          ;
  G4Point3D       Porg        ;
  G4Point3DVector PointList(4);

  G4Point3D   LocalOrigin(0.0, 0.0, 0.0);  
  G4Vector3D  TmpAxis(1, 0, 0);
  G4Vector3D  Axis(0, 0, 1);

  double RMIN0     = 2;
  double RMIN1     = 1;
  double Length    = 2;
  double PartAngle = pi/2;

  PointList[0] = LocalOrigin + (RMIN0 * TmpAxis);
  PointList[1] = LocalOrigin + (Length*Axis) + (RMIN1 * TmpAxis);
  
  TmpAxis.rotateZ(PartAngle);
  
  PointList[2] = LocalOrigin + (Length*Axis) + (RMIN1 * TmpAxis);
  PointList[3] = LocalOrigin + (RMIN0 * TmpAxis);	  
  
  // Messages output
  for(int i=0; i<4; i++)
    G4cout<<"\n Pt"<<i<<" :"
	  <<" x= "<<PointList[i].x()
	  <<" y= "<<PointList[i].y()
	  <<" z= "<<PointList[i].z();
  

  G4FPlane SurfaceVec1( &PointList );
  
  G4FPlane SurfaceVec2
    (  PointList[0] -  PointList[1]                                      ,
      (PointList[3] -  PointList[0]).cross(PointList[0] -  PointList[1]) ,
       PointList[0]                                                        );

  G4Axis2Placement3D Pplace1 = SurfaceVec1.GetPplace();
  G4Axis2Placement3D Pplace2 = SurfaceVec2.GetPplace();
  
  G4Vector3D Dir1  =  Pplace1.GetRefDirection() ;
  G4Vector3D Ax1   =  Pplace1.GetAxis()         ;
  G4Point3D  Porg1 =  Pplace1.GetLocation()     ;
 
  G4Vector3D Dir2  =  Pplace2.GetRefDirection() ;
  G4Vector3D Ax2   =  Pplace2.GetAxis()         ;
  G4Point3D  Porg2 =  Pplace2.GetLocation()     ;

  // Messages output
  G4cout<<"\n\n Dir1  : x= "<<Dir1.x()<<" y= "<<Dir1.y()<<" z= "<<Dir1.z();
  G4cout<<"\n Dir2  : x= "<<Dir2.x()<<" y= "<<Dir2.y()<<" z= "<<Dir2.z();

  G4cout<<"\n\n Ax1   : x= "<<Ax1.x()<<" y= "<<Ax1.y()<<" z= "<<Ax1.z();
  G4cout<<"\n Ax2   : x= "<<Ax2.x()<<" y= "<<Ax2.y()<<" z= "<<Ax2.z();

  G4cout<<"\n\n Porg1 : x= "<<Porg1.x()<<" y= "<<Porg1.y()<<" z= "<<Porg1.z();
  G4cout<<"\n Porg2 : x= "<<Porg2.x()<<" y= "<<Porg2.y()<<" z= "<<Porg2.z(); 

  G4cout<<"\n\n coordinate axis 1 : PX= "<<Pplace1.GetPX()
	<<" PY= "<<Pplace1.GetPY()
	<<" PZ= "<<Pplace1.GetPZ() ;

  G4cout<<"\n coordinate axis 2 : PX= "<<Pplace2.GetPX()
	<<" PY= "<<Pplace2.GetPY()
	<<" PZ= "<<Pplace2.GetPZ() ;


  G4Plane Plane1 = SurfaceVec1.GetPplane();
  G4Plane Plane2 = SurfaceVec2.GetPplane();

  G4cout<<"\n\n Plane1 : a= "<<Plane1.a
	<<" b= "<<Plane1.b
	<<" c= "<<Plane1.c
	<<" d= "<<Plane1.d   ;

  G4cout<<"\n Plane2 : a= "<<Plane2.a
	<<" b= "<<Plane2.b
	<<" c= "<<Plane2.c
	<<" d= "<<Plane2.d  ;

  G4Ray* Normal1 = SurfaceVec1.Norm();
  G4Ray* Normal2 = SurfaceVec2.Norm();

  G4Point3D  start1 = (*Normal1).GetStart(); 
  G4Point3D  start2 = (*Normal2).GetStart();
  G4Vector3D dir1   = (*Normal1).GetDir();
  G4Vector3D dir2   = (*Normal2).GetDir();

  G4cout<<"\n\n Normal 1" ;
  G4cout<<"\n Start : x= "<<start1.x()<<" y= "<<start1.y()<<" z= "<<start1.z();
  G4cout<<"\n Dir   : x= "<<dir1.x()<<" y= "<<dir1.y()<<" z= "<<dir1.z();

  G4cout<<"\n\n Normal 2" ;
  G4cout<<"\n Start : x= "<<start2.x()<<" y= "<<start2.y()<<" z= "<<start2.z();
  G4cout<<"\n Dir   : x= "<<dir2.x()<<" y= "<<dir2.y()<<" z= "<<dir2.z();

  // Now, I test the function ClosestDistanceToPoint
  G4Point3D Pout(2, 2, 4);
  G4double dist1 =  SurfaceVec1.ClosestDistanceToPoint(Pout);
  G4double dist2 =  SurfaceVec2.ClosestDistanceToPoint(Pout);
  
  G4cout<<"\n\n Distance 1 ="<<dist1;
  G4cout<<"\n Distance 2 ="<<dist2;



//////////////////////////////////////////////////////////////////////
//
// This test show that the creation of the plane not depend on the 
// points in the plane
//

  G4cout<<"\n\n//////////////////////////////////////////////////////";
  
  G4Point3D P0 = PointList[0];
  G4Point3D P1 = PointList[1];
  G4Point3D P2 = PointList[2];
  G4Point3D P3 = PointList[3];
  
  PointList[0] = P0;
  PointList[1] = P1; 
  PointList[2] = P2; 
  PointList[3] = P3;
  G4FPlane SurfaceA( &PointList );
  
  PointList[0] = P1;
  PointList[1] = P2; 
  PointList[2] = P3; 
  PointList[3] = P0;
  G4FPlane SurfaceB( &PointList ); 

  PointList[0] = P2;
  PointList[1] = P3; 
  PointList[2] = P0; 
  PointList[3] = P1;
  G4FPlane SurfaceC( &PointList ); 

  PointList[0] = P3;
  PointList[1] = P0; 
  PointList[2] = P1; 
  PointList[3] = P2;
  G4FPlane SurfaceD( &PointList );

  G4Plane plan1 = SurfaceA.GetPplane();
  G4Plane plan2 = SurfaceB.GetPplane();
  G4Plane plan3 = SurfaceC.GetPplane();
  G4Plane plan4 = SurfaceD.GetPplane();

  G4cout<<"\n\n Plan1 : a= "<<plan1.a
	<<" b= "<<plan1.b
	<<" c= "<<plan1.c
	<<" d= "<<plan1.d   ;
  G4cout<<"\n Plan2 : a= "<<plan2.a
	<<" b= "<<plan2.b
	<<" c= "<<plan2.c
	<<" d= "<<plan2.d   ;
  G4cout<<"\n Plan3 : a= "<<plan3.a
	<<" b= "<<plan3.b
	<<" c= "<<plan3.c
	<<" d= "<<plan3.d   ;
  G4cout<<"\n Plan4 : a= "<<plan4.a
	<<" b= "<<plan4.b
	<<" c= "<<plan4.c
	<<" d= "<<plan4.d   ;

  G4double d1 =  SurfaceA.ClosestDistanceToPoint(Pout);
  G4double d2 =  SurfaceB.ClosestDistanceToPoint(Pout);
  G4double d3 =  SurfaceC.ClosestDistanceToPoint(Pout);
  G4double d4 =  SurfaceD.ClosestDistanceToPoint(Pout);

  G4cout<<"\n\n Distance 1 ="<<d1;
  G4cout<<"\n Distance 2 ="<<d2;
  G4cout<<"\n Distance 3 ="<<d3;
  G4cout<<"\n Distance 4 ="<<d4;



//////////////////////////////////////////////////////////////////////
//
// Test for the function EvaluateIntersection
// 
//

  G4cout<<"\n\n//////////////////////////////////////////////////////";

  G4Point3D  Pdep1 (0, 0, 2);
  G4Vector3D DirRay(1 ,1, 0);
  G4Ray      Rayref(Pdep1, DirRay);

  G4cout<<"\n\nPdep of the ray :"
	<<"\n   x="<<Rayref.GetStart().x()
	<<"\n   y="<<Rayref.GetStart().y()
	<<"\n   z="<<Rayref.GetStart().z();

  G4cout<<"\n\nDirection of the ray :"
	<<"\n   x="<<Rayref.GetDir().x()
	<<"\n   y="<<Rayref.GetDir().y()
	<<"\n   z="<<Rayref.GetDir().z();

  int intersec = SurfaceVec1.Intersect(Rayref);

  if(intersec)
  {
    G4cout<<"\n\nIntersection founded at point :"
	  <<"\n   x="<<SurfaceVec1.hitpoint.x()
	  <<"\n   y="<<SurfaceVec1.hitpoint.y()
	  <<"\n   z="<<SurfaceVec1.hitpoint.z();

    if (  ( SurfaceVec1.hitpoint.x()*SurfaceVec1.GetPplane().a +
	    SurfaceVec1.hitpoint.y()*SurfaceVec1.GetPplane().b +
	    SurfaceVec1.hitpoint.z()*SurfaceVec1.GetPplane().c   < 
	    SurfaceVec1.GetPplane().d + kCarTolerance              ) &&
	  ( SurfaceVec1.hitpoint.x()*SurfaceVec1.GetPplane().a +
	    SurfaceVec1.hitpoint.y()*SurfaceVec1.GetPplane().b +
	    SurfaceVec1.hitpoint.z()*SurfaceVec1.GetPplane().c   > 
	    SurfaceVec1.GetPplane().d - kCarTolerance              )    )
       G4cout<<"\n\nPlain contain the hit point";
    else
      G4cout<<"\n\nPlain do not contain the hit point";

    G4cout<<"\n\nSquared distance from the Pdep to the hit point  :"
	  <<"\n   distance="<<SurfaceVec1.GetDistance();
  }
  else
    G4cout<<"\n\nNo Intersection"
	  <<"\n   distance="<<SurfaceVec1.GetDistance()<<G4endl;

  G4cout<<G4endl;
  return EXIT_SUCCESS;
}
