//////////////////////////////////////////////////////////////////////////
// $Id: CurveTest.cc,v 1.3 2000-08-28 08:58:03 gcosmo Exp $
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
// L. Broglia
// October 26, 1998
//
// Test the basic functionnality of the different curves :
//   - line
//   - circular curve
//   - ellipse
//   - parabola
//   - hyperbola
//
// Test the functionality of Axis2Placement3D and Project
//
// Use the function included in TestFunction.hh
//


#include "CurveTestFunction.hh"


G4int main()
{
  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4Axis2Placement3D";
  G4cout<<"\n\n";

  G4Axis2Placement3D p;
  p.Init( G4Vector3D(-1, 0, 0), 
	  G4Vector3D( 1, 2, 3),
	  G4Point3D ( 1, 1, 1) );
  TestPlacement(&p);
  
  G4Axis2Placement3D p2;
  p2.Init( G4Vector3D( 1, 0, 0), 
	   G4Vector3D(-1, 2, 3),
	   G4Point3D (-1, 1, 1) );
  TestPlacement(&p2);
  


  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4Line";
  G4cout<<"\n\n";
  
  G4Line l;
  l.Init(G4Point3D(2, 3, 4), G4Vector3D(-1, -1, -1));

  l.SetBounds(-10, +10);
  G4cout << "G4Line " << l.GetPnt() << " " << l.GetDir() << G4endl; 
  TestCurve(&l);  



  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4CircularCurve";
  G4cout<<"\n\n";
  
  G4CircularCurve c;
  c.Init(p, 5);

  c.SetBounds(0, 0);
  G4cout << "G4CircularCurve with radius=" << c.GetRadius() << G4endl;
  TestCurve(&c);
  
  c.SetBounds(pi/4, 3*pi/4);
  G4cout << "G4CircularCurve again" << G4endl;
  TestCurve(&c);
  
  

  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4Ellipse";
  G4cout<<"\n\n";

  G4Ellipse e;
  e.Init(p, 2, 10);

  e.SetBounds(0, 0);
  G4cout << "G4Ellipse " << e.GetSemiAxis1()
	 << " " << e.GetSemiAxis2() << G4endl;
  TestCurve(&e);
 
  e.SetBounds(100, 100+8*pi);
  G4cout << "G4Ellipse again" << G4endl;
  TestCurve(&e);
  
  e.SetBounds(3.21, 4.1);
  G4cout << "G4Ellipse again" << G4endl;
  TestCurve(&e);    



  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4Parabola";
  G4cout<<"\n\n";

  G4Parabola par;
  par.Init(p, 6);

  par.SetBounds(-100, 100);
  G4cout << "G4Parabola " << par.GetFocalDist() << G4endl;
  TestCurve(&par);
  
  par.SetBounds(-1, 10);
  G4cout << "G4Parabola again" << G4endl;
  TestCurve(&par);
  
  par.SetBounds(-6.2, -5.1);
  G4cout << "G4Parabola again" << G4endl;
  TestCurve(&par);
  


  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of G4Hyperbola";
  G4cout<<"\n\n";
  
  G4Hyperbola h;
  h.Init(p, 5, 3);

  h.SetBounds(-2, 2);
  G4cout << "G4Hyperbola " << h.GetSemiAxis()
	 << " " << h.GetSemiImagAxis() << G4endl;
  TestCurve(&h);
  
  h.SetBounds(1, 2);
  G4cout << "G4Hyperbola again" << G4endl;
  TestCurve(&h);

  h.SetBounds(-1, 0);
  G4cout << "G4Hyperbola again" << G4endl;
  TestCurve(&h);



  G4cout<<"\n\n//////////////////////////////////////////////////////////////";
  G4cout<<"\n\n Test the basic functionality of Projection";
  G4cout<<"\n\n";

 TestProject(  &e, G4Rotate3D( pi/2, G4Vector3D(1, 1, 1) )  ); 



 return EXIT_SUCCESS;
}
