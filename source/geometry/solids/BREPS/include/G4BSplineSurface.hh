// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BSplineSurface.hh,v 1.2 1999-05-21 18:40:16 japost Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __BSPLINESURFACE_H
#define __BSPLINESURFACE_H


#include "G4Point3D.hh"
#include "G4PointRat.hh"
#include "G4Surface.hh"
#include "G4ProjectedSurface.hh"



//#ifdef WIN32
//#  include "G4ios.hh"
//#else
//#  include <stream.h>
//#endif



class G4BSplineSurface : public G4Surface
{
public:

  G4BSplineSurface();
  G4BSplineSurface(char*, G4Ray&);
  G4BSplineSurface(const  G4BSplineSurface &tmp);
  G4BSplineSurface(G4int, G4int, G4KnotVector&,  G4KnotVector&, 
		   G4ControlPoints&);    
  ~G4BSplineSurface();

  int Intersect(const G4Ray&);
  void CalcBBox();
  
  G4double GetUHit()  { return Hit->u; }  
  G4double GetVHit()  { return Hit->v; } 
    	 
  inline int MyType()const {return 2;}
  
  G4double ClosestDistanceToPoint(const G4Point3D&);

  inline void Reset()
  {
    active=1;
    bezier_list.EmptyList();
    projected_list.EmptyList();
    Intersected=0;
    distance = kInfinity;
  }

  // get for controlpoints
  G4int     GetRows()                         { return ctl_points->GetRows(); }
  G4int     GetCols()                         { return ctl_points->GetCols(); }
  G4Point3D GetControlPoint(G4int a, G4int b) { return ctl_points->Get3D(a,b);}
    
    
private:
  
  G4SurfaceList bezier_list;
  G4SurfaceList projected_list;
  short dir;
  int order[2];
  G4KnotVector *u_knots;
  G4KnotVector *v_knots;
  G4KnotVector *tmp_knots;
  G4ControlPoints *ctl_points;
  G4UVHit* Hit;
  G4UVHit* first_hit;
  int ord;
  int k_index;
  G4double param;
  int Rational;
  
  void FindIntersections(const G4Ray&);

  inline int GetOrder(int direction)             { return order[direction]; }
  inline void PutOrder(int direction, int value) { order[direction]=value;  }

  void AddHit(G4double u, G4double v);
  void ProjectNURBSurfaceTo2D( const G4Plane& ,const G4Plane&,
			       G4ProjectedSurface*);

  G4ProjectedSurface* CopyToProjectedSurface(const G4Ray&);
  G4Point3D  FinalIntersection();

  // L. Broglia
  // Because  G4BSplineSurface::Evaluate hides the virtual function 
  // G4Surface::Evaluate(const G4Ray&), I modified the function name
  // G4Point3D  Evaluate();  
  G4Point3D  BSEvaluate();

  G4PointRat& InternalEvalCrv(int i, G4ControlPoints *crv);
  
  G4Point3D   Evaluation(const G4Ray&);

  G4Vector3D  SurfaceNormal(const G4Point3D& Pt)const
  {
    return G4Vector3D(0,0,0);
  }
  
}; 

#endif











