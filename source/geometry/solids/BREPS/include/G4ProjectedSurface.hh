// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ProjectedSurface.hh,v 1.1 1999-01-07 16:07:34 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4ProjectedSurface_h
#define __G4ProjectedSurface_h 1

#include "G4BezierSurface.hh"

class G4ProjectedSurface  : public G4Surface
{
  friend class G4BSplineSurface;
  
  friend void CopySurface(G4ProjectedSurface& proj);
  
public:
  //Default constructor
  G4ProjectedSurface();
  ~G4ProjectedSurface();
	
  // Copy-constructor 
  G4ProjectedSurface(const G4ProjectedSurface &tmp);
	
  // Test variables
  static int Splits;
	
  void CalcBBox();
  G4ControlPoints *ctl_points;
  
  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt)const
  {
    return G4Vector3D(0,0,0);
  }
        
  
private:
  short dir;
  G4KnotVector *u_knots;
  G4KnotVector *v_knots;
  
  void CopySurface();
	
  void ConvertToBezier ( G4SurfaceList&, G4SurfaceList&);
	
  inline int GetOrder(int direction)
  {
    return order[direction];
  }
	
  inline void PutOrder(int direction, int value)		
  {
    order[direction]=value;
  }
  
  G4SurfaceList* projected_list;
  G4SurfaceList* bezier_list;   
	
  int order[2];
  G4KnotVector *new_knots;
  int ord;
  int lower,upper;
  
  G4OsloMatrix* oslo_m;
  G4Point3D vmin;
  G4Point3D vmax;
  void SplitNURBSurface();
  int  CheckBezier();
	
  void CalcOsloMatrix();
  void MapSurface(G4ProjectedSurface* srf);
  inline int Amax(int i, int j) 
  {
    return( (i) > (j) ? (i) : (j) );
  }
  
  inline int Amin(int i, int j) 
  {
    return( (i) < (j) ? (i) : (j) );
  }
  
  inline int AhIndex(int j,int t, int iorder)
  {
    return(( (j) * ((j)+1)/2) + (t) - ((iorder-1) - (j)));
  }
};	

#endif

