// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4FPlane.hh,v 1.3 1999-01-27 16:13:09 broglia Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// L. Broglia
// 
// A G4FPlane is a plane created by 3 points or by an origin, an axis and 
// a direction. The plane created is a G4Plane, where his coefficient a, b, 
// c and d are stored. Be carreful that the equation of the plane is :
//                ax + by + cz = d
// 
// This class contain 2 intersection functions :
//       - closest intersection 
//       - intersection by a ray
//
//

#ifndef __PLANESURFACE_H
#define __PLANESURFACE_H


#include "G4Axis2Placement3D.hh"
#include "G4Plane.hh"
#include "G4Surface.hh"


class G4FPlane:public G4Surface
{

public:

  // Default constructor - destructor
  G4FPlane();
  ~G4FPlane() { delete NormalX; }

  // Normal constructor
  G4FPlane( const G4Vector3D& direction, 
	    const G4Vector3D& axis     ,	   
	    const G4Point3D&  Pt0       );
  
  // Constructor used by G4BREPSolidBox and G4BREPSolidPolyhedra
  G4FPlane(const G4Point3DVector* pVec, 
	   const G4Point3DVector* iVec= 0);

  // hit point of the ray on the surface
  G4Point3D hitpoint;

  // calculate the intersection of the plane and a ray
  int Intersect(const G4Ray& G4Rayref);
  //int Evaluate(const G4Ray& ray) { return Intersect(ray); }
  
  // Calculate bounding box
  void CalcBBox();

  // Calculate the projection of the plane
  void Project();
  
  // return the type, used in G4BREPSolid
  inline int MyType()const { return 1; }
  
  // return the convexity or not
  int GetConvex() { return Convex; }

  // is convex ?
  int IsConvex(); 
  
  // deactive, used in G4Surface
  inline void Deactivate() { active=0; }

  // get the number of the points on the surface boundary
  inline int GetNumberOfPoints() 
  {
    return (surfaceBoundary.GetNumberOfPoints());
  }

  // get the location point
  G4Point3D GetSrfPoint() { return pplace.GetLocation(); }
  
  // get a surface boundary point
  inline const G4Point3D& GetPoint(const int Count)
  {
    return surfaceBoundary.GetPoint(Count);
  }  
  
  void CalcNormal();
  
  // return the normal, used in BREPSolid
  G4Ray* Norm() { return NormalX; }

  G4Vector3D SurfaceNormal(const G4Point3D& Pt)const 
  {
    if( !sameSense )
      return -NormalX->GetDir();
    else
      return NormalX->GetDir();
  }
  
  virtual char *Name() const { return "G4FPlane"; }

  G4double ClosestDistanceToPoint(const G4Point3D& Pt);

  // L. Broglia : create this Surface function 
  virtual G4double HowNear( const G4Vector3D& x ) const ;

  inline G4Axis2Placement3D GetPplace() const { return  pplace; }

  inline G4Plane GetPplane() const { return  Pl; }

private:

  G4Axis2Placement3D pplace;
  G4Plane            Pl;
  G4Ray              *NormalX;
  int                Convex;
  G4SurfaceBoundary* projectedBoundary;

  inline int Sign(const G4double a)
  {
    register int i=1; 
    if(a<0) 
      i= -1;
    
    return i;
  } 
  
  
protected:

  // P. Urban
  virtual void InitBounded();
  
};

#endif




