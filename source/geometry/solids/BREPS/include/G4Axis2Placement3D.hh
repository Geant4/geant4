#ifndef __G4Placement3D_h
#define __G4Placement3D_h 1

#include "G4Point3D.hh"
#include "G4Vector3D.hh"
#include "G4Transform3D.hh"
#include "G4PointRat.hh"
#include "G4Ray.hh"

class G4Axis2Placement3D 
{
public:
  G4Axis2Placement3D();
  ~G4Axis2Placement3D();
  
  G4Axis2Placement3D(const G4Axis2Placement3D& place);
  
  
  //inline void Project (G4ThreeVec& Coord, const G4ThreeVec& Pt2, 
  //                     const G4Plane& Pl1, const G4Plane& Pl2)
  // {
  //   Coord.X(Pt2.X()*Pl1.a + Pt2.Y()*Pl1.b + Pt2.Z()*Pl1.c - Pl1.d);
  //   Coord.Y(Pt2.X()*Pl2.a + Pt2.Y()*Pl2.b + Pt2.Z()*Pl2.c - Pl2.d);
  //   Coord.Z(0);
  // }

  // Get/Set for geometric data
  void Init( const G4Vector3D& refDirection0 , 
	     const G4Vector3D& axis0         ,
	     const G4Point3D&  location0      );

  G4Axis2Placement3D( const G4Vector3D& refDirection0 ,
		      const G4Vector3D& axis0         , 
		      const G4Point3D&  location0      );

  G4Point3D  GetLocation() const;
  G4Vector3D GetAxis() const;
  G4Vector3D GetRefDirection() const;

  // placement coordinate axes
  G4Vector3D GetPX() const;
  G4Vector3D GetPY() const;
  G4Vector3D GetPZ() const;

  // transformation from/to the placement coordinate system
  const G4Transform3D& GetToPlacementCoordinates() const;
  const G4Transform3D& GetFromPlacementCoordinates() const; 
  
  virtual G4bool operator==(const G4Axis2Placement3D& other) const 
  {
    return (this==&other) ? true : false; 
  }


private:

  // geometric data
  G4Point3D location;
  G4Vector3D axis;
  G4Vector3D refDirection;  

  // placement coordinate axes
  G4Vector3D pX, pY, pZ;

  G4Transform3D toPlacementCoordinates;
  G4Transform3D fromPlacementCoordinates;

};

#include "G4Axis2Placement3D.icc"

#endif




