#ifndef __G4BoundingBox3D_h
#define __G4BoundingBox3D_h 1

#include "G4Ray.hh"
#include "G4Point3D.hh"
#include "G4Vector3D.hh"

class G4BoundingBox3D
{
public:    

  G4BoundingBox3D();

  G4BoundingBox3D(const G4Point3D&);
  G4BoundingBox3D(const G4Point3D&, const G4Point3D&);    
  ~G4BoundingBox3D();
    
  void Init(const G4Point3D&);
  void Init(const G4Point3D&, const G4Point3D&);
  void Extend(const G4Point3D&);

  G4Point3D GetBoxMin() const;
  G4Point3D GetBoxMax() const;

  G4double GetDistance() const;
  void SetDistance(G4double distance0);

  int GetTestResult() const;
  int Test(const G4Ray&);

  // this function return 1 if the point is inside and on the bbox,
  // 0 if the point is outside the bbox
  G4int Inside(const G4Point3D&);

  static const G4BoundingBox3D space;

    
private:

  G4Point3D box_min;
  G4Point3D box_max;
  G4double distance;

  int test_result;

  G4Point3D MiddlePoint;
  G4Vector3D GeantBox;    

  int BoxIntersect(const G4Point3D&, 
		   const G4Point3D&, 
		   const G4Vector3D&) const;

  G4double DistanceToIn(const G4Point3D&,
			const G4Vector3D&) const;			  
};


#include "G4BoundingBox3D.icc"

#endif





