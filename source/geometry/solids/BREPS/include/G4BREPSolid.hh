#ifndef __SOLID_H
#define __SOLID_H
#include "G4VSolid.hh"	      
#include "G4VisExtent.hh"      
#include "G4Surface.hh"
#include "G4Axis2Placement3D.hh"
#include "G4PointRat.hh"	 
#include "G4BoundingBox3D.hh"	 

class STEPentity;
class InstMgr;
class G4Ray;


class G4BREPSolid : public G4VSolid
{

public:

  G4BREPSolid(const G4String name);
  G4BREPSolid(const G4String, G4Surface**, G4int);
  ~G4BREPSolid();

  virtual G4String GetEntityType() const {return "Closed_Shell";}  
  virtual void Initialize();
  G4int CreateSTEPData(); // not yet implemented 

  G4bool CalculateExtent(const EAxis              pAxis      ,
			 const G4VoxelLimits&     pVoxelLimit,
			 const G4AffineTransform& pTransform ,
			 G4double&                pMin       , 
			 G4double&                pMax        ) const;

  virtual EInside Inside(register const G4ThreeVector&) const;

  virtual G4ThreeVector SurfaceNormal(const G4ThreeVector&) const;

  virtual G4double DistanceToIn(const G4ThreeVector&) const;
  virtual G4double DistanceToIn(register const G4ThreeVector&,
				register const G4ThreeVector&) const;

  virtual G4double DistanceToOut(const G4ThreeVector&) const;
  virtual G4double DistanceToOut(register const G4ThreeVector&,
				 register const G4ThreeVector&,
				 const G4bool  calcNorm=false , 
				 G4bool        *validNorm=0   , 
				 G4ThreeVector *n=0             ) const;
 

  G4Point3D Scope(); // ???

  void DescribeYourselfTo (G4VGraphicsScene& scene) const;
  G4VisExtent   GetExtent        () const;
  G4Polyhedron* CreatePolyhedron () const;
  G4NURBS*      CreateNURBS      () const;

  G4int Intersect(register const G4Ray&)const;

  inline G4double IntersectionDistance()const{return intersectionDistance;}
  void IntersectionDistance(const G4double d)const 
  {
    ((G4BREPSolid*)this)->intersectionDistance=d;
  }   

  G4Surface* GetSurface(G4int nr)
  {
    return SurfaceVec[nr];
  }
  
  inline void Active(const G4int x)const 
  {
    ((G4BREPSolid*)this)->active=x;
  }
  
  inline G4int Active() const {return active;}
  
  virtual inline void Reset() const
  {
    ((G4BREPSolid*)this)->active=1;
    ((G4BREPSolid*)this)->intersectionDistance=kInfinity;
    ((G4BREPSolid*)this)->startInside=0; 
     
    for(register G4int a=0;a<nb_of_surfaces;a++)
      SurfaceVec[a]->Reset();

    ShortestDistance = kInfinity;
  }

  static G4int NumberOfSolids;
  static InstMgr InstanceList;

  G4double GetShortestDistance() const {return ShortestDistance;}

  G4int GetId() const {return Id;}

  void SetId(G4int id) {Id = id;}

  G4String GetName() const {return solidname;}

  void SetName(G4String name) {solidname = name;}

  G4int NumberOfFaces() const {return nb_of_surfaces;} 

  // Add by L. Broglia
  G4Axis2Placement3D* GetPlace() { return place; }
  G4BoundingBox3D*    GetBBox()  { return bbox;  } 
 
protected:

  G4bool  IsConvex();
  virtual void CalcBBoxes();
  void    CheckSurfaceNormals();
  void    RemoveHiddenFaces(register const G4Ray& G4Rayref, G4int)const;   
  void    TestSurfaceBBoxes(register const G4Ray&) const;

  inline G4int StartInside() const 
  {
    return startInside;
  }
  
  inline void StartInside(const G4int si) const 
  {
    ((G4BREPSolid*)this)->startInside=si;
  }  


private:

  G4int IsBox();
  G4int FinalEvaluation(register const G4Ray&, const G4int =0) const;


protected:
  G4Axis2Placement3D* place;
  static G4Ray        Track;
  static G4double     ShortestDistance;
  G4int               Box, Convex, AxisBox, PlaneSolid;
  G4BoundingBox3D*    bbox;   
  G4double            intersectionDistance;
  G4int               active;
  G4int               startInside;
  G4int               nb_of_surfaces;
  G4Point3D           intersection_point;
  G4Surface**         SurfaceVec;
  G4double            RealDist;
  G4String            solidname; 
  G4int               Id;
   

  void QuickSort( register G4Surface** SrfVec, 
		  register G4int left, register G4int right) const
  {
    register G4int i=left;
    register G4int j=right;
    register G4Surface* elem1;
    register G4Surface* elem2 = SrfVec[(left+right)/2];
    register G4double tmpdistance;
    do
    {
      tmpdistance = elem2->Distance();
      while ( SrfVec[i]->Distance() < tmpdistance  && i < right ) i++;
      while (tmpdistance < SrfVec[j]->Distance()  && j > left ) j--;
      if(i<=j)
      {
	elem1 = SrfVec[i];
	SrfVec[i] = SrfVec[j];
	SrfVec[j] = elem1;
	i++;j--;
      }
     } while (i<=j);
      
     if( left < j  ) QuickSort(SrfVec,left, j );
     if( i < right ) QuickSort(SrfVec,i, right);    
    }

};

#endif









