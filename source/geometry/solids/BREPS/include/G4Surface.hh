#ifndef __surface_h
#define __surface_h 1

#include "geomdefs.hh"
#include "G4CurveVector.hh"
#include "G4PointRat.hh"
#include "G4Ray.hh"
#include "G4BoundingBox3D.hh"
#include "G4STEPEntity.hh"
#include "G4SurfaceBoundary.hh"

// This is the combined G4Surface class
class G4Surface: public G4STEPEntity
{
public:

  G4Surface();
  virtual ~G4Surface();

  // sets the boundaries of the surface.
  // The curves in the CurveVector must be non-intersecting
  // closed curves.
  void SetBoundaries(G4CurveVector*);
  // It calls InitBounded -- empty by default

protected:

  virtual void InitBounded() { }

public:

  // type information, needed for STEP output (see STEPinterface)
  virtual G4String GetEntityType(){return G4String("Surface");}

  // The origin should move to the derived classes
  int operator==( const G4Surface& s ) { return origin == s.origin; }

  // such a function is needed
  // (see G4VSolid::DistanceToIn(const G4ThreeVector&) )
  // but the G4surface implementation is useless.
  // Overriding functions don't take the surface
  // boundary into account.
  virtual G4double HowNear( const G4Vector3D& x ) const; 

  //virtual G4double distanceAlongRay( int which_way, const G4Ray* ry,
  //				       G4Vector3D& p ) const;

  // unnecessary -- origin should move to descendants
  G4Vector3D GetOrigin() const { return origin; }

  // Gerep members
  // bad function names -- use Set and Get
  // ??
  inline G4double Distance()    { return distance; }
  inline void Distance(const G4double Dist)  { distance=Dist; }    

  // a boolean flag, not used by the surfaces themselves
  virtual inline int  Active(){return active;}
  virtual inline void Active(const int act){active=act;}

  // Isn't this the same as HowNear? (This one is used by G4BREPSolid.)
  virtual G4double ClosestDistanceToPoint(const G4Point3D&);

  // uhit and vhit are never set.
  // Only BSplineSurface overrides.
  // There is a G4UVHit class.
  virtual G4double GetUHit() { return uhit; }
  virtual G4double GetVHit() { return vhit; }  

  // Intersection with a ray. the result is put into
  // some data members.
  virtual int Intersect(const G4Ray&);

  // Surface normal calculation.
  virtual G4Vector3D Normal( const G4Vector3D& p ) const;

  // Bounding box calculation.
  virtual void CalcBBox();

  // For NURBS, there is a two pass intersection algorithm.
  // Sometimes, the result of the cheap one tells us 
  // that execution of the expensive one is not necessary.
  // Evaluation (Evaluate?) is one of them.
  // better names wanted!
  virtual G4Point3D Evaluation(const G4Ray& G4Rayref);
  virtual int Evaluate(register const G4Ray& Rayref);  

  // There is Active(int) instead.
  virtual inline void Deactivate(){active=0;}
  
  // Distance(kInfinity); bbox->SetDistance(kInfinity);};

  virtual inline void Reset(){Intersected=0;active = 1; distance = kInfinity;};

  // one function for type info (GetEntityType) should be enough
  virtual char *Name() const { return (char*)("G4Surface"); }
  virtual int MyType() const { return Type; }  

  // To be replaced by a CLHEP vector operation
  inline static void Project (G4double& Coord, const G4Point3D& Pt2, 
			      const G4Plane& Pl1                     )
  {
    Coord = Pt2.x()*Pl1.a + Pt2.y()*Pl1.b + Pt2.z()*Pl1.c - Pl1.d;
  }

  // Used by BREPSolid. Thus it's probably needed.
  virtual void Project(){}

  // Only in G4FPlane. Should be private to that class?
  virtual void CalcNormal(){}  

  // Only in G4FPlane. BREPSolid::IsConvex uses it.
  // But who uses BREPSolid::IsConvex?
  // Thus: probably not needed. But knowing
  // if the surface is convex could be used for optimization. 
  virtual int IsConvex(){return -1;}

  // Only in G4FPlane, but G4BREPSolid uses them.
  virtual int GetConvex(){return 0;}

  virtual int GetNumberOfPoints(){return 0;}

  virtual const G4Point3D& GetPoint(const int Count)
  {
    const G4Point3D* tmp= new G4Point3D(0,0,0);
    return *tmp;
  }

  // L. Broglia
  void   SetSameSense(G4int sameSense0) { sameSense = sameSense0; }
  G4int  GetSameSense()                 { return sameSense      ; }
  
  G4BoundingBox3D* GetBBox() { return bbox; }

  // there is Normal as well -- so what do these do?
  virtual G4Ray* Norm(){return (G4Ray*)0;}
  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt) const =0;  

  // should be at least protected, but BREPSolid uses these data members.
  // So why not a Get function?


public:

  G4BoundingBox3D* bbox;
  G4Point3D closest_hit;

protected:

  // The boundaries of the surface.
  G4SurfaceBoundary surfaceBoundary;

  // BSplineSurface anf FPlane sets it, no one gets it
  int Intersected;

  // see Get... members
  G4Vector3D origin;	// origin of Surface
  int Type;
  int AdvancedFace;
  int active;
  G4double distance;
  G4double uhit,vhit;

  // L. Broglia
  G4int sameSense;

protected:

  // Maybe kInfinity instead?
  const G4double FLT_MAXX;

  // Maybe kCarTolerance instead?
  const G4double FLT_EPSILO;

  // temporary solution so that G4SurfaceList sees this member
  // but G4SurfaceList should go.


public:
  G4Surface* next;

};

#endif





