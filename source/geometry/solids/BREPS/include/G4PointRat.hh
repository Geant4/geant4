// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PointRat.hh,v 1.3 1999-11-08 09:50:31 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// Modif 8 oct 98 : A.Floquet
//      G4PointRat datas are made of
// 	  . a point 3D
//	  . a additional value : the scale factor which is set to 1 by default
//

#ifndef __G4POINT_RAT
#define __G4POINT_RAT

#include "geomdefs.hh"              // For kInfinity

#include "G4Point3D.hh"
#include "G4Plane3D.hh"

// L. Broglia
// Before included in G4Point.hh
#include "STEPaggregate.h"
#include "G4Plane.hh"
#include "G4UVHit.hh"
#define SQRT_SMALL_FASTF 1.0e-18
#define SMALL                   SQRT_SMALL_FASTF       
#define ROW 0
#define COL 1
// const G4double  INFINITY = 9.0e+99;
const G4Point3D PINFINITY(kInfinity, kInfinity, kInfinity);

class G4Plane;

class G4PointRat
{
public:
  G4PointRat();

  G4PointRat(const G4Point3D&);

  ~G4PointRat();

  void CopyRationalValue(const RealNode& Rnode);

  int GetType(void)const { return 4; } // This function should be removed
				       // if calls to this are also removed
  void operator=(const G4Point3D&); 

  void operator=(const G4PointRat&);

  inline G4double x() const {return pt3d.x();}

  inline void setX (const G4double Value) { pt3d.setX ( Value );}

  inline G4double y() const {return pt3d.y();}

  inline void setY (const G4double Value) { pt3d.setY ( Value );}

  inline G4double z() const {return pt3d.z();}

  inline void setZ (const G4double Value) { pt3d.setZ ( Value );}

  inline G4double w() const {return s;}

  inline void setW(const G4double Value) {s=Value;}

  inline G4Point3D pt() const { return pt3d; }


  inline G4double PlaneDistance(const G4Plane3D& Pl)
    {
      return ((Pl.a()*pt3d.x() + Pl.b()*pt3d.y() + Pl.c()*pt3d.z()) - Pl.d());
    }
  
private:
   G4Point3D pt3d;
   G4double  s;


public :

// L. Broglia
/*
  inline G4Point3D Min(const G4Point3D& p)
  {
    if(pt3d.x() < p.x()) pt3d.setX(p.x()); 
    if(pt3d.y() < p.y()) pt3d.setY(p.y()); 
    if(pt3d.z() < p.z()) pt3d.setZ(p.z()); 
  }     

  inline G4Point3D Max(const G4Point3D& p)     
  {
    if(pt3d.x() > p.x()) pt3d.setX(p.x()); 
    if(pt3d.y() > p.y()) pt3d.setY(p.y()); 
    if(pt3d.z() > p.z()) pt3d.setZ(p.z()); 
  }
*/


/*
  inline G4Point3D Min(const G4Vector3D& v)
  {
    if(pt3d.x() < v.x()) pt3d.setX(v.x()); 
    if(pt3d.y() < v.y()) pt3d.setY(v.y()); 
    if(pt3d.z() < v.z()) pt3d.setZ(v.z()); 
  }     

  inline G4Point3D Max(const G4Vector3D& v)     
  {
    if(pt3d.x() > v.x()) pt3d.setX(v.x()); 
    if(pt3d.y() > v.y()) pt3d.setY(v.y()); 
    if(pt3d.z() > v.z()) pt3d.setZ(v.z()); 
  }
*/
};

#endif

