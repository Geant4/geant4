// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ToroidalSurface.hh,v 1.2 1999-12-15 14:49:58 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef __G4TOROIDALSURAFCE
#define __G4TOROIDALSURAFCE

#include "G4FPlane.hh"
#include "G4OsloMatrix.hh"

class G4ToroidalSurface:public G4Surface
{
 public:

  G4ToroidalSurface();

  G4ToroidalSurface(const G4Vector3D&, 
		    const G4Vector3D&,  
		    const G4Vector3D&,  
		    const G4double, 
		    const G4double);    

  ~G4ToroidalSurface();

  G4String GetEntityType(){return G4String("Toroidal_Surface");}

  int Intersect(const G4Ray&);
  
  void CalcBBox();
  
  inline G4Vector3D  	GetDirection(){return Placement.GetRefDirection();}
  
  inline G4Vector3D 	GetAxis()     {return Placement.GetAxis();}
  
  inline G4Point3D	GetLocation() {return Placement.GetLocation();}
  
  inline G4double   	GetMinRadius(){return MinRadius;} 
  
  inline G4double   	GetMaxRadius(){return MaxRadius;} 
  
  G4double ClosestDistanceToPoint(const G4Point3D&);
  
  G4Vector3D SurfaceNormal(const G4Point3D& Pt)const 
                                         {return G4Vector3D(0,0,0);}


  inline void 	MultiplyPointByMatrix(G4Point3D& Base)
  {
    Base.setX((Base.x() * TransMatrix->get(0,0)) + 
	      (Base.y() * TransMatrix->get(1,0)) + 
	      (Base.z() * TransMatrix->get(2,0)));
    Base.setY((Base.x() * TransMatrix->get(0,1)) + 
	      (Base.y() * TransMatrix->get(1,1)) + 
	      (Base.z() * TransMatrix->get(2,1)));
    Base.setZ((Base.x() * TransMatrix->get(0,2)) + 
	      (Base.y() * TransMatrix->get(1,2)) + 
	      (Base.z() * TransMatrix->get(2,2)));
  }

  inline void	MultiplyVectorByMatrix(G4Vector3D& DCos)
  {
    G4double w;
    DCos.setX((DCos.x() * TransMatrix->get(0,0)) + 
	      (DCos.y() * TransMatrix->get(1,0)) + 
	      (DCos.z() * TransMatrix->get(2,0)) + TransMatrix->get(3,0));

    DCos.setY((DCos.x() * TransMatrix->get(0,1)) + 
	      (DCos.y() * TransMatrix->get(1,1)) + 
	      (DCos.z() * TransMatrix->get(2,1)) + TransMatrix->get(3,1));

    DCos.setY((DCos.x() * TransMatrix->get(0,2)) + 
	      (DCos.y() * TransMatrix->get(1,2)) + 
	      (DCos.z() * TransMatrix->get(2,2)) + TransMatrix->get(3,2));
  
    w      = ((DCos.x() * TransMatrix->get(0,3)) + 
	      (DCos.y() * TransMatrix->get(1,3)) + 
	      (DCos.z() * TransMatrix->get(2,3)) + TransMatrix->get(3,3));
    
    if (w != 0.0) 
    {
      DCos.setX(DCos.x() / w);  
      DCos.setY(DCos.y() / w); 
      DCos.setZ(DCos.z() / w);
    }
  }
  

private:

  G4Axis2Placement3D Placement;
  G4double           MinRadius;
  G4double	     MaxRadius;
  Matrix*            TransMatrix;   // transformation matrix  
  G4Point3D          hitpoint;
  const G4double     EQN_EPS;

  int SolveQuartic(G4double c[], G4double s[]);	

  inline int IsZero(G4double x)
  {
    if((x) > -EQN_EPS && (x) < EQN_EPS) 
      return 1;
    else return 0;
  }

  int SolveCubic(G4double c[], G4double s[]);
  
  int SolveQuadric(G4double c[], G4double s[]);

};


#endif












