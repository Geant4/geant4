// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SphericalSurface.hh,v 1.2 1999-12-15 14:49:57 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __G4SpheShell_H
#define __G4SpheShell_H

#include "G4Surface.hh"
#include "G4ThreeMat.hh"
// #include "G4Vector3D.hh"	already included in G4ThreeMat


class G4SphericalSurface: public G4Surface
{

protected:
  G4Vector3D x_axis;  // direction (unit vector) of axis of G4SphericalSurface
                      // which defines azimuthal angle of zero

  G4Vector3D z_axis;  // direction (unit vector) of axis of G4SphericalSurface
                      // which defines polar angle of zero
	
  G4double radius;    // radius of G4SphericalSurface

  G4double phi_1;     // lower azimuthal angle limit of G4SphericalSurface
		      // (in radians).  Allowed range 0 <= phi_1 < 2*PI
  
  G4double phi_2;     // upper azimuthal angle limit of G4SphericalSurface
		      // (in radians).  Allowed range 
		      // phi_1 < phi_2 <= phi_1 + 2*PI

  G4double theta_1;   // lower polar angle limit of G4SphericalSurface
                      // (in radians).  Allowed range 0 <= theta_1 < PI
	
  G4double theta_2;   // upper polar angle limit of G4SphericalSurface
		      // (in radians).  Allowed range
		      // theta_1 < theta_2 <= theta_1 + PI

public:

  G4SphericalSurface();
  G4SphericalSurface( const G4Vector3D& o, 
		      const G4Vector3D& xhat, const G4Vector3D& zhat,
		      G4double r, 
		      G4double ph1, G4double ph2,
		      G4double th1, G4double th2                       ); 
  ~G4SphericalSurface() {}

  G4String GetEntityType() { return G4String("Spherical_Surface"); }

  // G4SphericalSurface( const G4SphericalSurface& s ): G4Surface( s.origin )
  //		{ x_axis = s.x_axis;
  //		  z_axis = s.z_axis;
  //		  radius = s.radius;
  //		  phi_1 = s.phi_1;
  //		  phi_2 = s.phi_2;
  //		  theta_1 = s.theta_1;
  //		  theta_2 = s.theta_2; }			                               
  

  int Intersect(const G4Ray&);
  
  void CalcBBox();
  
  inline void Comp( G4Vector3D& v, G4Point3D& min , G4Point3D& max)
  {
    // Compares the x,y and z values of v and min
    // / v and max. min/max-values are replaced if 
    // greater/smaller than v-values.
    
    if(v.x() > max.x()) max.setX(v.x());
    if(v.y() > max.y()) max.setY(v.y());
    if(v.z() > max.z()) max.setZ(v.z());

    if(v.x() < min.x()) min.setX(v.x());
    if(v.y() < min.y()) min.setY(v.y());
    if(v.z() < min.z()) min.setZ(v.z());
  }

  virtual char *NameOf() const { return "G4SphericalSurface"; }
  virtual void PrintOn( G4std::ostream& os = G4cout ) const;
  
  int operator==( const G4SphericalSurface& s )
  { return origin  == s.origin  &&  
      x_axis  == s.x_axis  &&
      z_axis  == s.z_axis  &&
      radius  == s.radius  && 
      phi_1   == s.phi_1   &&
      phi_2   == s.phi_2   &&
      theta_1 == s.theta_1 &&
      theta_2 == s.theta_2;   }   
  
  virtual G4double HowNear( const G4Vector3D& x ) const;
        
  //virtual G4double distanceAlongRay( int which_way, const G4Ray* ry,
  //	G4ThreeVec& p ) const;
  //	virtual G4double distanceAlongHelix( int which_way, const Helix* hx,
  //		G4ThreeVec& p ) const;
  //	virtual G4Vector3D Normal( const G4Point3D& p ) const;
  
  virtual G4Vector3D Normal( const G4Vector3D& p ) const;
  virtual G4Vector3D SurfaceNormal( const G4Point3D& p ) const;
	      
  virtual int Inside( const G4Vector3D& x ) const;
  virtual int WithinBoundary( const G4Vector3D& x ) const;

  virtual G4double Scale() const;
  virtual G4double Area() const; 

  virtual void resize( G4double r, G4double ph1, G4double ph2, 
		       G4double th1, G4double th2);

  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, G4ThreeMat& m, int inverse ); 
  //	virtual void rotate( G4double alpha, G4double beta, 
  //		G4double gamma, int inverse ); 
  //

  G4Vector3D GetXAxis() const { return x_axis; }
  G4Vector3D GetZAxis() const { return z_axis; }

  G4double   GetRadius() const { return radius; }

  G4double   GetPhi1() const  { return phi_1; }
  G4double   GetPhi2() const  { return phi_2; }

  G4double   GetTheta1() const { return theta_1; }
  G4double   GetTheta2() const { return theta_2; }
  
private:
//	virtual G4double gropeAlongHelix( const Helix* hx ) const;



//
//  Description of functions -----------------------------------------
//
//  default constructor
//----->G4SphericalSurface();
//
//  Normal constructor: first argument is the origin of the G4SphericalSurface
//		       second argument is the axis of the G4SphericalSurface
//				which defines azimuthal angle equals zero
//		        third argument is the axis of the G4SphericalSurface
//				which defines polar angle equals zero
//		       fourth argument is the radius of the G4SphericalSurface
//			fifth argument is the lower azimuthal angle limit of
//				the G4SphericalSurface
//			sixth argument is the upper azimuthal angle limit of
//				the G4SphericalSurface
//		      seventh argument is the lower polar angle limit of
//				the G4SphericalSurface
//			eigth argument is the upper polar angle limit of
//				the G4SphericalSurface
//----->G4SphericalSurface( const G4ThreeVec& o, const G4ThreeVec& xhat, 
//----->           const G4ThreeVec& zhat,
//----->           G4double r, G4double ph1, G4double ph2,
//----->	   G4double th1, G4double th2 ); 
//
//  destructor 
//----->virtual ~G4SphericalSurface() {}
//
//  copy constructor
//----->G4SphericalSurface( const G4SphericalSurface& s ): Surface( s.origin )
//----->	{ x_axis = s.X()_axis;
//----->	  z_axis = s.Z()_axis;
//----->	  radius = s.radius;
//----->	  phi_1 = s.phi_1;
//----->	  phi_2 = s.phi_2;
//----->	  theta_1 = s.theta_1;
//----->	  theta_2 = s.theta_2; }			                               
//
//  function to return class name
//----->virtual char *NameOf() const { return "G4SphericalSurface"; }
//
//  printing function
//----->virtual void PrintOn( G4std::ostream& os = G4cout ) const;
//
//  equality operator
//----->int operator==( const G4SphericalSurface& s )
//----->	{ return origin  == s.origin  &&  
//----->		 x_axis  == s.X()_axis  &&
//----->		 z_axis  == s.Z()_axis  &&
//----->		 radius  == s.radius  && 
//----->		 phi_1   == s.phi_1   &&
//----->		 phi_2   == s.phi_2   &&
//----->		 theta_1 == s.theta_1 &&
//----->		 theta_2 == s.theta_2;   }   
//
//  function which returns the distance from a point to a G4SphericalSurface
//		the (input) argument is the point x
//		the distance is positive if the point is Inside,
//		negative if it is outside
//----->virtual G4double HowNear( const G4ThreeVec& x ) const;
//
//  function which returns the distance along a Ray to enter or leave a
//  G4SphericalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Ray
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongRay( int which_way, const Ray* ry,
//----->	G4ThreeVec& p ) const;
//
//  function which returns the distance along a Helix to enter or leave a
//  G4SphericalSurface.  
//		the first (input) argument is +1 to leave or -1 to enter
//		the second (input) argument is a pointer to the Helix
//		the third (output) argument returns the intersection point
//----->virtual G4double distanceAlongHelix( int which_way, const Helix* hx,
//----->	G4ThreeVec& p ) const;
//
//  function which returns the Normal unit vector to a G4SphericalSurface at a point p
//  on (or nearly on) the G4SphericalSurface
//----->virtual G4ThreeVec Normal( const G4ThreeVec& p ) const;
//
//  function which returns true (1) if the point x is Inside the 
//  G4SphericalSurface, returns false (0) otherwise
//----->virtual int Inside( const G4ThreeVec& x ) const;
//
//  function which returns true (1) if the point x is within the boundary,
//  false (0) otherwise.
//----->virtual int WithinBoundary( const G4ThreeVec& x ) const;
//
//  function which returns the radius, unless it is zero, in which case it
//  returns 1.  Used for Scale-invariant tests of surface thickness.
//----->virtual G4double Scale() const;
//
//  function to calculate the Area of a G4SphericalSurface
//----->virtual G4double Area() const; 
//
//  function to resize the G4SphericalSurface to new radius and angle limits
//		         first argument is the radius of the G4SphericalSurface
//			second argument is the lower azimuthal angle limit of
//				the G4SphericalSurface
//			 third argument is the upper azimuthal angle limit of
//				the G4SphericalSurface
//			fourth argument is the lower polar angle limit of
//				the G4SphericalSurface
//			 fifth argument is the upper polar angle limit of
//				the G4SphericalSurface
//----->virtual void resize( G4double r, G4double ph1, G4double ph2, 
//----->			       G4double th1, G4double th2);
//
//  function to rotate the G4SphericalSurface (4 input arguments) 
//		 first about global x_axis by angle alpha,
//		second about global y-axis by angle beta,
//		 third about global z_axis by angle gamma
//		the angles are assumed to be given in radians
//		the fourth (output) argument gives the calculated rotation
//			matrix
//		the fifth (input) argument is an integer flag which if
//			non-zero reverses the order of the rotations
//----->virtual void rotate( G4double alpha, G4double beta, 
//----->	G4double gamma, G4ThreeMat& m, int inverse ); 
//
//  function to rotate the G4SphericalSurface (4 input arguments) 
//		 first about global x_axis by angle alpha,
//		second about global y-axis by angle beta,
//		 third about global z_axis by angle gamma
//		the angles are assumed to be given in radians
//		the fourth (input) argument is an integer flag which if
//			non-zero reverses the order of the rotations
//----->virtual void rotate( G4double alpha, G4double beta, 
//----->	G4double gamma, int inverse ); 
//
//  functions to return the axes, radius, and angles of the G4SphericalSurface
//----->direction GetXAxis() const { return x_axis; }
//----->direction GetZAxis() const { return z_axis; }
//----->G4double   GetRadius() const { return radius; }
//----->G4double   GetPhi1() const  { return phi_1; }
//----->G4double   GetPhi2() const  { return phi_2; }
//----->G4double   GetTheta1() const { return theta_1; }
//----->G4double   GetTheta2() const { return theta_2; }
//
//
//  Private function to use a crude technique to find the intersection
//  of a Helix with a G4SphericalSurface.  It returns the turning angle 
//  along the Helix at which the intersection occurs or -1.0 if no intersection
//  point is found.  The argument to the call is the pointer to the Helix.
//----->virtual G4double gropeAlongHelix( const Helix* hx ) const;
};

#endif
