// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4BezierSurface.hh,v 1.2 1999-12-15 14:49:56 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
#ifndef __BEZIERSURFACE_H
#define __BEZIERSURFACE_H

#include "G4Ray.hh"
#include "G4ControlPoints.hh"
#include "G4SurfaceList.hh"
#include "G4PointRat.hh"
#include "G4OsloMatrix.hh"
#include "G4KnotVector.hh"     


class G4ProjectedSurface;

class G4BezierSurface : public G4Surface
{
  friend class G4BSplineSurface;
  friend class G4ProjectedSurface; 
   
public:
  // Test variables 
  static int Clips;
  static int Splits;    
  
  G4BezierSurface();
  G4BezierSurface(const G4BezierSurface &tmp);
  ~G4BezierSurface();
  
  friend void CopySurface(G4BezierSurface& bez);
  static G4double Tolerance;
  
  inline G4Point3D AveragePoint()               { return average_pt;         };
  inline void SetAveragePoint(G4Point3D p)      { average_pt=p;              }

  inline G4double UAverage()                    { return average_u;          }
  inline G4double VAverage()                    { return average_v;          }

  inline void Dir(int d)                        { dir=d;                     }
  inline void ChangeDir()                       { dir=!dir;                  }

  inline G4double SMin()                        {return smin;                }
  inline G4double SMax()                        {return smax;                }

  inline int GetOrder(int direction)            { return order[direction];   }
  inline void PutOrder(int direction, int value){ order[direction]=value;    }

  inline G4double GetU()                        { return (u_min + u_max)/2.0;}
  inline G4double GetV()                        { return (v_min + v_max)/2.0;}

  void CalcBBox();

  // L. Broglia
  // Because G4BezierSurface::Intersect hides the virtual function 
  // G4Surface::Intersect(const G4Ray&), I changed the name of this
  // function
  // G4int Intersect(G4SurfaceList&); 
  G4int BIntersect(G4SurfaceList&);

  G4SurfaceList* bezier_list;
  int ClipBothDirs();
  void ClipSurface();

  virtual G4Vector3D SurfaceNormal(const G4Point3D& Pt)const
  {
    return G4Vector3D(0,0,0);
  }


private:

  int             order[2];
  G4double        smin;
  G4double        smax;
  G4Point3D       line;
  G4double        average_u;
  G4double        average_v;
  G4Point3D       average_pt;
  int             dir;    
  G4KnotVector    *u_knots;
  G4KnotVector    *v_knots;
  G4ControlPoints *ctl_points;

  void CalcAverage();
  void CalcDistance(const G4Point3D&);
  void SetValues();

  inline void LocalizeClipValues()
  {
    if ( dir == ROW)
    {
      smin = (1.0 - smin) * u_knots->GetKnot(0) +
	smin * u_knots->GetKnot(u_knots->GetSize() - 1);
      smax = (1.0 - smax) * u_knots->GetKnot(0) +
	smax * u_knots->GetKnot(u_knots->GetSize() - 1);
    }
    else
    {
      smin = (1.0 - smin) * v_knots->GetKnot(0) +
	smin * v_knots->GetKnot(v_knots->GetSize() - 1);
      smax = (1.0 - smax) * v_knots->GetKnot(0) +
	smax * v_knots->GetKnot(v_knots->GetSize() - 1);
    }  
  }
  
  G4KnotVector     *new_knots;
  int              ord;
  G4OsloMatrix     * oslo_m;
  int              lower,upper;
  G4double         u[2];
  G4double         v[2];
  G4double         u_min;
  G4double         u_max;
  G4double         v_min;
  G4double         v_max;
  G4ControlPoints* old_points;  

  void SplitNURBSurface();
  void GetClippedRegionFromSurface();
  void RefineSurface();
  void CalcOsloMatrix();
  void MapSurface(G4Surface*);
    
  // For ClipSurface...
  inline G4double Findzero(G4double x0,G4double x1,G4double y0,G4double y1) 
  {
    return(x0 - y0 * ( x1 - x0) / (y1-y0));
  };
  
  inline int Sign(G4double a)	
  {
    return((a < 0.0)? -1 : 1)	;
  };

  // For calc_G4OsloMatrix...
  inline int Amax(int i, int j)    {return( (i) > (j) ? (i) : (j) );};
  inline int Amin(int i, int j)    {return( (i) < (j) ? (i) : (j) );};
	
  inline int AhIndex(int j,int t, int iorder)   
  {
    return(( (j) * ((j)+1)/2) + (t) - ((iorder-1) - (j)));
  };

};

#endif
































