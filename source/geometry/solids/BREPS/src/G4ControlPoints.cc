// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ControlPoints.cc,v 1.5 2000-11-08 14:22:09 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4ControlPoints.cc
//
// ----------------------------------------------------------------------

#include "G4ControlPoints.hh"


G4ControlPoints::G4ControlPoints()
{
  nr=nc=0;
  data=(G4PointRat**)0;
}


G4ControlPoints::G4ControlPoints( G4int rows, G4int columns)
{
  nr=rows; 
  nc=columns; 
  data = (G4PointRat**) new G4PointRat *[nr*nc];
  
  for(G4int a =0; a<nr*nc;a++) 
    data[a]=new G4PointRat;
}


G4ControlPoints::G4ControlPoints( G4int point_type, G4int rows, G4int columns)
{

//     point_type is maintained only for compatibility 
//     G4ControlPoints is now a array of G4pointRat only

      nr=rows;
      nc=columns;
      data = (G4PointRat**)new G4PointRat *[nr*nc];

      for(G4int a = 0; a < nr*nc ; a++ )
	data[a]=new G4PointRat;
}


G4ControlPoints::G4ControlPoints(const G4ControlPoints& old_points)
{ 
   // copy constructor
  
  for( G4int i = 0; i < nr*nc; i++)
    delete data[i];
  delete[] data;

  nr   = old_points.nr;
  nc   = old_points.nc;
  data = (G4PointRat**)new G4PointRat *[nr*nc];
  
  G4int a, b;
  
  for (a = 0; a < nr*nc ; a++ )
    data[a] = new G4PointRat;
      
  for ( a = 0; a < nr ; a++ )
    for ( b = 0; b < nc ; b++ )
      put( a, b, old_points.GetRat(a,b));
}
    

G4ControlPoints::~G4ControlPoints()
{
  for( G4int a = 0; a < nr*nc; a++)
    delete data[a];
  
  delete[] data;
}


G4ControlPoints& G4ControlPoints::operator=(const G4ControlPoints& c)
{ 
   // assignment operator

  if (&c == this) return *this;

  for( G4int i = 0; i < nr*nc; i++)
    delete data[i];
  delete[] data;

  nr   = c.nr;
  nc   = c.nc;
  data = (G4PointRat**)new G4PointRat *[nr*nc];
  
  G4int a, b;
  
  for (a = 0; a < nr*nc ; a++ )
    data[a] = new G4PointRat;
      
  for ( a = 0; a < nr ; a++ )
    for ( b = 0; b < nc ; b++ )
      put( a, b, c.GetRat(a,b));

  return *this;
}


void G4ControlPoints::SetWeights(G4double* weights)
{
  for ( G4int a = 0; a < nr*nc; a++ )
    (data[a])->setW(weights[a]);
}


void G4ControlPoints::CalcValues ( G4double k1, G4double param, 
				   G4PointRat& pts1, G4double k2, 
				   G4PointRat& pts2               )
{
  pts2.setX(Calc(k1,param,pts1.x(),k2,pts2.x()));
  pts2.setY(Calc(k1,param,pts1.y(),k2,pts2.y()));
  pts2.setZ(Calc(k1,param,pts1.z(),k2,pts2.z()));   
  pts2.setW(Calc(k1,param,pts1.w(),k2,pts2.w()));
}
 
 
void G4ControlPoints::CalcValues(G4double k1, G4double param, G4Point3D& pts1,
				 G4double k2, G4Point3D& pts2)
{ 		
  pts2.setX(Calc(k1,param,pts1.x(),k2,pts2.x()));
  pts2.setY(Calc(k1,param,pts1.y(),k2,pts2.y()));
  pts2.setZ(Calc(k1,param,pts1.z(),k2,pts2.z()));
}


G4double G4ControlPoints::ClosestDistanceToPoint( const G4Point3D& Pt)
{
  // Square distance 
  
  G4double  PointDist=1.e20; 
  G4double  TmpDist;
  G4Point3D Pt2;
  
  for(G4int a=0;a<nr;a++)
    for(G4int b=0;b<nc;b++)
    {
      Pt2       = Get3D(a,b);
      TmpDist   = Pt.distance2(Pt2);
      PointDist = ( PointDist > TmpDist ) ? TmpDist : PointDist;
    }
  
  return PointDist;
}
