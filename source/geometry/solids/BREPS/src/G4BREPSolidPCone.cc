//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: G4BREPSolidPCone.cc,v 1.24 2002-01-22 22:45:38 radoone Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// ----------------------------------------------------------------------
// GEANT 4 class source file
//
// G4BREPSolidPCone.cc
//
// ----------------------------------------------------------------------
// The polyconical solid G4BREPSolidPCone is a shape defined by a set of 
// inner and outer conical or cylindrical surface sections and two planes 
// perpendicular to the Z axis. Each conical surface is defined by its 
// radius at two different planes perpendicular to the Z-axis. Inner and 
// outer conical surfaces are defined using common Z planes. 
// ----------------------------------------------------------------------

#include "g4std/strstream"

#include "G4BREPSolidPCone.hh"
#include "G4FCylindricalSurface.hh"
#include "G4FConicalSurface.hh"
#include "G4CircularCurve.hh"
#include "G4FPlane.hh"

typedef enum {
  EInverse = 0,
  ENormal = 1
} ESurfaceSense;

G4BREPSolidPCone::G4BREPSolidPCone(const G4String& name,
				   G4double start_angle,
				   G4double opening_angle,
				   G4int    num_z_planes, // sections,
				   G4double z_start,		   
				   G4double z_values[],
				   G4double RMIN[],
				   G4double RMAX[] )
 : G4BREPSolid(name)

{
  G4int sections= num_z_planes-1;
  nb_of_surfaces = 2*sections+2;
  SurfaceVec = new G4Surface*[nb_of_surfaces];
  G4ThreeVector Axis(0,0,1);
  G4ThreeVector Origin(0,0,z_start);    
  G4double Length;
  G4ThreeVector LocalOrigin(0,0,z_start);    
  G4int a, b = 0;

  G4ThreeVector PlaneAxis(0, 0, 1);
  G4ThreeVector PlaneDir (0, 1, 0);   

  ///////////////////////////////////////////////////
  // Test the validity of the R values
  
  // RMIN[0] and RMIN[num_z_planes-1] cannot be = 0
  // when RMIN[0] or RMIN[num_z_planes-1] are = 0
  if(  ((RMIN[0] == 0)              && (RMAX[0] == 0))             || 
       ((RMIN[num_z_planes-1] == 0) && (RMAX[num_z_planes-1] == 0))   )
      G4Exception("G4BREPSolidPCone::G4BREPSolidPCone() - RMIN at the extremities cannot be 0 when RMAX = 0 !");  
  
  // only RMAX[0] and RMAX[num_z_planes-1] can be = 0
  for(a = 1; a < num_z_planes-1; a++)
    if (RMAX[a] == 0)
      G4Exception("G4BREPSolidPCone::G4BREPSolidPCone() - RMAX inside the solid cannot be 0 !");
  
  // RMAX[a] must be greater than RMIN[a] 
  for(a = 1; a < num_z_planes-1; a++)
    if (RMIN[a] >= RMAX[a])
      G4Exception("G4BREPSolidPCone::G4BREPSolidPCone() - RMAX must be greater than RMIN in the middle Z planes !");
  
  if( (RMIN[num_z_planes-1] > RMAX[num_z_planes-1] )
      || (RMIN[0] > RMAX[0]) ) 
      G4Exception("G4BREPSolidPCone::G4BREPSolidPCone() - RMAX must be greater or equal than RMIN at the ends !");

  // Create surfaces
  for( a = 0; a < sections; a++)
  {
    // Surface length
    Length = z_values[a+1] - z_values[a];

    if( Length == 0 )
    {
      // We need to create planar surface(s)
      if( RMAX[a] != RMAX[a+1] && RMIN[a] != RMIN[a+1] )
      {
        // We can have the 8 following configurations here:
        //
        // 1.     2.     3.     4.    
        // --+      +--  --+      +--   
        // xx|->  <-|xx  xx|      |xx   
        // xx+--  --+xx  --+      +--   
        // xxxxx  xxxxx    |      |     
        // xxxxx  xxxxx    +--  --+   
        // xx+--  --+xx    |xx  xx|   
        // xx|->  <-|xx    +--  --+   
        // --+      +--  
        // -------------------------- Z axis
        //
        //////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////
        //
        // 5.     6.     7.     8.
        // --+      +--  --+      +--
        // xx|->  <-|xx  xx|->  <-|xx
        // --+--  --+--  xx+--  --+xx
        // <-|xx  xx|->  xxxxx  xxxxx
        //   +--  --+    --+xx  xx+--
        //               <-|xx  xx|->
        //                 +--  --+  
        // -------------------------- Z axis
        //
        // NOTE: The pictures shows only one half of polycone!
        //       The arrows show the expected surface normal direction.
        //       The configuration No. 3 and 4 are not valid solids!
        
        // Eliminate the invalid cases 3 and 4.
        // At this point is guaranteed that each RMIN[i] < RMAX[i]
        // where i in in interval 0 < i < num_z_planes-1. So:
        if( RMIN[a] > RMAX[a+1] || RMAX[a] < RMIN[a+1] ) {
          char msgbuf[512];
          G4std::ostrstream os(msgbuf,512);
          os << G4endl << "G4BREPSolidPCone::G4BREPSolidPCone() - The values "
                       << "of RMIN[" << a << "] & RMAX[" << a+1 << "] or RMAX[" << a << "] & RMIN[" << a+1 << "] "
                       << "make an invalid configuration of G4BREPSolidPCone " << name.c_str() << "!" << G4endl;
          G4Exception( msgbuf );
        }
        
        ESurfaceSense UpSurfSense, LowSurfSense;
        
        // We need to clasify all the cases in order to figure out the planar surface sense
        if( RMAX[a] > RMAX[a+1] ) {
          // Cases 1, 5, 7
          if( RMIN[a] < RMIN[a+1] ) {
            // Case 1
            UpSurfSense  = ENormal;
            LowSurfSense = ENormal;
          } else if( RMAX[a+1] != RMIN[a]) {
            // Case 7
            UpSurfSense  = ENormal;
            LowSurfSense = EInverse;
          } else {
            // Case 5
            UpSurfSense  = ENormal;
            LowSurfSense = EInverse;
          }
        } else {
          // Cases 2, 6, 8
          if( RMIN[a] > RMIN[a+1] ) {
            // Case 2
            UpSurfSense  = EInverse;
            LowSurfSense = EInverse;
          } else if( RMIN[a+1] != RMAX[a] ) {
            // Case 8
            UpSurfSense  = EInverse;
            LowSurfSense = ENormal;
          } else {
            // Case 6
            UpSurfSense  = EInverse;
            LowSurfSense = ENormal;
          }
        }
            
        SurfaceVec[b]  = ComputePlanarSurface( RMAX[a], RMAX[a+1], LocalOrigin, PlaneAxis, PlaneDir, UpSurfSense );
        SurfaceVec[b++]->SetSameSense( UpSurfSense );
        
        SurfaceVec[b]  = ComputePlanarSurface( RMIN[a], RMIN[a+1], LocalOrigin, PlaneAxis, PlaneDir, LowSurfSense );
        SurfaceVec[b++]->SetSameSense( LowSurfSense );
      }
      else
      {
        // The original code creating single planar surface
        // in case where only either RMAX or RMIN have changed at the same Z value
        // e.g.:
        //
        //     RMAX          RMIN
        //    change        change
        //
        // 1      2      3      4
        // --+      +--  -----  -----
        // 00|->  <-|00  00000  00000
        // 00+--  --+00  --+00  00+--
        // 00000  00000  <-|00  00|->
        //                 +--  --+
        // --------------------------- Z axis
        //
        // NOTE: The picture shows only one half of polycone!
        
        G4double R1, R2;
        ESurfaceSense SurfSense;
        
        // test where is the plane surface
//         if( RMAX[a] != RMAX[a+1] ) {
// 	        R1 = RMAX[a];
//           R2 = RMAX[a+1];
//         } else if(RMIN[a] != RMIN[a+1]) {
// 	        R1 = RMIN[a];
//           R2 = RMIN[a+1];
//         } else   {
//           G4cerr << "Error in construction of G4BREPSolidPCone. \n"
//                  << "Exactly the same z, rmin and rmax given for \n"
//                  << "consecutive indices, " << a << " and " << a+1 << G4endl;
//           continue; 
//         }
        if( RMAX[a] != RMAX[a+1] ) {
	        // Cases 1, 2
          R1 = RMAX[a];
          R2 = RMAX[a+1];
          if( R1 > R2 ) {
            // Case 1
            SurfSense = ENormal;
          } else {
            // Case 2
            SurfSense = EInverse;
          }
        } else if(RMIN[a] != RMIN[a+1]) {
	        // Cases 1, 2
	        R1 = RMIN[a];
          R2 = RMIN[a+1];
          if( R1 > R2 ) {
            // Case 3
            SurfSense = EInverse;
          } else {
            // Case 4
            SurfSense = ENormal;
          }
        } else   {
          G4cerr << "Error in construction of G4BREPSolidPCone. \n"
                 << "Exactly the same z, rmin and rmax given for \n"
                 << "consecutive indices, " << a << " and " << a+1 << G4endl;
          continue; 
        }
        
        SurfaceVec[b] = ComputePlanarSurface( R1, R2, LocalOrigin, PlaneAxis, PlaneDir, SurfSense );
        SurfaceVec[b++]->SetSameSense( SurfSense );
        
        //      SurfaceVec[b]->SetSameSense(1);
        nb_of_surfaces--;
      }
    }
    else
    {
      // The surface to create is conical or cylindrical

      // Inner PCone
      if(RMIN[a] != RMIN[a+1])
      {
	      // Create cone
	      if(RMIN[a] > RMIN[a+1]) {
	        G4Vector3D ConeOrigin = G4Vector3D(LocalOrigin) ; 
	        SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis, Length, RMIN[a+1], RMIN[a]);
	        SurfaceVec[b]->SetSameSense(0);
	      } else {      
	        G4Vector3D Axis2 = G4Vector3D( -1*Axis );
	        G4Vector3D LocalOrigin2 = G4Vector3D( LocalOrigin + (Length*Axis) );
	        G4Vector3D ConeOrigin = LocalOrigin2; 
	        SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis2, Length, RMIN[a], RMIN[a+1]); 
	        SurfaceVec[b]->SetSameSense(0);
	      }
	      b++;  
      }
      else
      {
	      if( RMIN[a] == 0 ) {
	        // Do not create any surface and decrease nb_of_surfaces
	        nb_of_surfaces--;
	      } else {
	        // Create cylinder
	        G4Vector3D CylOrigin = G4Vector3D( LocalOrigin ); 
	        SurfaceVec[b] = new G4FCylindricalSurface(CylOrigin, Axis, RMIN[a], Length );
	        SurfaceVec[b]->SetSameSense(0);
	        b++;
	      }	    
      }

      // Outer PCone
      if(RMAX[a] != RMAX[a+1])
      {
	      // Create cone
	      if(RMAX[a] > RMAX[a+1]) {
	        G4Vector3D ConeOrigin = G4Vector3D( LocalOrigin );
	        SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis, Length, RMAX[a+1], RMAX[a]);
	        SurfaceVec[b]->SetSameSense(1);
	      } else {
	        G4Vector3D Axis2 = G4Vector3D( -1*Axis );
	        G4Vector3D LocalOrigin2 = G4Vector3D( LocalOrigin + (Length*Axis) );
	        G4Vector3D ConeOrigin = LocalOrigin2 ;

	        SurfaceVec[b] = new G4FConicalSurface(ConeOrigin, Axis2, Length, RMAX[a], RMAX[a+1]);
	        SurfaceVec[b]->SetSameSense(1);
	      }
	      b++;
      }
      else
      {
	      // Create cylinder
	      G4Vector3D CylOrigin = G4Vector3D( LocalOrigin );

	      if (RMAX[a] == 0) {
	        // Do not create any surface and decrease nb_of_surfaces
	        nb_of_surfaces--;
	      } else {
	        // Create cylinder
	        G4Vector3D CylOrigin = G4Vector3D( LocalOrigin ); 
	        SurfaceVec[b] = new G4FCylindricalSurface(CylOrigin, Axis, RMAX[a], Length );
	        SurfaceVec[b]->SetSameSense(1);
	        b++;
	      }	  	
      }
    }

    // Move surface origin to next section
    LocalOrigin = LocalOrigin + (Length*Axis);
  }
  
  ///////////////////////////////////////////////////
  // Create two end planes  

  G4CurveVector cv;
  G4CircularCurve* tmp;

  if(RMIN[0] < RMAX[0]) {

    // Create start G4Plane & boundaries    
    G4Point3D ArcStart1a = G4Point3D( Origin + (RMIN[0]*PlaneDir) );
    G4Point3D ArcStart1b = G4Point3D( Origin + (RMAX[0]*PlaneDir) );
 
    tmp = new G4CircularCurve;
    tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, Origin), RMIN[0]);
    tmp->SetBounds(ArcStart1a, ArcStart1a);
    tmp->SetSameSense(0);
    cv.push_back(tmp);
  
    tmp = new G4CircularCurve;
    tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, Origin), RMAX[0]);
    tmp->SetBounds(ArcStart1b, ArcStart1b);
    tmp->SetSameSense(1);
    cv.push_back(tmp);

    SurfaceVec[nb_of_surfaces-2] = new G4FPlane(PlaneDir, -PlaneAxis, Origin);
    SurfaceVec[nb_of_surfaces-2]->SetBoundaries(&cv);
    SurfaceVec[nb_of_surfaces-2]->SetSameSense(0);
  } else {
    // RMIN[0] == RMAX[0] so no surface is needed, it is a line!
    nb_of_surfaces--;    
  }

  if(RMIN[sections] < RMAX[sections]) {

    // Create end G4Plane & boundaries    
    G4Point3D ArcStart2a = G4Point3D( LocalOrigin+(RMIN[sections]*PlaneDir) );  
    G4Point3D ArcStart2b = G4Point3D( LocalOrigin+(RMAX[sections]*PlaneDir) );    
  
    cv.clear();

    tmp = new G4CircularCurve;
    tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), 
	      RMIN[sections]);
    tmp->SetBounds(ArcStart2a, ArcStart2a);
    tmp->SetSameSense(0);
    cv.push_back(tmp);
    
    tmp = new G4CircularCurve;
    tmp->Init(G4Axis2Placement3D(PlaneDir, PlaneAxis, LocalOrigin), 
	      RMAX[sections]);
    tmp->SetBounds(ArcStart2b, ArcStart2b);
    tmp->SetSameSense(1);
    cv.push_back(tmp);
  
    SurfaceVec[nb_of_surfaces-1]= new G4FPlane(PlaneDir, PlaneAxis, LocalOrigin);
    SurfaceVec[nb_of_surfaces-1]->SetBoundaries(&cv);
  
    // set sense of the surface
    SurfaceVec[nb_of_surfaces-1]->SetSameSense(0);
  } else {
    // RMIN[0] == RMAX[0] so no surface is needed, it is a line!
    nb_of_surfaces--;    
  }

  active=1;
  Initialize();

  // Store the original parameters, to be used in visualisation
  original_parameters.Start_angle= start_angle;
  original_parameters.Opening_angle= opening_angle;		   
  original_parameters.Num_z_planes= num_z_planes; 
  // original_parameters.z_start= z_start;		   
  original_parameters.Z_values= new G4double[num_z_planes];
  original_parameters.Rmin= new G4double[nb_of_surfaces];
  original_parameters.Rmax= new G4double[nb_of_surfaces];
  
  for(int is=0;is<num_z_planes;is++)
  {
    original_parameters.Z_values[is]= z_values[is]; 
    original_parameters.Rmin[is]= RMIN[is];
    original_parameters.Rmax[is]= RMAX[is];
  }

  // z_values[0]  should be equal to z_start, for consistency 
  //   with what the constructor does.
  // Otherwise the z_values that are given are used 
  //   shifted by   z_values[0] - z_start: 
  //  (because z_values are only used in 
  //    line 26:      Length = z_values[a+1] - z_values[a]; 
  //  )                                                      // JA Apr 2, 97
  
  /*
  if( z_values[0] != z_start )
  {
    G4cerr << "ERROR in creating G4BREPSolidPCone: "  
           << " z_values[0]= " << z_values[0] << " is not equal to " 
           << " z_start= " , z_start; 
    // G4Exception(" Error in creating G4BREPSolidPCone: z_values[0] must be equal to z_start" );
    original_parameters.Z_values[0]= z_start;
  }
  */ 
  
}

G4BREPSolidPCone::~G4BREPSolidPCone()
{
  delete [] original_parameters.Z_values;
  delete [] original_parameters.Rmin;
  delete [] original_parameters.Rmax;
}

void G4BREPSolidPCone::Initialize()
{ 
  // Calc bounding box for solids and surfaces
  // Convert concave planes to convex     
  ShortestDistance=1000000;
  CheckSurfaceNormals();
  
  if(!Box || !AxisBox)
    IsConvex();
  
  CalcBBoxes();
}

EInside G4BREPSolidPCone::Inside(register const G4ThreeVector& Pt) const
{
  // This function find if the point Pt is inside, 
  // outside or on the surface of the solid

  //  G4Vector3D v(1, 0, 0.01);
  //G4Vector3D v(1, 0, 0);
  G4Vector3D v(0, 0, 1);
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(v);
  G4Ray r(Pttmp, Vtmp);
  
  // Check if point is inside the PCone bounding box
  if( !GetBBox()->Inside(Pttmp) )
    return kOutside;

  // Set the surfaces to active again
  Reset();
  
  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  G4double dist = kInfinity;
  G4bool isIntersected = false;
  G4int WhichSurface = 0;

  // Chech if the point is on the surface, otherwise
  // find the nearest intersected suface. If there are not intersections the
  // point is outside

  for(G4int a=0; a < nb_of_surfaces; a++)
  {
    G4Surface* surface = SurfaceVec[a];
    
    if( surface->IsActive() ) {
      G4double hownear = surface->HowNear(Pt);
      
      if( fabs( hownear ) < kCarTolerance )
        return kSurface;

      if( surface->Intersect(r) ) {
        isIntersected = true;
        hownear = surface->GetDistance();
        
        if ( fabs( hownear ) < dist ) {
	        dist         = hownear;
	        WhichSurface = a;
        }
      }
    }
  }
  
  if ( !isIntersected )
    return kOutside; 

  // Find the point of intersection on the surface and the normal
  // !!!! be carefull the distance is sqrt(dist) !!!!

  dist = sqrt( dist );
  G4Vector3D IntersectionPoint = Pttmp + dist*Vtmp;
  G4Vector3D Normal = SurfaceVec[WhichSurface]->SurfaceNormal( IntersectionPoint );
  
  G4double dot = Normal*Vtmp;
  
  return( (dot > 0) ? kInside : kOutside );
}

G4ThreeVector G4BREPSolidPCone::SurfaceNormal(const G4ThreeVector& Pt) const
{  
  // This function calculates the normal of the surface
  // at a point on the surface
  // Note : the sense of the normal depends on the sense of the surface 

  G4ThreeVector   n(0,0,0);
  G4int        iplane;
  G4int        normflag = 0;
    
  G4Vector3D norm;
  G4Ray r( Pt, G4Vector3D(1, 0, 0) );

  // Find on which surface the point is
  for(iplane = 0; iplane < nb_of_surfaces; iplane++)
  {
    if( fabs(SurfaceVec[iplane]->HowNear(Pt)) < kCarTolerance) {
      // the point is on this surface
      normflag = 1;
      break;
    }
  }
  
  // calcul of the normal at this point
  if ( normflag ) {
    norm = SurfaceVec[iplane]->SurfaceNormal(Pt);

    n = G4ThreeVector ( norm.x(), norm.y(), norm.z() );
    n = n.unit();

    return n;
  } else {
    G4cout << "Warning ... PCone not able to return normal .. " << G4endl;
    return ( G4ThreeVector(1,0,0));
  }
}

G4double G4BREPSolidPCone::DistanceToIn(const G4ThreeVector& Pt) const
{
  // Calculates the shortest distance ("safety") from a point
  // outside the solid to any boundary of this solid.
  // Return 0 if the point is already inside.


  G4double *dists = new G4double[nb_of_surfaces];
  G4int a;

  // Set the surfaces to active again
  Reset();
  
  // calcul of the shortest distance of the point to each surfaces
  // Be carreful : it's a signed value
  for(a=0; a< nb_of_surfaces; a++)  
    dists[a] = SurfaceVec[a]->HowNear(Pt);
     
  G4double Dist = kInfinity;
  
  // if dists[] is positive, the point is outside
  // so take the shortest of the shortest positive distances
  // dists[] can be equal to 0 : point on a surface
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is inside to return 0, utilize the Inside
  //   function. But I don`t know if it is really needed because dToIn is 
  //   called only if the point is outside )
  for(a = 0; a < nb_of_surfaces; a++)
    if( fabs(Dist) > fabs(dists[a]) ) 
      //if( dists[a] >= 0)
	Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    // the point is inside the solid or on a surface
    return 0;
  else 
    //return Dist;
    return fabs(Dist);
}

G4double G4BREPSolidPCone::DistanceToIn(register const G4ThreeVector& Pt, 
					register const G4ThreeVector& V) const
{
  // Calculates the distance from a point outside the solid
  // to the solid`s boundary along a specified direction vector.
  // 	
  // Note : Intersections with boundaries less than the 
  //	    tolerance must be ignored if the direction 
  // 	    is away from the boundary
  
  G4int a;
  
  // Set the surfaces to active again
  Reset();
  
  G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
  
  for(a=0; a< nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->IsActive())
    {
      // test if the ray intersect the surface
      G4Vector3D Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);
      if( (Norm * Vtmp) < 0 && fabs(SurfaceVec[a]->HowNear(Pt)) < halfTolerance )
	return 0;
      if( (SurfaceVec[a]->Intersect(r)) ) {

	// if more than 1 surface is intersected,
	// take the nearest one
	if( SurfaceVec[a]->GetDistance() < ShortestDistance )
	  if( SurfaceVec[a]->GetDistance() > halfTolerance )
	    {
	      ShortestDistance = SurfaceVec[a]->GetDistance();
	    }
      }
    }
  }
  
  // Be careful !
  // SurfaceVec->Distance is in fact the squared distance
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  else
    // no intersection, return kInfinity
    return kInfinity; 
}

G4double G4BREPSolidPCone::DistanceToOut(register const G4ThreeVector& Pt, 
					 register const G4ThreeVector& V, 
					 const G4bool calcNorm, 
					 G4bool *validNorm, 
					 G4ThreeVector *n            ) const
{
  // Calculates the distance from a point inside the solid
  // to the solid`s boundary along a specified direction vector.
  // Return 0 if the point is already outside.
  //
  // Note : If the shortest distance to a boundary is less 
  // 	    than the tolerance, it is ignored. This allows
  // 	    for a point within a tolerant boundary to leave
  //	    immediately

  // Set the surfaces to active again
  Reset();

  const G4double halfTolerance = kCarTolerance*0.5;    
  G4Vector3D Ptv = G4Vector3D( Pt );
  G4int a;

  // I don`t understand this line
  if(validNorm)
    *validNorm=false;

  G4Vector3D Pttmp(Pt);
  G4Vector3D Vtmp(V);   
  
  G4Ray r(Pttmp, Vtmp);

  // Test if the bounding box of each surface is intersected
  // by the ray. If not, the surface become deactive.
  TestSurfaceBBoxes(r);
  
  ShortestDistance = kInfinity;
 
  for(a=0; a< nb_of_surfaces; a++)
  {
    if(SurfaceVec[a]->IsActive())
    {
      G4Vector3D Norm = SurfaceVec[a]->SurfaceNormal(Pttmp);
      if( (Norm * Vtmp) > 0 && fabs(SurfaceVec[a]->HowNear(Pt)) < halfTolerance )
        return 0;
      // test if the ray intersect the surface
      if( (SurfaceVec[a]->Intersect(r)) )
      {
	// if more than 1 surface is intersected,
	// take the nearest one
	if( SurfaceVec[a]->GetDistance() < ShortestDistance )
	  if( SurfaceVec[a]->GetDistance() > halfTolerance )
	  {
	    ShortestDistance = SurfaceVec[a]->GetDistance();
	  }
      }
    }
  }

  // Be carreful !
  // SurfaceVec->Distance is in fact the squared distance
  if(ShortestDistance != kInfinity)
    return sqrt(ShortestDistance);
  else
    // if no intersection is founded, the point is outside
    // so return 0
    return 0; 
}

G4double G4BREPSolidPCone::DistanceToOut(const G4ThreeVector& Pt) const
{
  // Calculates the shortest distance ("safety") from a point
  // inside the solid to any boundary of this solid.
  // Return 0 if the point is already outside.	

  G4double *dists = new G4double[nb_of_surfaces];
  G4int a;

  // Set the surfaces to active again
  Reset();
  
  // calcul of the shortest distance of the point to each surfaces
  // Be carreful : it's a signed value
  for(a=0; a< nb_of_surfaces; a++)
    dists[a] = SurfaceVec[a]->HowNear(Pt);  

  G4double Dist = kInfinity;
  
  // if dists[] is negative, the point is inside
  // so take the shortest of the shortest negative distances
  // dists[] can be equal to 0 : point on a surface
  // ( Problem with the G4FPlane : there is no inside and no outside...
  //   So, to test if the point is outside to return 0, utilize the Inside
  //   function. But I don`t know if it is really needed because dToOut is 
  //   called only if the point is inside )

  for(a = 0; a < nb_of_surfaces; a++)
    if( fabs(Dist) > fabs(dists[a]) ) 
      //if( dists[a] <= 0)
	Dist = dists[a];
  
  delete[] dists;

  if(Dist == kInfinity)
    // the point is ouside the solid or on a surface
    return 0;
  else
    // return Dist;
    return fabs(Dist);
}

void G4BREPSolidPCone::Reset() const
{
  Active(1);
  ((G4BREPSolidPCone*)this)->intersectionDistance=kInfinity;
  StartInside(0);
  for(register int a=0;a<nb_of_surfaces;a++)
    SurfaceVec[a]->Reset();
  ShortestDistance = kInfinity;
}

G4Surface* G4BREPSolidPCone::ComputePlanarSurface( G4double r1, G4double r2,
                                                         G4ThreeVector& origin,
                                                         G4ThreeVector& planeAxis,
                                                         G4ThreeVector& planeDirection,
                                                         G4int surfSense )
{
  // The planar surface to return
  G4Surface* planarFace = 0;
  
  G4CurveVector    cv1;
  G4CircularCurve  *tmp1, *tmp2;
  
  // Create plane surface
  G4Point3D ArcStart1 = G4Point3D( origin + (r1 * planeDirection) );
  G4Point3D ArcStart2 = G4Point3D( origin + (r2 * planeDirection) );

  if(r1 != 0)  
  {
	  tmp1 = new G4CircularCurve;
	  tmp1->Init(G4Axis2Placement3D( planeDirection, planeAxis, origin), r1);
	  tmp1->SetBounds(ArcStart1, ArcStart1);
    
	  if( surfSense )
	    tmp1->SetSameSense(1);
	  else
	    tmp1->SetSameSense(0);

	  cv1.push_back(tmp1);
  }

  if(r2 != 0)  
  {
	  tmp2 = new G4CircularCurve;
	  tmp2->Init(G4Axis2Placement3D( planeDirection, planeAxis, origin), r2);
	  tmp2->SetBounds(ArcStart2, ArcStart2);
    
	  if( surfSense )
	    tmp2->SetSameSense(0);
	  else
	    tmp2->SetSameSense(1);
    
	  cv1.push_back(tmp2);
  }

  planarFace   = new G4FPlane( planeDirection, planeAxis, origin);
  
  planarFace->SetBoundaries(&cv1);

  return planarFace;
}              

//  In graphics_reps:   

#include "G4Polyhedron.hh"   

G4Polyhedron* G4BREPSolidPCone::CreatePolyhedron() const
{
  return new G4PolyhedronPcon( original_parameters.Start_angle, 
			       original_parameters.Opening_angle, 
			       original_parameters.Num_z_planes, 
			       original_parameters.Z_values,
			       original_parameters.Rmin,
			       original_parameters.Rmax);
}
