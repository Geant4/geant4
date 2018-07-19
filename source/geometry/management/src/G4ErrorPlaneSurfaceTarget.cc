//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4ErrorPlaneSurfaceTarget.cc 107408 2017-11-10 13:35:09Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class implementation file 
// --------------------------------------------------------------------

#include "G4ErrorPlaneSurfaceTarget.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh" //for verbosity checking
#endif

#include "G4Point3D.hh"
#include "G4ThreeVector.hh"

//---------------------------------------------------------------------

G4ErrorPlaneSurfaceTarget::
G4ErrorPlaneSurfaceTarget(G4double aa, G4double ab, G4double ac, G4double ad)
  : G4Plane3D( aa, ab, ac, ad ) 
{
  theType = G4ErrorTarget_PlaneSurface;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 )
  { 
    Dump( " $$$ creating G4ErrorPlaneSurfaceTarget from parameters");
  }
#endif
}

//---------------------------------------------------------------------

G4ErrorPlaneSurfaceTarget::
G4ErrorPlaneSurfaceTarget(const G4Normal3D &norm, const G4Point3D &pt)
  : G4Plane3D( norm, pt ) 
{
  theType = G4ErrorTarget_PlaneSurface;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 )
  { 
    Dump( " $$$ creating G4ErrorPlaneSurfaceTarget from point and normal");
  }
#endif
}

//---------------------------------------------------------------------

G4ErrorPlaneSurfaceTarget::
G4ErrorPlaneSurfaceTarget(const G4Point3D &p1,
                          const G4Point3D &p2,
                          const G4Point3D &p3)
  : G4Plane3D( p1, p2, p3 )
{
  theType = G4ErrorTarget_PlaneSurface;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 )
  { 
    Dump( " $$$ creating G4ErrorPlaneSurfaceTarget from three points");
  }
#endif
}

//---------------------------------------------------------------------

G4ErrorPlaneSurfaceTarget::~G4ErrorPlaneSurfaceTarget()
{
}

//---------------------------------------------------------------------

G4ThreeVector G4ErrorPlaneSurfaceTarget::
Intersect( const G4ThreeVector& pt, const G4ThreeVector& dir ) const
{
  G4double lam = GetDistanceFromPoint( pt, dir );
  G4Point3D inters = pt + lam * dir;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 4 )
  { 
    G4cout << " $$$ creating G4ErrorPlaneSurfaceTarget::Intersect "
           << inters << G4endl;
  }
#endif

  return inters;
}

//---------------------------------------------------------------------

G4double G4ErrorPlaneSurfaceTarget::
GetDistanceFromPoint( const G4ThreeVector& pt, const G4ThreeVector& dir ) const
{
  if( std::fabs( dir.mag() -1. ) > 1.E-6 )
  {
    std::ostringstream message;
    message << "Direction is not a unit vector: " << dir << " !";
    G4Exception("G4ErrorPlaneSurfaceTarget::GetDistanceFromPoint()",
                "GeomMgt1002", JustWarning, message);
  }
  G4double dist = -(a_ * pt.x() + b_ * pt.y() + c_ * pt.z() + d_)
                 / (a_ * dir.x() + b_ * dir.y() + c_ * dir.z() );

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 )
  {
    G4cout << " G4ErrorPlaneSurfaceTarget::GetDistanceFromPoint()" << G4endl
           << "   Point: " << pt << ", Direction: " << dir << G4endl
           << "   Distance: " << dist << G4endl;
  }
#endif
  
  return dist;
}

//---------------------------------------------------------------------

G4double G4ErrorPlaneSurfaceTarget::
GetDistanceFromPoint( const G4ThreeVector& pt ) const
{
  G4ThreeVector vec = point() - pt;
  G4double dist = std::fabs(vec * normal() / normal().mag());
  
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 )
  {
    G4cout << " G4ErrorPlaneSurfaceTarget::GetDistanceFromPoint()" << G4endl
           << "   Point: " << pt << G4endl
           << "   Distance: " << dist << G4endl;
  }
#endif

  return dist;
}

//---------------------------------------------------------------------

G4Plane3D G4ErrorPlaneSurfaceTarget::
GetTangentPlane( const G4ThreeVector& ) const
{
  return *this;
}


void G4ErrorPlaneSurfaceTarget::Dump( const G4String& msg ) const
{
  G4cout << msg << " point = " << point()
                << " normal = " << normal() << G4endl;
}
