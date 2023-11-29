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
// G4ErrorCylSurfaceTarget class implementation
// 
// Created: P.Arce, September 2004
// --------------------------------------------------------------------

#include "G4ErrorCylSurfaceTarget.hh"
#include "G4GeometryTolerance.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh" // for verbosity checking
#endif

#include "geomdefs.hh"
#include "G4Normal3D.hh"
#include "G4Plane3D.hh"

//---------------------------------------------------------------------

G4ErrorCylSurfaceTarget::
G4ErrorCylSurfaceTarget( const G4double& radius,
                         const G4ThreeVector& trans,
                         const G4RotationMatrix& rotm )
  : fradius(radius)
{
  theType = G4ErrorTarget_CylindricalSurface;

  ftransform = G4AffineTransform( rotm.inverse(), -trans );
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 )
  { 
    Dump( " $$$ creating G4ErrorCylSurfaceTarget ");
  }
#endif
}

//---------------------------------------------------------------------

G4ErrorCylSurfaceTarget::
G4ErrorCylSurfaceTarget( const G4double& radius,
                         const G4AffineTransform& trans )
  : fradius(radius), ftransform(trans.Inverse())
{
  theType = G4ErrorTarget_CylindricalSurface;

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 2 )
  { 
    Dump( " $$$ creating G4ErrorCylSurfaceTarget ");
  }
#endif
}

//---------------------------------------------------------------------

G4double G4ErrorCylSurfaceTarget::
GetDistanceFromPoint( const G4ThreeVector& point,
                      const G4ThreeVector& dir ) const
{
  if( dir.mag() == 0. )
  {
    G4Exception("G4ErrorCylSurfaceTarget::GetDistanceFromPoint()",
                "GeomMgt0003", FatalException, "Direction is zero !");
  }

  //----- Get intersection point
  G4ThreeVector localPoint = ftransform.TransformPoint( point );
  G4ThreeVector localDir = ftransform.TransformAxis( dir ); 
  G4ThreeVector inters = IntersectLocal(localPoint, localDir);

  G4double dist = (localPoint-inters).mag();
  
#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 )
  { 
    G4cout << " G4ErrorCylSurfaceTarget::GetDistanceFromPoint():" << G4endl
           << " Global point " << point << " dir " << dir << G4endl
           << " Intersection " << inters << G4endl
           << " Distance " << dist << G4endl;
    Dump( " CylSurface: " );
  }
#endif

  return dist;
}


//---------------------------------------------------------------------

G4double G4ErrorCylSurfaceTarget::
GetDistanceFromPoint( const G4ThreeVector& point ) const
{
  G4ThreeVector localPoint = ftransform.TransformPoint( point );

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 3 )
  { 
    G4cout << " G4ErrorCylSurfaceTarget::GetDistanceFromPoint:" << G4endl
           << " Global point " << point << G4endl
           << " Distance " << fradius - localPoint.perp() << G4endl;
    Dump( " CylSurface: " );
  }
#endif

  return fradius - localPoint.perp();
}

//---------------------------------------------------------------------

G4ThreeVector G4ErrorCylSurfaceTarget::
IntersectLocal( const G4ThreeVector& localPoint,
                const G4ThreeVector& localDir ) const
{
  G4double eqa = localDir.x()*localDir.x()+localDir.y()*localDir.y();
  G4double eqb = 2*(localPoint.x()*localDir.x()+localPoint.y()*localDir.y());
  G4double eqc = -fradius*fradius+localPoint.x()*localPoint.x()
                 +localPoint.y()*localPoint.y();
  G4int inside = (localPoint.perp() > fradius) ? -1 : 1;
  G4double lambda;

  if( eqa*inside > 0. )
  {
    lambda = (-eqb + std::sqrt(eqb*eqb-4*eqa*eqc) ) / (2.*eqa);
  }
  else if( eqa*inside < 0. )
  {
    lambda = (-eqb - std::sqrt(eqb*eqb-4*eqa*eqc) ) / (2.*eqa);
  }
  else
  {
    if( eqb != 0. )
    {
      lambda = -eqc/eqb;
    }
    else
    {
      std::ostringstream message;
      message << "Intersection not possible !" << G4endl
              << "          Point: " << localPoint << ", direction: "
              << localDir;
      Dump( " CylSurface: " ); 
      G4Exception("G4ErrorCylSurfaceTarget::IntersectLocal()",
                  "GeomMgt1002", JustWarning, message);
      lambda = kInfinity;
    }
  }

  G4ThreeVector inters = localPoint + lambda*localDir/localDir.mag();

#ifdef G4VERBOSE
  if(G4ErrorPropagatorData::verbose() >= 4 ) { 
    G4cout << " G4ErrorCylSurfaceTarget::IntersectLocal " << inters << " "
           << inters.perp() << " localPoint " << localPoint << " localDir "
           << localDir << G4endl;
  } 
#endif

  return inters;
}

//---------------------------------------------------------------------

G4Plane3D G4ErrorCylSurfaceTarget::
GetTangentPlane( const G4ThreeVector& point ) const
{
  G4ThreeVector localPoint = ftransform.TransformPoint( point );

  // check that point is at cylinder surface
  //
  if( std::fabs( localPoint.perp() - fradius )
      > 1000.*G4GeometryTolerance::GetInstance()->GetSurfaceTolerance() )
  {
    std::ostringstream message;
    message << "Local point not at surface !" << G4endl
            << "          Point: " << point << ", local: " << localPoint
            << G4endl
            << "          is not at surface, but far away by: "
            << localPoint.perp() - fradius << " !";
    G4Exception("G4ErrorCylSurfaceTarget::GetTangentPlane()",
                "GeomMgt1002", JustWarning, message);
  }

  G4Normal3D normal = localPoint - ftransform.NetTranslation();

  return G4Plane3D( normal, point );
}


//---------------------------------------------------------------------

void G4ErrorCylSurfaceTarget::Dump( const G4String& msg ) const
{
  G4cout << msg << " radius " << fradius
                << " centre " << ftransform.NetTranslation()
                << " rotation " << ftransform.NetRotation() << G4endl;
}
