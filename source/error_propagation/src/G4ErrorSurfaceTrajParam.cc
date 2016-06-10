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
// $Id: G4ErrorSurfaceTrajParam.cc 69014 2013-04-15 09:42:51Z gcosmo $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorSurfaceTrajParam.hh"
#include <iomanip>

#include "G4ThreeVector.hh"
#include "G4GeometryTolerance.hh"

//------------------------------------------------------------------------
G4ErrorSurfaceTrajParam::
G4ErrorSurfaceTrajParam( const G4Point3D& pos, const G4Vector3D& mom,
                         const G4Vector3D& vecV, const G4Vector3D& vecW )
{
  SetParameters( pos, mom, vecV, vecW );
}


//------------------------------------------------------------------------
G4ErrorSurfaceTrajParam::
G4ErrorSurfaceTrajParam( const G4Point3D& pos, const G4Vector3D& mom,
                         const G4Plane3D& plane )
{
  SetParameters( pos, mom, plane );
}

//------------------------------------------------------------------------
void G4ErrorSurfaceTrajParam::
SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
               const G4Plane3D& plane )
{
  //--- Get two perpendicular vectors: first parallel X
  //    (unless normal is parallel to X, then take Y)

  G4double kCarTolerance =
    G4GeometryTolerance::GetInstance()->GetSurfaceTolerance();

  G4ThreeVector Xvec(1.,0.,0.);
  G4Vector3D vecV = -Xvec.cross(plane.normal());
  if( vecV.mag() < kCarTolerance )
  {
    G4ThreeVector Zvec(0.,0.,1.);
    vecV = Zvec.cross(plane.normal());
  }

  G4Vector3D vecW = plane.normal().cross( vecV );

  SetParameters( pos, mom, vecV, vecW );
}


//------------------------------------------------------------------------
void G4ErrorSurfaceTrajParam::
SetParameters( const G4Point3D& pos, const G4Vector3D& mom,
               const G4Vector3D& vecV, const G4Vector3D& vecW )
{
  if( mom.mag() > 0. ) {
    fDir = mom;
    fDir /= mom.mag();
  } else {
    fDir = G4Vector3D(0.,0.,0.);
  }
  fVectorV = vecV / vecV.mag();
  fVectorW = vecW / vecW.mag();
  fInvP = 1./mom.mag();
  G4ThreeVector momv(mom);
  //check 3 vectors are ortogonal and right handed

  // now all 4 scalar memeber variables retain the signs
  //  fPV = momv.project( vecV ).mag();
  //  fPW = momv.project( vecW ).mag();
  fPV = momv.dot( vecV );
  fPW = momv.dot( vecW );


  G4ThreeVector posv(pos);
  // fV = posv.project( vecV ).mag();
  // fW = posv.project( vecW ).mag();
  fV = posv.dot( vecV );
  fW = posv.dot( vecW );
}


//------------------------------------------------------------------------
std::ostream& operator<<(std::ostream& out, const G4ErrorSurfaceTrajParam& tp)
{
  //  long mode = out.setf(std::ios::fixed,std::ios::floatfield);
  
  //  out << tp.theType;
  //  out << std::setprecision(5) << std::setw(10);
  out << " InvP= " << tp.fInvP << " PV= " << tp.fPV
      << " PW= " << tp.fPW << " V= " << tp.fV << " W= " << tp.fW << G4endl;
  out << " vectorV direction= " << tp.fVectorV
      << " vectorW direction= " << tp.fVectorW << G4endl;
    
  return out;
}
