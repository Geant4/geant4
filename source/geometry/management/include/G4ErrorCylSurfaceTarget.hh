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
// $Id: G4ErrorCylSurfaceTarget.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class header file 
// --------------------------------------------------------------------
//
// Class Description:
//
// G4ErrorTarget class: limits step when track reaches a cylindrical surface.

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorCylSurfaceTarget_hh
#define G4ErrorCylSurfaceTarget_hh

#include "globals.hh"
#include "G4ErrorSurfaceTarget.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4AffineTransform.hh"
#include "G4Plane3D.hh"

class G4ErrorCylSurfaceTarget : public G4ErrorSurfaceTarget
{
  public:  // with description

    G4ErrorCylSurfaceTarget( const G4double& radius,
                             const G4ThreeVector& trans=G4ThreeVector(),
                             const G4RotationMatrix& rotm=G4RotationMatrix() );
      // Constructs cylindrical surface by radius, translation and rotation

    G4ErrorCylSurfaceTarget( const G4double& radius,
                             const G4AffineTransform& trans );
      // Constructs cylindrical surface by radius and affine transformation

    ~G4ErrorCylSurfaceTarget();

    virtual G4ThreeVector IntersectLocal( const G4ThreeVector& point,
                                          const G4ThreeVector& direc ) const;
      // Intersects the cylindrical surface with the line in local (cylinder)
      // coordinates, given by point and direction

    virtual G4double GetDistanceFromPoint( const G4ThreeVector& point,
                                           const G4ThreeVector& direc ) const;
      // Distance from a point to the cylindrical surface in a given direction

    virtual G4double GetDistanceFromPoint( const G4ThreeVector& point ) const;
      // Minimal distance from a point to the cylindrical surface in any
      // direction

    virtual G4Plane3D GetTangentPlane( const G4ThreeVector& point ) const;
      // Get plane tangent to cylindrical surface at a given point

    virtual void Dump( const G4String& msg ) const;
      // Dump cylindrical surface parameter

  private:

     G4double fradius;
     G4AffineTransform ftransform;
};

#endif
