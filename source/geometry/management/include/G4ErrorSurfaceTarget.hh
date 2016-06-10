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
// $Id: G4ErrorSurfaceTarget.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class header file 
// --------------------------------------------------------------------
//
// Class Description:
//
// Base class for G4ErrorTarget classes that are surfaces.

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorSurfaceTarget_hh
#define G4ErrorSurfaceTarget_hh

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorTanPlaneTarget.hh"
#include "G4Plane3D.hh"

class G4ErrorSurfaceTarget : public G4ErrorTanPlaneTarget
{
  public:  // with description

    G4ErrorSurfaceTarget();
    virtual ~G4ErrorSurfaceTarget();

    virtual double GetDistanceFromPoint( const G4ThreeVector& point,
                                         const G4ThreeVector& direc ) const = 0;
      // Distance from a point to the surface in a given direction

    virtual double GetDistanceFromPoint( const G4ThreeVector& point ) const = 0;
      // Minimal distance from a point to the surface in any direction

    virtual G4Plane3D GetTangentPlane( const G4ThreeVector& point ) const = 0;
      // Get tangent plane at point

    virtual void Dump( const G4String& msg ) const = 0;
      // Dump surface
};

#endif
