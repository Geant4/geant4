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
// G4ErrorTarget
//
// Class Description:
//
// Base class for all error propagation targets.

// Author: Pedro Arce (CIEMAT), September 2004
// --------------------------------------------------------------------
#ifndef G4ERRORTARGET_HH
#define G4ERRORTARGET_HH

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Step;

enum G4ErrorTargetType{ G4ErrorTarget_PlaneSurface,
                        G4ErrorTarget_CylindricalSurface,
                        G4ErrorTarget_GeomVolume,
                        G4ErrorTarget_TrkL };
/**
 * @brief G4ErrorTarget is base class for all error propagation targets.
 */

class G4ErrorTarget
{
  public:

    /**
     * Default Constructor and Destructor.
     */
    G4ErrorTarget() = default;
    virtual ~G4ErrorTarget() = default;

    /**
     * Methods to compute the distance from the target volume.
     */
    virtual G4double GetDistanceFromPoint( const G4ThreeVector&,
                                           const G4ThreeVector& ) const;
    virtual G4double GetDistanceFromPoint( const G4ThreeVector& ) const;

    /**
     * Returns true if the target surface / target track length is reached.
     */
    virtual G4bool TargetReached(const G4Step*);
      // 

    /**
     * Dumps parameters to standard output.
     */
    virtual void Dump( const G4String& msg ) const = 0;

    /**
     * Returns the type ID of the target.
     */
    inline G4ErrorTargetType GetType() const { return theType; }

  protected:

    G4ErrorTargetType theType{G4ErrorTarget_GeomVolume};
};

#endif
