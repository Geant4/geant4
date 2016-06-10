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
// $Id: G4ErrorTarget.hh 66356 2012-12-18 09:02:32Z gcosmo $
//
//
// --------------------------------------------------------------------
//      GEANT 4 class header file 
// --------------------------------------------------------------------
//
// Class Description:
//
// Base class for all error propagation targets.

// History:
// - Created. P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorTarget_hh
#define G4ErrorTarget_hh

#include "globals.hh"
#include "G4ThreeVector.hh"

class G4Step;

enum G4ErrorTargetType{ G4ErrorTarget_PlaneSurface,
                        G4ErrorTarget_CylindricalSurface,
                        G4ErrorTarget_GeomVolume,
                        G4ErrorTarget_TrkL };
class G4ErrorTarget
{
  public:  // with description

    G4ErrorTarget();
    virtual ~G4ErrorTarget();

    virtual G4double GetDistanceFromPoint( const G4ThreeVector&,
                                           const G4ThreeVector& ) const;
    virtual G4double GetDistanceFromPoint( const G4ThreeVector& ) const;
      // For the target volume

    virtual G4bool TargetReached(const G4Step*);
      // For the target surface and target TrackLength

    virtual void Dump( const G4String& msg ) const = 0;

    // Access methods

    inline G4ErrorTargetType GetType() const;

  protected:

    G4ErrorTargetType theType;
};

// Inline methods

inline G4ErrorTargetType G4ErrorTarget::GetType() const
{
  return theType;
}

#endif
