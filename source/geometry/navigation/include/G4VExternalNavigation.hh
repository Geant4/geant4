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
// G4VExternalNavigation
//
// Class description:
//
// Pure virtual class to be specialised by the user for tracking with
// an external navigation

// Authors: V.Vlachoudis, G.Cosmo - CERN, 2019
// --------------------------------------------------------------------
#ifndef G4VEXTERNALNAVIGATION_HH
#define G4VEXTERNALNAVIGATION_HH

#include "G4LogicalVolume.hh"
#include "G4NavigationHistory.hh"
#include "G4ThreeVector.hh"
#include "G4VNavigation.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VSolid.hh"

class G4VExternalNavigation : public G4VNavigation
{
  public:  // with description

    G4VExternalNavigation();
      // Constructor
   
    virtual ~G4VExternalNavigation();
      // Destructor

    virtual G4VExternalNavigation* Clone() = 0;

    // Optional methods - may be necessary under particular circumstances

    virtual EInside Inside( const G4VSolid*      solid,
                            const G4ThreeVector& position,
                            const G4ThreeVector& direction );
     // Special 'Inside' call that includes direction of next motion
     //   provided for potential optimisations.

    virtual void RelocateWithinVolume( G4VPhysicalVolume*  motherPhysical,
                                       const G4ThreeVector& localPoint );
     //   Update any relevant internal state to take account that
     //      - the location has been moved to 'localPoint'
     //      - it remains in the current (mother) physical volume 'motherPhysical'
};

#endif
