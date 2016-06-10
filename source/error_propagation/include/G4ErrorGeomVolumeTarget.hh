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
// $Id: G4ErrorGeomVolumeTarget.hh 66892 2013-01-17 10:57:59Z gunter $
//
// Class Description:
//
// G4ErrorTarget class: limits step when volume is reached.

// History:
// - Created:  P. Arce, September 2004
// --------------------------------------------------------------------

#ifndef G4ErrorGeomVolumeTarget_HH
#define G4ErrorGeomVolumeTarget_HH

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4ErrorTarget.hh"
#include "G4Plane3D.hh"
  
class G4Step;
class G4String;

class G4ErrorGeomVolumeTarget : public G4ErrorTarget
{
 public:  // with description

  G4ErrorGeomVolumeTarget( const G4String& name );
  virtual ~G4ErrorGeomVolumeTarget(){}

  virtual G4bool TargetReached(const G4Step* aStep);
    // return true when particle is entering the volume

  virtual void Dump( const G4String& msg ) const;

 private:

  G4String  theName;
};

#endif
