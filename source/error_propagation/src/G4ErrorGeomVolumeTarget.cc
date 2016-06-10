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
// $Id: G4ErrorGeomVolumeTarget.cc 66892 2013-01-17 10:57:59Z gunter $
//
// ------------------------------------------------------------
//      GEANT 4 class implementation file 
// ------------------------------------------------------------
//

#include "G4ErrorGeomVolumeTarget.hh"
#include "G4Point3D.hh"
#include "G4ThreeVector.hh"
#include "G4Step.hh"

#ifdef G4VERBOSE
#include "G4ErrorPropagatorData.hh" //for verbosity checking
#endif

//------------------------------------------------------------------------
G4ErrorGeomVolumeTarget::G4ErrorGeomVolumeTarget( const G4String& name )
{
  theType = G4ErrorTarget_GeomVolume;
  theName = name; 
}


//------------------------------------------------------------------------
G4bool G4ErrorGeomVolumeTarget::TargetReached( const G4Step* aStep )
{
  if( aStep->GetTrack()->GetNextVolume() != 0 ){
#ifdef G4VERBOSE
    if(G4ErrorPropagatorData::verbose() >= 3 ) { 
      G4cout << " G4ErrorGeomVolumeTarget::TargetReached( "
             << aStep->GetTrack()->GetNextVolume()->GetName()
             << " =? " <<  theName  << G4endl;
    }
#endif
    if( aStep->GetTrack()->GetNextVolume()->GetName() == theName ){
      return 1;
    } else {
      return 0;
    }
  } else {
    return 0;
  }
}


//------------------------------------------------------------------------
void G4ErrorGeomVolumeTarget::Dump( const G4String& msg ) const
{
  G4cout << msg << " G4ErrorGeomVolumeTarget:  Volume " << theName << G4endl;
}
