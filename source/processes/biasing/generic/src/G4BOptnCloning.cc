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
#include "G4BOptnCloning.hh"


G4BOptnCloning::G4BOptnCloning(G4String name)
  : G4VBiasingOperation( name    ),
    fClone1W           ( -1.0    ),
    fClone2W           ( -1.0    ),
    fParticleChange(),
    fCloneTrack        ( nullptr )
{}

G4BOptnCloning::~G4BOptnCloning()
{}

G4VParticleChange*  G4BOptnCloning::GenerateBiasingFinalState( const G4Track* track,
                                                               const G4Step*       )
{
  fParticleChange.Initialize(*track);
  fParticleChange.ProposeParentWeight( fClone1W );
  fParticleChange.SetSecondaryWeightByProcess(true);
  fParticleChange.SetNumberOfSecondaries(1);
  fCloneTrack = new G4Track( *track );
  fCloneTrack->SetWeight( fClone2W );
  fParticleChange.AddSecondary( fCloneTrack );
  return &fParticleChange;
}
