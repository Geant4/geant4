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
// (Why only antiproton-proton, when the antiproton-nucleus is made? - M.K.)
// 17.02.2009 M.Kossov, now it is recommended to use the G4QCollision process
#ifndef G4ProtonAntiProtonReaction_h
#define G4ProtonAntiProtonReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4AntiProton.hh"

class G4ProtonAntiProtonReaction : public G4HadronicInteraction
{
  public: 
    virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
    G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
};

inline
G4VParticleChange * G4ProtonAntiProtonReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4AntiProton::AntiProton())
  {
    throw G4HadronicException(__FILE__, __LINE__, "Calling G4ProtonAntiProtonReaction with particle other than p-bar!!!");
  }
  if(aTargetNucleus.GetZ() != 1)
  {
    throw G4HadronicException(__FILE__, __LINE__, "Calling G4ProtonAntiProtonReaction for target other than Hydrogen!!!");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

#endif
