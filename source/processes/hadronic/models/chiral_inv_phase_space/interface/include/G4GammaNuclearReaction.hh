//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4GammaNuclearReaction.hh,v 1.7 2002-05-22 11:43:15 mkossov Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// GEANT4 physics class: G4GammaNuclearReaction -- header file
// Created: M.V. Kossov, CERN/ITEP(Moscow), 10-OCT-01
// The last update: M.V. Kossov, CERN/ITEP (Moscow) 21-May-02
//
#ifndef G4GammaNuclearReaction_h
#define G4GammaNuclearReaction_h

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4ChiralInvariantPhaseSpace.hh"
#include "G4Gamma.hh"

class G4GammaNuclearReaction : public G4HadronicInteraction
{
  public: 
    virtual ~G4GammaNuclearReaction()
    {
    }
    
    G4VParticleChange * ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus);

  private:
    G4ChiralInvariantPhaseSpace theModel;
};

inline
G4VParticleChange * G4GammaNuclearReaction::
ApplyYourself(const G4Track& aTrack, G4Nucleus& aTargetNucleus)
{
  if(aTrack.GetDynamicParticle()->GetDefinition() != G4Gamma::GammaDefinition())
  {
    G4Exception("Called G4GammaNuclearReaction for particle other than gamma");
  }
  return theModel.ApplyYourself(aTrack, aTargetNucleus);
}

#endif
