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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
#ifndef G4BinaryLightIonReaction_h
#define G4BinaryLightIonReaction_h

#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadFinalState.hh"
#include "G4ExcitationHandler.hh"

class G4BinaryLightIonReaction : public G4HadronicInteraction 
{
  public:
    G4BinaryLightIonReaction();
    virtual ~G4BinaryLightIonReaction(){}
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                              G4Nucleus& theNucleus);
  
  private:
    G4BinaryCascade theModel;
    G4ExcitationHandler theHandler;
    G4PreCompoundModel theProjectileFragmentation;
    G4HadFinalState theResult;
    G4bool EnergyAndMomentumCorrector(G4ReactionProductVector* products,
    				G4LorentzVector& TotalCollisionMom);
};

#endif
