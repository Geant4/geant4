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
#ifndef G4BinaryLightIonReaction_h
#define G4BinaryLightIonReaction_h 1

#include "G4BinaryCascade.hh"
#include "G4PreCompoundModel.hh"
#include "G4HadFinalState.hh"
#include "G4ExcitationHandler.hh"

class G4BinaryLightIonReaction : public G4HadronicInteraction
{
  public:
    G4BinaryLightIonReaction(G4VPreCompoundModel* ptr = 0);
    virtual ~G4BinaryLightIonReaction();
    G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack,
                                              G4Nucleus& theNucleus);
    inline void SetPrecompound(G4VPreCompoundModel* ptr);
    inline void SetDeExcitation(G4ExcitationHandler* ptr);

    virtual void ModelDescription(std::ostream&) const ;

  private:
    G4bool EnergyAndMomentumCorrector(G4ReactionProductVector* products,
               G4LorentzVector& TotalCollisionMom);
    G4bool SetLighterAsProjectile(G4LorentzVector & mom,const G4LorentzRotation & toBreit);
    G4ReactionProductVector * FuseNucleiAndPrompound(const G4LorentzVector & mom);
    G4ReactionProductVector * Interact(G4LorentzVector & mom, const G4LorentzRotation & );
    G4double GetProjectileExcitation();
    void DeExciteSpectatorNucleus(G4ReactionProductVector * spectators, G4ReactionProductVector * cascaders,
                           G4double theStatisticalExEnergy, G4LorentzVector & momentum);
    G4LorentzVector SortResult(G4ReactionProductVector * result,G4ReactionProductVector * spectators,G4ReactionProductVector * cascaders);

    G4BinaryCascade* theModel;
    G4ExcitationHandler* theHandler;
    G4VPreCompoundModel* theProjectileFragmentation;
    G4HadFinalState theResult;
    G4int pA, pZ, tA, tZ,spectatorA,spectatorZ;
    G4Fancy3DNucleus * projectile3dNucleus, * target3dNucleus;
    G4FermiMomentum theFermi;
    G4LorentzVector pInitialState, pFinalState;

    G4bool debug_G4BinaryLightIonReactionResults;
    static G4int theBLIR_ID;
#ifdef G4MULTITHREADED
    static G4Mutex BLIRMutex;
#endif

};
inline void G4BinaryLightIonReaction::SetPrecompound(G4VPreCompoundModel* ptr)
{
  if(ptr) { theProjectileFragmentation = ptr; }
  theHandler = theProjectileFragmentation->GetExcitationHandler();
}
inline void G4BinaryLightIonReaction::SetDeExcitation(G4ExcitationHandler* ptr)
{
  theProjectileFragmentation->SetExcitationHandler(ptr);
  theHandler = ptr;
}

#endif
