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
//
// -------------------------------------------------------------------
//      GEANT4 Class file
//
//
//      File name:    G4BinaryCascade.hh 
//
//      Author: Alessandro Brunengo (Alessandro.Brunengo@ge.infn.it)       
// 
//      Creation date: 8 June 2000
// -----------------------------------------------------------------------------

#ifndef G4BinaryCascade_hh
#define G4BinaryCascade_hh

#include "G4VIntraNuclearTransportModel.hh"
#include "G4ReactionProductVector.hh"
#include "G4KineticTrackVector.hh"
#include "G4ListOfCollisions.hh"
#include "G4V3DNucleus.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4Fragment.hh"
#include "G4VFieldPropagation.hh"
#include "G4VScatterer.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"

#include "G4BCDecay.hh"
#include "G4BCLateParticle.hh"
#include "G4BCAction.hh"

class G4CollisionManager;

class G4Track;
class G4KineticTrack;
class G43DNucleus;
class G4Scatterer;

class G4BinaryCascade : public G4VIntraNuclearTransportModel
{
public:

  G4BinaryCascade();
  G4BinaryCascade(const G4BinaryCascade & right);

  virtual ~G4BinaryCascade();

  const G4BinaryCascade& operator=(G4BinaryCascade & right);
  G4int operator==(G4BinaryCascade& right) {return (this == &right);}
  G4int operator!=(G4BinaryCascade& right) {return (this != &right);}

  G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, 
                                              G4Nucleus& theNucleus);
  virtual G4ReactionProductVector * Propagate(G4KineticTrackVector *,
					      G4V3DNucleus *);

private:

  G4int GetTotalCharge(std::vector<G4KineticTrack *> & aV)
  {
    G4int result = 0;
    std::vector<G4KineticTrack *>::iterator i;
    for(i = aV.begin(); i != aV.end(); ++i)
    {
       result += G4lrint((*i)->GetDefinition()->GetPDGCharge());
    }
    return result;
  }
  G4int GetTotalBaryonCharge(std::vector<G4KineticTrack *> & aV)
  {
    G4int result = 0;
    std::vector<G4KineticTrack *>::iterator i;
    for(i = aV.begin(); i != aV.end(); ++i)
    {
       if ( (*i)->GetDefinition()->GetBaryonNumber() != 0 ){ 
             result += G4lrint((*i)->GetDefinition()->GetPDGCharge());
       }
    }
    return result;
  }
  void PrintWelcomeMessage();
  void BuildTargetList();
  void FindCollisions(G4KineticTrackVector *);
  void FindDecayCollision(G4KineticTrack *);
  void FindLateParticleCollision(G4KineticTrack *);
  G4bool ApplyCollision(G4CollisionInitialState *);
  G4bool Capture(G4bool verbose=false);
  G4bool Absorb();
  G4bool CheckPauliPrinciple(G4KineticTrackVector *);
  G4double GetExcitationEnergy();
  G4bool CheckDecay(G4KineticTrackVector *);
  void CorrectFinalPandE();
  void UpdateTracksAndCollisions(G4KineticTrackVector * oldSecondaries,
				 G4KineticTrackVector * oldTarget,
				 G4KineticTrackVector * newSecondaries);
  G4bool DoTimeStep(G4double timeStep);
  G4KineticTrackVector* CorrectBarionsOnBoundary(G4KineticTrackVector *in, 
                                               G4KineticTrackVector *out);
  G4Fragment * FindFragments();
  void StepParticlesOut();
  G4LorentzVector GetFinal4Momentum();
  G4LorentzVector GetFinalNucleusMomentum();
  G4ReactionProductVector * Propagate1H1(G4KineticTrackVector *,
					 G4V3DNucleus *);
  G4double GetIonMass(G4int Z, G4int A);
  
// utility methods
  G4ThreeVector GetSpherePoint(G4double r, const G4LorentzVector & momentumdirection);
  void ClearAndDestroy(G4KineticTrackVector * ktv);
  void ClearAndDestroy(G4ReactionProductVector * rpv);

// for debugging purpose
  void PrintKTVector(G4KineticTrackVector * ktv, std::string comment=std::string(""));
  void PrintKTVector(G4KineticTrack* kt, std::string comment=std::string(""));
  void DebugApplyCollision(G4CollisionInitialState * collision, 
                           G4KineticTrackVector *products);
  void DebugEpConservation(const G4HadProjectile & aTrack, G4ReactionProductVector* products);			   
			   
private:
  G4KineticTrackVector theProjectileList;
  G4KineticTrackVector theTargetList;
  G4KineticTrackVector theSecondaryList;
  G4KineticTrackVector theCapturedList;
  G4KineticTrackVector theFinalState;

  G4ExcitationHandler * theExcitationHandler;
  G4CollisionManager * theCollisionMgr;
  
  G4Scatterer * theH1Scatterer;

  std::vector<G4BCAction *> theImR;
  G4BCDecay * theDecay;
  G4BCLateParticle * theLateParticle;
  G4VFieldPropagation * thePropagator;
  G4double theCurrentTime;
  G4double theBCminP;
  G4double theCutOnP;
  G4double theCutOnPAbsorb;
  G4LorentzVector theInitial4Mom;
  G4int currentA, currentZ;
  G4double massInNucleus;
  G4double currentInitialEnergy;   // for debugging
  G4LorentzRotation precompoundLorentzboost;
  G4double theOuterRadius;
  G4bool thePrimaryEscape;
  G4ParticleDefinition * thePrimaryType;
  G4ThreeVector theMomentumTransfer;
  
  

};

#endif




