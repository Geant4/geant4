#ifndef PCTBinaryCascade_hh
#define PCTBinaryCascade_hh

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
#include "G4HadFinalState.hh"

class G4CollisionManager;

class G4HadProjectile;
class G4KineticTrack;
class G43DNucleus;

class PCTProjectile;


class PCTBinaryCascade : public G4VIntraNuclearTransportModel
{
public:

  PCTBinaryCascade();
  PCTBinaryCascade(const PCTBinaryCascade & right);

  ~PCTBinaryCascade();

  const PCTBinaryCascade& operator=(PCTBinaryCascade & right);
  G4int operator==(PCTBinaryCascade& right) {return (this == &right);}
  G4int operator!=(PCTBinaryCascade& right) {return (this != &right);}

  G4HadFinalState* ApplyYourself(const G4HadProjectile & track,
				 G4Nucleus & nucleus);
  G4ReactionProductVector * DeExcite(const PCTProjectile * theProjectile,
				     const G4int A, const G4int Z);
  G4ReactionProductVector * Propagate(G4KineticTrackVector * secondaries,
				      G4V3DNucleus * nucleus);

private:

  G4int GetTotalCharge(std::vector<G4KineticTrack *>& aV)
  {
    G4int result = 0;
    std::vector<G4KineticTrack*>::iterator i;
    for (i = aV.begin(); i != aV.end(); ++i)
      {
        if  ((*i)->GetDefinition() == G4Proton::Proton())
          {
            ++result;
          }
      }
    return result;
  }

  void PrintWelcomeMessage();
  void BuildTargetList();
  void FindCollisions(G4KineticTrackVector * secondaries);
  G4bool ApplyCollision(G4CollisionInitialState * collision);
  G4bool Capture(G4bool verbose=false);
  G4bool Absorb();
  G4bool CheckPauliPrinciple(G4KineticTrackVector * products);
  G4double GetExcitationEnergy();
  G4bool CheckDecay(G4KineticTrackVector * products);
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
  G4ReactionProductVector * Propagate1H1(G4KineticTrackVector * secondaries,
                                         G4V3DNucleus * nucleus);

// utility methods
  G4ThreeVector GetSpherePoint(G4double r, const G4LorentzVector & momentumdirection);
  void ClearAndDestroy(G4KineticTrackVector * ktv);
  void ClearAndDestroy(G4ReactionProductVector * rpv);

// for debugging purpose
  void PrintKTVector(G4KineticTrackVector * ktv, std::string comment=std::string(""));

private:
  G4KineticTrackVector theProjectileList;
  G4KineticTrackVector theTargetList;
  G4KineticTrackVector theSecondaryList;
  G4KineticTrackVector theCapturedList;
  G4KineticTrackVector theFinalState;

  G4ExcitationHandler * theExcitationHandler;
  G4CollisionManager * theCollisionMgr;
  G4VScatterer * theScatterer;
  G4VFieldPropagation * thePropagator;
  G4double theCurrentTime;
  G4double theCutOnP;
  G4double theCutOnPAbsorb;
  G4LorentzVector theInitial4Mom;
  G4int currentA, currentZ;
  G4double massInNucleus;
  G4LorentzRotation precompoundLorentzboost;
  G4double theOuterRadius;
  G4bool thePrimaryEscape;
  const G4ParticleDefinition * thePrimaryType;
  G4ThreeVector theMomentumTransfer;


};

#endif




