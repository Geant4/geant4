#ifndef __MODEL_DEFINED__
#define __MODEL_DEFINED__

class G4Track;
class G4Nucleus;
class G4Fragment;
class G4VParticleChange;
class G4ExcitationHandler;

#include "G4HadronicInteraction.hh"
#include "globals.hh"
#include "INCModel.hh"
#include "G4Fragment.hh"

class Model : public G4HadronicInteraction
{
  G4ReactionProductVector* m_pResultVector;
  INCModel m_INCModel;
  G4ExcitationHandler* m_pExcitation;
public:
  Model(G4ExcitationHandler* pHandler): m_pResultVector(NULL),m_pExcitation(pHandler){};
  virtual ~Model();
  virtual G4VParticleChange* ApplyYourself(const G4Track& aTrack,G4Nucleus& aNucleus);
  G4ExcitationHandler* SetExcitationHandler(G4ExcitationHandler* phNew){
    G4ExcitationHandler* pRet = m_pExcitation;
    m_pExcitation = phNew;
    return pRet;
  }
  G4ExcitationHandler* GetExcitationHandler(){return m_pExcitation;};
private:
  double GetPotentialDiff(unsigned nPart);
  void ScaleRadial(unsigned nPart,double RefrCoef);
  double GetRefrCoef(unsigned nPart);
  double GetTau();
  void Emit(unsigned nPart);
  bool Propagate(unsigned nPart,double tau);
  void MainCycle();
  void PerformPreequilibrium();
  void DeExcite(G4Fragment& fragm);
  G4ReactionProductVector* DeExciteEq(G4Fragment& m_CurrFragm);
};

#endif
