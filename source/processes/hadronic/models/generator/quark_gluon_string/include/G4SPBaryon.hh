#ifndef G4SPBaryon_h
#define G4SPBaryon_h

#include "globals.hh"

#include "G4Proton.hh"
#include "G4Neutron.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaZero.hh"
#include "G4SigmaMinus.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4Lambda.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiProton.hh"
#include "G4AntiNeutron.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiLambda.hh"
#include "G4AntiOmegaMinus.hh"

#include "G4ParticleDefinition.hh"
#include "G4SPPartonInfo.hh"
#include "g4std/vector"
#include "globals.hh"

class G4SPBaryon
{
  public:
    G4SPBaryon(G4Proton * aProton);
    G4SPBaryon(G4Neutron * aNeutron);
    G4SPBaryon(G4Lambda * aLambda);
    G4SPBaryon(G4SigmaPlus * aSigmaPlus);
    G4SPBaryon(G4SigmaZero * aSigmaZero);
    G4SPBaryon(G4SigmaMinus * aSigmaMinus);
    G4SPBaryon(G4XiMinus * aXiMinus);
    G4SPBaryon(G4XiZero * aXiZero);
    G4SPBaryon(G4OmegaMinus * anOmegaMinus);
    
    G4SPBaryon(G4AntiProton * aAntiProton);
    G4SPBaryon(G4AntiNeutron * aAntiNeutron);
    G4SPBaryon(G4AntiLambda * aAntiLambda);
    G4SPBaryon(G4AntiSigmaPlus * aAntiSigmaPlus);
    G4SPBaryon(G4AntiSigmaZero * aAntiSigmaZero);
    G4SPBaryon(G4AntiSigmaMinus * aAntiSigmaMinus);
    G4SPBaryon(G4AntiXiMinus * aAntiXiMinus);
    G4SPBaryon(G4AntiXiZero * aAntiXiZero);
    G4SPBaryon(G4AntiOmegaMinus * anAntiOmegaMinus);
    
    G4SPBaryon(G4ParticleDefinition * aDefinition);
    
    G4bool operator == ( const G4SPBaryon & aBaryon) const
    {return this == &aBaryon; }
    
    G4ParticleDefinition * GetDefinition() {return theDefinition;}
    void SampleQuarkAndDiquark(G4int & quark, G4int & diQuark) const;
    void FindDiquark(G4int quark, G4int & diQuark) const;
    G4int FindQuark(G4int diQuark) const;
    G4double GetProbability(G4int diQuark) const;
    G4int MatchDiQuarkAndGetQuark(const G4SPBaryon & aBaryon, G4int & aDiQuark) const;
        
  private:
  
  G4ParticleDefinition * theDefinition;
  G4std::vector<G4SPPartonInfo *> thePartonInfo;
};

#endif
