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
#include "G4ParticleDefinition.hh"
#include "g4rw/tpordvec.h"
#include "globals.hh"

class G4SPBaryon
{
  public:
    G4SPBaryon(G4Proton * aProton);
    G4SPBaryon(G4Neutron * aNeutron);
    G4SPBaryon(G4SigmaPlus * aSigmaPlus);
    G4SPBaryon(G4SigmaZero * aSigmaZero);
    G4SPBaryon(G4SigmaMinus * aSigmaMinus);
    G4SPBaryon(G4XiMinus * aXiMinus);
    G4SPBaryon(G4XiZero * aXiZero);
    G4SPBaryon(G4Lambda * aLambda);
    G4SPBaryon(G4OmegaMinus * anOmegaMinus);
    G4SPBaryon(G4ParticleDefinition * aDefinition);
    
    G4bool operator == ( const G4SPBaryon & aBaryon)
    {return this == &aBaryon; }
    
    G4ParticleDefinition * GetDefinition() {return theDefinition;}
    void SampleQuarkAndDiquark(G4int & quark, G4int & diQuark) const;
    void FindDiquark(G4int quark, G4int & diQuark) const;
    
  private:

    class G4SPPartonInfo
    {
      public:
        G4SPPartonInfo(G4int diq, G4int q, G4double prob) 
        { diQuarkPDGCode = diq; quarkPDGCode = q; probability = prob; }
        G4int GetQuark() const {return quarkPDGCode;}
        G4int GetDiQuark() const {return diQuarkPDGCode;}
        G4double GetProbability() const {return probability;}      
        G4bool operator == (const G4SPPartonInfo & aInfo) const
        {return this == &aInfo;}
      private:      
        G4int quarkPDGCode;
        G4int diQuarkPDGCode;
        G4double probability;
    };
    
  private:
  
  G4ParticleDefinition * theDefinition;
  G4RWTPtrOrderedVector<G4SPPartonInfo> thePartonInfo;
};

#endif
