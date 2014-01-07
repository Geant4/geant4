#ifndef TstHistoSet_h
#define TstHistoSet_h

#include "G4LorentzVector.hh"
#include "TstReader.hh"

// fwd declaration
class G4VParticleChange;

class TstHistoSet
{

   public:
   
      TstHistoSet(G4String) : fDoResDecay(false) {}
      virtual ~TstHistoSet() {}
      
      virtual void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& ) = 0;
      virtual void Write( G4int, G4double ) = 0;
      
      void SetDoResDecay( G4bool rdec ) { fDoResDecay=rdec; return; }
      
   protected:
   
      void AccountForResDecay( G4VParticleChange* );

      G4bool fDoResDecay;

};

#endif
