#ifndef TstHistoSet_h
#define TstHistoSet_h

#include "G4LorentzVector.hh"
#include "TstReader.hh"

// fwd declaration
class G4VParticleChange;
class G4Track;
class G4DecayProducts;

class TstHistoSet
{

   public:

//      enum NuERange        { kNone=-1, kR_0_2=0, kR_2_5=1, kR_5_10=2, kR_10_20=3, kR_20_50=4 };
//      enum InteractionType { kOther=-1, kFromP=0, kFromPim=1, kFromPip=2 };
   
      TstHistoSet(G4String) : fDoResDecay(false) {}
      virtual ~TstHistoSet() {}
      
      virtual void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& ) = 0;
      virtual void Write( G4int, G4double ) = 0;
      
      void SetDoResDecay( G4bool rdec ) { fDoResDecay=rdec; return; }
      
   protected:
   
      void AccountForResDecay( G4VParticleChange* );
//      void AccountForPionDecay( const G4Track* );

      G4bool           fDoResDecay;
//      NuERange         fNuERange;
//      NuERange         fNubarERange;
//      G4DecayProducts* fPionDecay;

};

#endif
