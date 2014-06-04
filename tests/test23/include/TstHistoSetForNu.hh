#ifndef TstHistoSetForNu_h
#define TstHistoSetForNu_h

#include "G4LorentzVector.hh"
// #include "TstReader.hh"

// fwd declaration
class G4Track;
class G4DecayProducts;

class TstHistoSetForNu
{

   public:

      enum NuERange        { kNone=-1, kR_0_2=0, kR_2_5=1, kR_5_10=2, kR_10_20=3, kR_20_50=4 };
      enum InteractionType { kOther=-1, kFromP=0, kFromPim=1, kFromPip=2 };
   
      TstHistoSetForNu() : fPionDecay(0) {}
      virtual ~TstHistoSetForNu() {}
            
   protected:
   
      void AccountForPionDecay( const G4Track* );

      NuERange         fNuERange;
      NuERange         fNubarERange;
      G4DecayProducts* fPionDecay;

};

#endif
