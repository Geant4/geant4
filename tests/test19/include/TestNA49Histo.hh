#ifndef TestNA49Histo_h
#define TestNA49Histo_h

#include "TstHistoSet.hh"

#include <vector>

#include "TProfile.h"

class TestNA49Histo : public TstHistoSet
{

public:

   TestNA49Histo( G4String );
   ~TestNA49Histo();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:
   
   TH1D*              fHistoNSec;
   std::vector<TH1D*> fHistoSecProton; 
   std::vector<TH1D*> fHistoSecAntiProton;
   std::vector<TH1D*> fHistoSecPiMinus; 
   std::vector<TH1D*> fHistoSecPiPlus; 
   std::vector<TH1D*> fHistoSecNeutron;

};

#endif
