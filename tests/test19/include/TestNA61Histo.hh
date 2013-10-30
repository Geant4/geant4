#ifndef TestNA61Histo_h
#define TestNA61Histo_h

#include "TstHistoSet.hh"

#include <vector>

#include "TH1F.h"

class TestNA61Histo : public TstHistoSet
{

public:

   TestNA61Histo( G4String );
   ~TestNA61Histo();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:
   
   std::vector<TH1F*> fHistoSecProton; 
   std::vector<TH1F*> fHistoSecPiMinus; 
   std::vector<TH1F*> fHistoSecPiPlus;
   std::vector<TH1F*> fHistoSecKPlus;
   std::vector<TH1F*> fHistoSecPiPlus2; 

};

#endif
