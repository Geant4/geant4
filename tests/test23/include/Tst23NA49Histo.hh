#ifndef Tst23NA49Histo_h
#define Tst23NA49Histo_h

#include "TstHistoSet.hh"

#include <vector>

#include "TProfile.h"

class Tst23NA49Histo : public TstHistoSet
{

public:

   Tst23NA49Histo( G4String );
   ~Tst23NA49Histo();
   
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
      
   std::vector<TH1D*> fHistoSecProtonTot; 
   std::vector<TH1D*> fHistoSecAntiProtonTot;
   std::vector<TH1D*> fHistoSecPiMinusTot; 
   std::vector<TH1D*> fHistoSecPiPlusTot; 
   std::vector<TH1D*> fHistoSecNeutronTot;

};

#endif
