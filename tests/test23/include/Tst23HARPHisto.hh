#ifndef Tst23HARPHisto_h
#define Tst23HARPHisto_h

#include "TstHistoSet.hh"

#include <vector>

// #include "TProfile.h"
#include "TH1.h"

class Tst23HARPHisto : public TstHistoSet
{

public:

   Tst23HARPHisto( G4String );
   ~Tst23HARPHisto();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double xsec=1. );


private:
   
   TH1D*              fHistoNSec;
   
   std::vector<TH1D*> fHistoSecPiMinusFW; 
   std::vector<TH1D*> fHistoSecPiPlusFW; 
   std::vector<TH1D*> fHistoSecPiMinusLA; 
   std::vector<TH1D*> fHistoSecPiPlusLA; 
      
//   std::vector<TH1D*> fHistoSecPiMinusTot; 
//   std::vector<TH1D*> fHistoSecPiPlusTot; 
   
   G4int                fNThetaBinsFW;
   G4double             fThetaMinFW;
   G4double             fDeltaThetaFW;   
   G4int                fNThetaBinsLA;
   G4double             fThetaMinLA;
   G4double             fDeltaThetaLA;   

};

#endif
