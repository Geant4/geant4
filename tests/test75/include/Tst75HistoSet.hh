#ifndef Tst75HistoSet_h
#define Tst75HistoSet_h

#include "TstHistoSet.hh"

#include <vector>

#include "TProfile.h"

// class G4VParticleChange;

class Tst75HistoSet : public TstHistoSet
{

public:

   Tst75HistoSet( G4String, G4float );
   ~Tst75HistoSet();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );

private:
   
   // data members
   //
   // plots for secondary proton
   //
   TH1F* fProton45;
   TH1F* fProton60;
   TH1F* fProton72;
   TH1F* fProton84;
   TH1F* fProton90;
   TH1F* fProton135;
   TH1F* fProKE;
   TH1F* fProCosTh;
   //
   // plots for secondary pi-
   //
   TH1F* fPim28deg;
   TH1F* fPim44deg;
   //
   // plots for secondary pi+
   //
   TH1F* fPip28deg;
   TH1F* fPip44deg;
   
   G4float fMaxKE;
   
   // G4VParticleChange* fInteraction;

};

#endif
