#ifndef Tst23NA49Histo_h
#define Tst23NA49Histo_h

#include "TstHistoSet.hh"
#include "TstHistoSetForNu.hh"

#include <vector>
#include <map>

#include "TProfile.h"
#include "TH2.h"

class Tst23NA49Histo : public TstHistoSet, public TstHistoSetForNu
{

public:

   Tst23NA49Histo( G4String );
   ~Tst23NA49Histo();
   
   void FillEvt( G4VParticleChange*, const G4LorentzVector&, const G4LorentzVector& );

protected:

   virtual void Write( G4int stat=1, G4double scale=1. );


private:
   
   double CalculateBinWeight( const G4LorentzVector&, double, double, double, int, double );
   
   TH1D*              fHistoNSec;
   
   std::vector<TH1D*> fHistoSecProton; 
   std::vector<TH1D*> fHistoSecAntiProton;
   std::vector<TH1D*> fHistoSecPiMinus; 
   std::vector<TH1D*> fHistoSecPiPlus; 
   std::vector<TH1D*> fHistoSecNeutron;
   
   TH2D* fHistoPTvsXFProton;
   TH2D* fHistoPTvsXFAntiProton;
   TH2D* fHistoPTvsXFPiMinus;
   TH2D* fHistoPTvsXFPiPlus;
//    TH2D* fHistoPTvsXFNeutron;
   
   std::vector<TH1D*> fHistoPTPiMinus;
   std::vector<TH1D*> fHistoPTPiPlus;
   
   int                fNPiBinsXF;
   double*            fPiBinsXF;
   int                fNPiBinsPT;
   double*            fPiBinsPT;
      
//NOTE: this needs to be refined
//      in general, one can NOT merge it all together with the outcome
//      of the 1st interaction, because different secondaries of different
//      energies contribte to re-iteractions, thus there's no way to check
//      vs NA49 data that are specifically for 158GeV p on C
//      
   std::vector<TH1D*> fHistoSecProtonTot; 
   std::vector<TH1D*> fHistoSecAntiProtonTot;
   std::vector<TH1D*> fHistoSecPiMinusTot; 
   std::vector<TH1D*> fHistoSecPiPlusTot; 
   std::vector<TH1D*> fHistoSecNeutronTot;

// plots for differnt E-ranges of the final nu
//
   std::map< NuERange, std::vector<TH1D*> > fHistoPiMinusFrom1stInt;
   std::map< NuERange, std::vector<TH1D*> > fHistoPiPlusFrom1stInt;

   std::map< InteractionType, std::map< NuERange, std::vector<TH1D*> > > fHistoPiMinusFromReint;
   std::map< InteractionType, std::map< NuERange, std::vector<TH1D*> > > fHistoPiPlusFromReint;
   
};

#endif
