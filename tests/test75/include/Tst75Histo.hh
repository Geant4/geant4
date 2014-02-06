#ifndef Tst75Histo_h
#define Tst75Histo_h

#include "TstHisto.hh"

// fwd declaration
//class G4VParticleChange;
class TstReader;

class Tst75Histo : public TstHisto
{

public:
      
   // ctor & dtor
   //
   Tst75Histo( const TstReader* pset ) ;
   virtual ~Tst75Histo() {} ;
   
   virtual TFile* OpenHistoFile();

//private:
 
//    G4String fHistoDirName;

};

#endif
