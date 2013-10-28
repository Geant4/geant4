#ifndef Tst23Histo_h
#define Tst23Histo_h

#include "TstHisto.hh"

// fwd declaration
class G4VParticleChange;
class TstReader;

class Tst23Histo : public TstHisto
{

public:
      
   // ctor & dtor
   //
   Tst23Histo( const TstReader* pset ) ;
   virtual ~Tst23Histo() {} ;
   
   virtual TFile* OpenHistoFile();
   
private:

   G4String fHistoDirName;

};

#endif
