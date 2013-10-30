#ifndef Test19Histo_h
#define Test19Histo_h

//#include <string>
//#include "TFile.h"

//#include "G4LorentzVector.hh"

#include "TstHisto.hh"

// fwd declaration
class G4VParticleChange;
class TstReader;

class Test19Histo : public TstHisto
{

public:
      
   // ctor & dtor
   //
   Test19Histo( const TstReader* pset ) ;
   virtual ~Test19Histo() {} ;
   
   virtual TFile* OpenHistoFile();

private:
 
    G4String fHistoDirName;

};

#endif
