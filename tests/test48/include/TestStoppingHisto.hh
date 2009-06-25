#ifndef TestStoppingHisto_h
#define TestStoppingHisto_n

#include <string>
#include <vector>

#include "TFile.h"
#include "TH1F.h"

// fwd declaration
class G4VParticleChange;
class G4DynamicParticle;

class TestStoppingHisto
{

   public:
      
      // ctor & dtor
      //
      TestStoppingHisto( std::string beam, std::string target, std::string model) 
         : fBeam(beam), fTarget(target), fModel(model) { Init(); }
      ~TestStoppingHisto();
      
      void FillEvt( G4VParticleChange* );
      void Write( int ) ;

   private:
      
      void Init();
      void InitPionMinus();
      void InitAntiProton();
      
      // void FillEvtAntiProton( const G4DynamicParticle* );
      
      // data members
      //
      std::string        fBeam;
      std::string        fTarget;
      std::string        fModel;
      std::vector<TH1F*> fHisto; 

};

#endif
