// G4HadronValues.hh

#ifndef  G4HadronValues_h
#define  G4HadronValues_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Proton.hh"
#include "G4Neutron.hh"

#include "G4AntiNeutron.hh"
#include "G4AntiProton.hh"
#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

class G4HadronValues 
 {
   public:
      
       G4HadronValues() {;}
      ~G4HadronValues() {;}

   protected:  

       void GetHadronValues(const G4DynamicParticle * aHadron);

       G4double  HadrTot, HadrSlope, HadrReIm,  DDSect2, DDSect3,
                 MomentumCM;

  };

#endif

 
