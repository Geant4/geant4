// G4HadronValues.hh

#ifndef  G4HadronValues_h
#define  G4HadronValues_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4Lambda.hh"
#include "G4SigmaPlus.hh"
#include "G4SigmaMinus.hh"
#include "G4SigmaZero.hh"
#include "G4XiMinus.hh"
#include "G4XiZero.hh"
#include "G4OmegaMinus.hh"

#include "G4AntiNeutron.hh"
#include "G4AntiProton.hh"
#include "G4AntiLambda.hh"
#include "G4AntiSigmaPlus.hh"
#include "G4AntiSigmaMinus.hh"
#include "G4AntiSigmaZero.hh"
#include "G4AntiXiMinus.hh"
#include "G4AntiXiZero.hh"
#include "G4AntiOmegaMinus.hh"

#include "G4PionPlus.hh"
#include "G4PionMinus.hh"
#include "G4KaonPlus.hh"
#include "G4KaonMinus.hh"

#define MyPi      3.141593
#define MbToGeV2  2.568
#define GeV2ToMb  0.38939
#define MbToFm2   25.68

class G4HadronValues 
 {
   public:
      
       G4HadronValues() {;}
      ~G4HadronValues() {;}

//   protected:  

       void GetHadronValues(const G4DynamicParticle * aHadron);

       G4double  HadrTot, HadrSlope, HadrReIm,  DDSect2, DDSect3,
                 MomentumCM;

  };

#endif

 
