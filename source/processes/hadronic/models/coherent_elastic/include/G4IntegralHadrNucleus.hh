// G4IntegralHadrNucleus.hh

#ifndef  G4IntegralHadrNucleus_h
#define  G4IntegralHadrNucleus_h 1

#include "globals.hh"
//#include "G4HadronicProcess.hh"
//#include "G4Isotope.hh"
//#include "G4VDiscreteProcess.hh"
#include "G4DynamicParticle.hh"
#include "G4Element.hh"

class G4IntegralHadrNucleus 

 {

   public:

       G4double GetElasticCrossSection(const G4DynamicParticle*  aHadron,
                                             G4int               iNucleus); 

       G4double GetTotalCrossSection(const G4DynamicParticle*  aHadron,
                                           G4int               iNucleus); 

       G4double GetProductionCrossSection(const G4DynamicParticle*  aHadron,
                                                G4int               iNucleus); 

       G4double GetInelasticCrossSection(const G4DynamicParticle*  aHadron, 
                                               G4int               iNucleus); 

       G4double GetQuasyElasticCrossSection(const G4DynamicParticle*  aHadron,             
                                                  G4int               iNucleus); 

   private:  

       void GetIntegralCrSec(G4int A);

       void GetHadronValues(G4int       iHadron,
                            G4double    HadrEnergy);

       void TransformValues(const G4DynamicParticle* aHadron,
                                  G4int              aNucleus);

       G4int     Anucleus, iHadron;

       G4String  TypeHadron;
          
       G4double  TotalCrSec,  InelCrSec,   ProdCrSec,   
                 ElasticCrSec, QuasyElasticCrSec;

       G4double  HadrTot, HadrSlope, HadrReIm, DDSect1, DDSect2, DDSect3, 
                 sHadr, mHadr, HadrEnergy;

 static const G4String  HadronNames[6];
 };

#endif

 
