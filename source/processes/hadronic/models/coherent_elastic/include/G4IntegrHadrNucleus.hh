// G4IntegrHadrNucleus.hh

#ifndef  G4IntegrHadrNucleus_h
#define  G4IntegrHadrNucleus_h 1

#include "globals.hh"
#include "G4DynamicParticle.hh"
#include "G4Nucleus.hh"
#include "G4HadronValues.hh"

         class G4IntegrHadrNucleus : public G4HadronValues
 {

   public:

       G4IntegrHadrNucleus() : G4HadronValues() {;}
      ~G4IntegrHadrNucleus() {;}

       G4double GetElasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

       G4double GetTotalCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

       G4double GetProductionCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 


       G4double GetInelasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 


       G4double GetQuasyElasticCrossSection(
                              const G4DynamicParticle *  aHadron,
                                    G4Nucleus         * aNucleus); 

   private:  

       void GetIntegralCrSec(G4Nucleus * aNucleus);

       G4double  TotalCrSec,  InelCrSec,  ProdCrSec,  ElasticCrSec, 
                 QuasyElasticCrSec,  HadrEnergy; 

 };

#endif

 
