// G4IntegralHadrNucleus.hh

//#ifndef  G4IntegralHadrNucleus_h
//#define  G4IntegralHadrNucleus_h 1

//#include "globals.hh"
//#include "G4HadronicProcess.hh"
//#include "G4Isotope.hh"
//#include "G4VDiscreteProcess.hh"
//#include "G4DynamicParticle"

class G4IntegralHadrNucleus 

 {

   public:

//        G4IntegralHadrNucleus(const G4DynamicParticle*  aHadron,
//                              const G4Isotope*          aNucleus)
          G4IntegralHadrNucleus(char*   TypeHadron, double HadrEnergy,
                                int     aNucleus)
          {
//             G4ParticleDefinition*  aParticle;

//              aParticle  = aHadromn->GetDefinition();
//              TypeHadron = aParticle->GetParticleName();
//              mHadr      = aHadron->GetMass();  
//              HadrEnergy = aHadron->GetTotalEnergy();
            if( HadrEnergy<5)  cout << " This method is not applicabale at that energy ";

            if( TypeHadron!="proton" || TypeHadron!="anti_proton" ||
                TypeHadron!="pi+"    || TypeHadron!="pi-"         ||
                TypeHadron!="kaon+"  || TypeHadron!="kaon-") 
              cout << " This hadron is absent in this model ";

//                Anucleus   = aNucleus->GetN();
          }

       void   GetIntegralCrSec(int       Anucleus);

       double  TotalCrSec,  InelCrSec,   ProdCrSec,   
                ElasticCrSec, QuasyElasticCrSec;

//   private:

       void   GetHadronValues(char*       TypeHadron,
                              double       HadrEnergy);

//      G4int    Anucleus;
//      G4String TypeHadron;
//      G4double HadrTot, HadrSlope, HadrReIm, DDSect1, DDSect2, sHadr, 
//               mHadr, HadrEnergy;
          
//        int     Anucleus;
//        char*   TypeHadron;
        double  HadrTot,  HadrSlope, HadrReIm, DDSect1, DDSect2, HadrFnergy,
                 mHadr, sHadr, DDSect3;
        char*   HadronsName[6]={"proton","anti_proton","pi+","pi-",
                                "kaon+","kaon-");
 }

//#endif

 
