// G4ElasticHadrNucleusHe.hh

#ifndef G4ElasticHadrNucleusHE_h
#define G4ElasticHadrNucleusHE_h 1

#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "Randomize.hh"
#include "G4Nucleus.hh"
#include "G4IonTable.hh"
#include "G4DiffElasticHadrNucleus.hh"
#include "G4IntegrHadrNucleus.hh"
#include "G4HadronicInteraction.hh"
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"
#include "iostream.h"
#include "fstream.h"

#define   ONQ2     150      //  The number of steps on Q2
#define   ONE      5        //  The number of steps on E
#define   AreaNumb 6        //  The number of decimal areas on E
#define   ONQ2XE   ONQ2*ONE //  The dimension of a distr. func. array
#define   MaxN     10       //  The atomic number where the calculation
                            //  on the formula is changed on the integral
                            //  one

   class G4ElasticHadrNucleusHE : public G4DiffElasticHadrNucleus,
                                         G4HadronicInteraction 
   {

 public:
         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus);

         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4double             dEbeg,
                                      G4double             dEend,
                                      G4int                iNpoE,
                                      G4String          sNameFile);

         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4String          sNameFile);  

         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4double             dEbeg,
                                      G4double             dEend,
                                      G4int                iNpoE);

         G4ElasticHadrNucleusHE(const G4DynamicParticle * aHadron,
                                      G4Nucleus         * aNucleus,
                                      G4int                 iNpoE);  

        ~G4ElasticHadrNucleusHE() {;}

         G4VParticleChange * ApplyYourself( const G4Track   &aTrack,
                                                  G4Nucleus &aNucleus);

         G4double RandomElastic0( const G4DynamicParticle *   aHadron,
                                        G4Nucleus *           aNucleus);

         G4double RandomElastic1( const G4DynamicParticle *   aHadron,
                                        G4Nucleus *           aNucleus);
 private:
         G4String GetHadronName(const G4DynamicParticle * aHadron);

         G4double GetQ2limit(G4double R1);

         void  CreationArray(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         void  ArrayForHeavy(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         void  ArrayForLight(const G4DynamicParticle * aHadron,
                                   G4Nucleus         * aNucleus);

         G4double InterPol(G4double X1, G4double X2, G4double X3,
                           G4double Y1, G4double Y2, G4double Y3, 
                           G4double X);

//    +++++++++++++++++++++++++++++++++++++++++++++++++++++
   G4double Factorial1(G4int N)
     {
        G4double  Res;
                  Res = 1;
              if(N == 0) return Res;

          if(N < 100) for(G4int M = 1; M<=N; M++)  Res = Res*M;         

           else  Res = 2.50662827*exp(-N-1)*pow(N+1,N+0.5)*
                         (1+1/12/(N+1)+1/288/(N+1)/(N+1)-
                         139/51840/(N+1)/(N+1)/(N+1)-
                         571/2488320/(N+1)/(N+1)/(N+1)/(N+1));
              return Res;
      }

//     ++++++++++++++++++++++++++++++++++++++++++++++++++
         G4IonTable                * MyIonTable;
         G4DiffElasticHadrNucleus    aDiffElHadNcls;
 
         G4int     Nstep,         // The number of steps on Q2
                   iKindWork,     // 
                   iContr,        //
                   iPoE;          // The number of steps on E
         G4int     iTypeWork;
         G4double  aNucleon, 
                 * pTableCrSec,   //  The array of distr. func.
                                  //  at all energies
                 * pTableE;       //  The array of E values
         G4double  iQ2[ONQ2],     //  The array of Q2 values
                   iIntgr[ONQ2],  //  The array of distr. func.   
                                  //    at one energy
                   Factorials1[250]; // The array for factorials
         G4double  dEbeg1, dEend1, dQ2, maxQ2;

       };     //   The end of the class description
  
#endif
