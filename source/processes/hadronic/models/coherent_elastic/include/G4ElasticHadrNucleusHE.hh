//      G4ElasticHadrNucleusHe.hh

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
#include <rw/tpordvec.h>
#include <rw/cstring.h>
#include "globals.hh"
#include "G4ParticleDefinition.hh"
#include "G4Ions.hh"
#include "G4ParticleTable.hh"
#include "G4NucleiProperties.hh"

//#include "G4Integrator.hh"

   class G4ElasticHadrNucleusHE : public G4DiffElasticHadrNucleus,
                                         G4HadronicInteraction 
   {

 public:
         G4ElasticHadrNucleusHE() : G4DiffElasticHadrNucleus(),
                                    G4HadronicInteraction()
          {
            for(G4int ii=1; ii<149; ii++)
              {
                coefSimp[2*ii]   = 4;
                coefSimp[2*ii+1] = 2;
              }
                coefSimp[0]      = 1;
          }

        ~G4ElasticHadrNucleusHE() {;}

         G4ParticleChange * ApplyYourself( const G4Track &aTrack,
                                           G4Nucleus     &aNucleus);

 private:
         G4double RandomElastic( const G4DynamicParticle *   aHadron,
                                       G4Nucleus *           aNucleus);

//         G4IonTable                 aIonTable;
//         G4IntegrHadrNucleus         aIntegralHadNcls;
         G4DiffElasticHadrNucleus    aDiffElHadNcls; 
         G4int Nstep;
//         G4Integrator integral;
         G4double aNucleon;  
         G4double coefSimp[300], iQ2[200], iIntgr[200];
   };

#endif
