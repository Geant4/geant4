// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4HEAntiSigmaZeroInelastic_h
#define G4HEAntiSigmaZeroInelastic_h 1

#include "G4HEAntiLambdaInelastic.hh"

class G4HEAntiSigmaZeroInelastic : public G4HEInelastic  
{
 public: 
        G4HEAntiSigmaZeroInelastic() : G4HEInelastic()
           {
           }

        ~G4HEAntiSigmaZeroInelastic(){ };
         
        G4int verboseLevel;
        G4int MAXPART;
        G4int vecLength;

        void SetMaxNumberOfSecondaries(G4int maxnumber)
           { MAXPART = maxnumber; };
        void SetVerboseLevel(G4int verbose)
           { verboseLevel = verbose;};
        G4VParticleChange * ApplyYourself( const G4Track &aTrack, G4Nucleus &targetNucleus );
        G4int  GetNumberOfSecondaries()
               { return vecLength; }         
};
#endif                     
                                         

