// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEAntiSigmaZeroInelastic.hh,v 1.3 1999-12-15 14:52:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
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
                                         

