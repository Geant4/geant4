// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4HEAntiXiMinusInelastic.hh,v 1.3 1999-12-15 14:52:52 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// G4 Gheisha High Energy model class -- header file
// H. Fesefeldt, RWTH Aachen 23-October-1996
// Last modified: 10-December-1996

// A prototype of the Gheisha High Energy collision model.

#ifndef G4HEAntiXiMinusInelastic_h
#define G4HEAntiXiMinusInelastic_h 1

#include "G4HEInelastic.hh"

class G4HEAntiXiMinusInelastic : public G4HEInelastic  
{
 public: 
        G4HEAntiXiMinusInelastic() : G4HEInelastic()
           {
             G4int    vecLen       = 0;
           }

        ~G4HEAntiXiMinusInelastic(){ };
         
        G4int vecLength;
        
        G4VParticleChange * ApplyYourself( const G4Track &aTrack, G4Nucleus &targetNucleus );

        G4int  GetNumberOfSecondaries()
               { return vecLength; }         

 private:

        void   FirstIntInCasAntiXiMinus(G4bool &inElastic, const G4double availableEnergy,
                                        G4HEVector pv[],
                                        G4int &vecLen, 
                                        G4HEVector incidentParticle,
                                        G4HEVector targetParticle,
                                        const G4double atomicWeight);
};
#endif                     
                                         

