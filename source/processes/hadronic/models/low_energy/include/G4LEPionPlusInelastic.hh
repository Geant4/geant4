// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEPionPlusInelastic.hh,v 1.1 1999-01-07 16:12:41 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy PionPlus Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Nov-1996
 // Last modified: 27-Mar-1997
 
#ifndef G4LEPionPlusInelastic_h
#define G4LEPionPlusInelastic_h 1
 
#include "G4InelasticInteraction.hh"
 
 class G4LEPionPlusInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEPionPlusInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LEPionPlusInelastic() { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void Cascade(                               // derived from CASPIP
      G4FastVector<G4ReactionProduct,128> &vec,
      G4int &vecLen,
      const G4DynamicParticle *originalIncident,
      G4ReactionProduct &currentParticle,
      G4ReactionProduct &targetParticle,
      G4bool &incidentHasChanged, 
      G4bool &targetHasChanged,
      G4bool &quasiElastic );
    
 };
#endif
 
