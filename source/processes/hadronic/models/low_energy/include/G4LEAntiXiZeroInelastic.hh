// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEAntiXiZeroInelastic.hh,v 1.1 1999-01-07 16:12:39 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy AntiXiZero Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Feb-1997
 // Last modified: 27-Mar-1997
 
#ifndef G4LEAntiXiZeroInelastic_h
#define G4LEAntiXiZeroInelastic_h 1
 
#include "G4InelasticInteraction.hh"
 
 class G4LEAntiXiZeroInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEAntiXiZeroInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LEAntiXiZeroInelastic() { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void
     Cascade(                               // derived from CASAL0
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
 
