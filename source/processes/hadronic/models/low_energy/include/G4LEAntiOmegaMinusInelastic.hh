// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEAntiOmegaMinusInelastic.hh,v 1.1 1999-01-07 16:12:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy AntiOmegaMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Feb-1997
 // Last modified: 27-Mar-1997
 
#ifndef G4LEAntiOmegaMinusInelastic_h
#define G4LEAntiOmegaMinusInelastic_h 1
 
#include "G4InelasticInteraction.hh"
 
 class G4LEAntiOmegaMinusInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEAntiOmegaMinusInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LEAntiOmegaMinusInelastic() { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void
     Cascade(                               // derived from CASOM
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
 
