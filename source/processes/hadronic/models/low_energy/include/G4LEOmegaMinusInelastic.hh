// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEOmegaMinusInelastic.hh,v 1.2 1999-12-15 14:53:06 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy OmegaMinus Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Nov-1996
 // Last modified: 27-Mar-1997
 
#ifndef G4LEOmegaMinusInelastic_h
#define G4LEOmegaMinusInelastic_h 1
 
#include "G4InelasticInteraction.hh"
 
 class G4LEOmegaMinusInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEOmegaMinusInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LEOmegaMinusInelastic() { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void Cascade(                               // derived from CASOM
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
 
