 // Hadronic Process: Low Energy Neutron Inelastic Process
 // original by J.L. Chuma, TRIUMF, 04-Feb-1997
 
#ifndef G4LENeutronInelastic_h
#define G4LENeutronInelastic_h 1
 
#include "G4InelasticInteraction.hh"
 
 class G4LENeutronInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LENeutronInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LENeutronInelastic()
    { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void Cascade(                               // derived from CASN
      G4FastVector<G4ReactionProduct,128> &vec,
      G4int &vecLen,
      const G4DynamicParticle *originalIncident,
      G4ReactionProduct &currentParticle,
      G4ReactionProduct &targetParticle,
      G4bool &incidentHasChanged, 
      G4bool &targetHasChanged,
      G4bool &quasiElastic );
    
    void SlowNeutron(
     const G4DynamicParticle *originalIncident,
     G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &targetParticle,
     G4Nucleus & targetNucleus );
 };
 
#endif
 
