 // Hadronic Process: Low Energy Neutron Inelastic Process
 // original by J.L. Chuma, TRIUMF, 04-Feb-1997
 
#ifndef G4LENeutronInelastic_h
#define G4LENeutronInelastic_h 1
 
// Class Description
// Final state production model for Neutron inelastic scattering below 20 GeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

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
 
