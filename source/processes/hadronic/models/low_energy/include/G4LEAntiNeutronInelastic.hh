// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEAntiNeutronInelastic.hh,v 1.3 2000-12-14 09:12:43 hpw Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Low Energy AntiNeutron Inelastic Process
 // J.L. Chuma, TRIUMF, 19-Feb-1997
 // Last modified: 27-Mar-1997
 
#ifndef G4LEAntiNeutronInelastic_h
#define G4LEAntiNeutronInelastic_h 1
 
// Class Description
// Final state production model for AntiNeutron inelastic scattering below 20 GeV; 
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Class Description - End

#include "G4InelasticInteraction.hh"
 
 class G4LEAntiNeutronInelastic : public G4InelasticInteraction
 {
 public:
    
    G4LEAntiNeutronInelastic() : G4InelasticInteraction()
    {
      SetMinEnergy( 0.0 );
      SetMaxEnergy( 25.*GeV );
    }
    
    ~G4LEAntiNeutronInelastic()
    { }
    
    G4VParticleChange *ApplyYourself( const G4Track &aTrack,
                                      G4Nucleus &targetNucleus );
    
 private:
    
    void
     Cascade(                               // derived from CASNB
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
 
