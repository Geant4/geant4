// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4InelasticInteraction.hh,v 1.2 1999-12-15 14:52:07 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Inelastic Interaction 
 // This class is an abstract base class, since the pure virtual
 // function ApplyYourself has not been defined yet.
 // original by H.P. Wellisch
 // Modified by J.L. Chuma, TRIUMF, 22-Nov-1996
 // Modified by J.L. Chuma  27-Mar-1997
 // Modified by J.L. Chuma  30-Apr-1997
 // Modified by J.L. Chuma  05-Aug-1997  to pass the original incident particle to
 //                                      CalculateMomenta
 // Modified by J.L. Chuma  05-Jun-1998  to include quasiElastic flag to allow for
 //                                      TwoBody to be called directly, bypassing
 //                                      TwoCluster, and allowing TwoCluster to be
 //                                      called with no secondaries
 
#ifndef G4InelasticInteraction_h
#define G4InelasticInteraction_h 1

#include "globals.hh"
#include "G4FastVector.hh"
#include "G4HadronicInteraction.hh"
#include "G4ReactionProduct.hh"
#include "G4ParticleTypes.hh" 
#include "Randomize.hh"
 
 class G4InelasticInteraction : public G4HadronicInteraction
 {
 public:
    
    G4InelasticInteraction() : G4HadronicInteraction()
    { }
    
    virtual ~G4InelasticInteraction()
    { }
    
 protected:
    
    G4double Pmltpc( G4int np, G4int nm, G4int nz, G4int n,
                     G4double b, G4double c );
    
    G4bool MarkLeadingStrangeParticle( const G4ReactionProduct &currentParticle,
                                       const G4ReactionProduct &targetParticle,
                                       G4ReactionProduct &leadParticle );
    
    void SetUpPions( const G4int np, const G4int nm, const G4int nz,
                     G4FastVector<G4ReactionProduct,128> &vec,
                     G4int &vecLen );
    
    void GetNormalizationConstant( const G4double availableEnergy,
                                   G4double &n,
                                   G4double &anpn );
    
    void CalculateMomenta( G4FastVector<G4ReactionProduct,128> &vec,
                           G4int &vecLen,
                           const G4DynamicParticle *originalIncident,
                           const G4DynamicParticle *originalTarget,
                           G4ReactionProduct &modifiedOriginal,
                           G4Nucleus &targetNucleus,
                           G4ReactionProduct &currentParticle,
                           G4ReactionProduct &targetParticle,
                           G4bool &incidentHasChanged,
                           G4bool &targetHasChanged,
                           G4bool quasiElastic );
    
    void SetUpChange( G4FastVector<G4ReactionProduct,128> &vec,
                      G4int &vecLen,
                      G4ReactionProduct &currentParticle,
                      G4ReactionProduct &targetParticle,
                      G4bool &incidentHasChanged );
 };
 
#endif
 
