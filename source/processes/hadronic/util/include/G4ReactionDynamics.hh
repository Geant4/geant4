// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4ReactionDynamics.hh,v 1.1 1999-01-07 16:13:49 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Reaction Dynamics
 // original by H.P. Wellisch
 // Modified by J.L.Chuma 19-Nov-96
 // Modified by J.L.Chuma 27-Mar-97
 // Modified by J.L.Chuma 30-Apr-97
 // Modified by J.L.Chuma 06-Aug-97  to include the original incident particle
 //                                  before Fermi motion and evaporation effects
 
#ifndef G4ReactionDynamics_h
#define G4ReactionDynamics_h 1

#include "G4ParticleTypes.hh"
#include "G4DynamicParticle.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4FastVector.hh"
 
 class G4ReactionDynamics 
 {
 public:
    
    G4ReactionDynamics() {}
    
    ~G4ReactionDynamics() {}
    
    virtual G4double FindInelasticity()
    { return 0.0; }
    
    virtual G4double FindTimeDelay()
    { return 0.0; }
    
    G4bool GenerateXandPt(                   // derived from GENXPT
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen,
     G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4DynamicParticle *originalIncident,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged, 
     G4bool &targetHasChanged,
     G4bool leadFlag,
     G4ReactionProduct &leadingStrangeParticle );
    
    void SuppressChargedPions(
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen,
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged,
     G4bool &targetHasChanged );
      
    G4bool TwoCluster(                       // derived from TWOCLU
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen,
     G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4DynamicParticle *originalIncident,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged, 
     G4bool &targetHasChanged,
     G4bool leadFlag,
     G4ReactionProduct &leadingStrangeParticle );
    
    void TwoBody(                         // derived from TWOB
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen,
     G4ReactionProduct &modifiedOriginal,
     const G4DynamicParticle *originalTarget,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4Nucleus &targetNucleus,
     G4bool &targetHasChanged );
    
    G4int Factorial( G4int n );
    
    G4double GenerateNBodyEvent(                // derived from PHASP
     const G4double totalEnergy,
     const G4bool constantCrossSection,
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen );
    
    void ProduceStrangeParticlePairs(
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen,
     const G4ReactionProduct &modifiedOriginal,
     const G4DynamicParticle *originalTarget, 
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4bool &incidentHasChanged,
     G4bool &targetHasChanged );
    
    void NuclearReaction(                     // derived from NUCREC
     G4FastVector<G4ReactionProduct,3> &vec,
     G4int &vecLen,
     const G4DynamicParticle *originalIncident,
     const G4Nucleus &aNucleus,
     const G4double theAtomicMass,
     const G4double *massVec );
    
 private:
    
    void Rotate(
     const G4double numberofFinalStateNucleons,
     const G4ThreeVector &temp,
     const G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4DynamicParticle *originalIncident,
     const G4Nucleus &targetNucleus,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen );
    
    void Defs1(
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen );
    
    void AddBlackTrackParticles(
     const G4double epnb,
     const G4int npnb,
     const G4double edta,
     const G4int ndta,
     const G4double sprob,
     const G4double kineticMinimum,
     const G4double kineticFactor,
     const G4ReactionProduct &modifiedOriginal,
     G4double spall,
     const G4Nucleus &aNucleus,
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen );
    
    void MomentumCheck(
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,128> &vec,
     G4int &vecLen );
    
    G4double normal();
    
    G4int Poisson( G4double x );
   
 };
 
#endif
 
