// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4LEDeuteronInelastic.cc,v 1.2 1999-12-15 14:53:08 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
 // Hadronic Process: Deuteron Inelastic Process
 // J.L. Chuma, TRIUMF, 25-Feb-1997
 // Last modified: 27-Mar-1997

#include "G4LEDeuteronInelastic.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

 G4VParticleChange *
  G4LEDeuteronInelastic::ApplyYourself( const G4Track &aTrack,
                                        G4Nucleus &targetNucleus )
  {
    theParticleChange.Initialize( aTrack );
    
    const G4DynamicParticle *originalIncident = aTrack.GetDynamicParticle();
    
    if( verboseLevel > 1 )
    {
      G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEDeuteronInelastic::ApplyYourself called" << G4endl;
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
    }
    
    // Work-around for lack of model above 100 MeV
    if (originalIncident->GetKineticEnergy()/MeV > 100. ||
        originalIncident->GetKineticEnergy() <= 0.1*MeV) return &theParticleChange;

    G4double N = targetNucleus.GetN();
    G4double Z = targetNucleus.GetZ();
    G4double theAtomicMass = targetNucleus.AtomicMass( N, Z )-(Z)*G4Electron::Electron()->GetPDGMass();
    G4double massVec[9];
    massVec[0] = targetNucleus.AtomicMass( N+2.0, Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[1] = targetNucleus.AtomicMass( N+1.0, Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[2] = targetNucleus.AtomicMass( N+1.0, Z     )-(Z)*G4Electron::Electron()->GetPDGMass();
    massVec[3] = theAtomicMass;
    massVec[4] = targetNucleus.AtomicMass( N-1.0, Z     )-(Z)*G4Electron::Electron()->GetPDGMass();
    massVec[5] = targetNucleus.AtomicMass( N-2.0, Z-1.0 )-(Z-1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[6] = targetNucleus.AtomicMass( N,     Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[7] = massVec[3];
    massVec[8] = targetNucleus.AtomicMass( N,     Z-1.0 )-(Z-1.0)*G4Electron::Electron()->GetPDGMass();
    
    G4FastVector<G4ReactionProduct,3> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    theReactionDynamics.NuclearReaction( vec, vecLen, originalIncident,
                                         targetNucleus, theAtomicMass, massVec );
    
    G4double p = originalIncident->GetTotalMomentum();
    theParticleChange.SetMomentumChange( originalIncident->GetMomentum() * (1.0/p) );
    theParticleChange.SetEnergyChange( originalIncident->GetKineticEnergy() );
    
    theParticleChange.SetNumberOfSecondaries( vecLen );
    G4DynamicParticle *pd;
    for( G4int i=0; i<vecLen; ++i )
    {
      pd = new G4DynamicParticle();
      pd->SetDefinition( vec[i]->GetDefinition() );
      pd->SetMomentum( vec[i]->GetMomentum() );
      theParticleChange.AddSecondary( pd );
    }
    return &theParticleChange;
  }
 
 /* end of file */
 
