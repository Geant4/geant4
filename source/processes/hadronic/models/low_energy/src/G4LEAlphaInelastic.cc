//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//
 // Hadronic Process: Alpha Inelastic Process
 // J.L. Chuma, TRIUMF, 25-Feb-1997
 // Last modified: 27-Mar-1997
 // J.L. Chuma, 08-May-2001: Update original incident passed back in vec[0]
 //                          from NuclearReaction
 //
#include "G4LEAlphaInelastic.hh"
#include "Randomize.hh"
#include "G4Electron.hh"
 
 G4VParticleChange *
  G4LEAlphaInelastic::ApplyYourself( const G4Track &aTrack,
                                     G4Nucleus &targetNucleus )
  {
    theParticleChange.Initialize( aTrack );
    G4double A = targetNucleus.GetN();
    G4double Z = targetNucleus.GetZ();
    
    const G4DynamicParticle *originalIncident = aTrack.GetDynamicParticle();
    
    if( verboseLevel > 1 )
    {
      G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEAlphaInelastic::ApplyYourself called" << G4endl;
      G4cout << "kinetc energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << G4endl;
    }
    
    // Work-around for lack of model above 100 MeV
    if (originalIncident->GetKineticEnergy()/MeV > 100. ||
        originalIncident->GetKineticEnergy() <= 0.1*MeV) return &theParticleChange;

    G4double theAtomicMass = targetNucleus.AtomicMass( A, Z )-Z*G4Electron::Electron()->GetPDGMass();
    G4double massVec[9];
    massVec[0] = targetNucleus.AtomicMass( A+4.0, Z+2.0 )-(Z+2.0)*G4Electron::Electron()->GetPDGMass();
    massVec[1] = targetNucleus.AtomicMass( A+3.0, Z+2.0 )-(Z+2.0)*G4Electron::Electron()->GetPDGMass();
    massVec[2] = targetNucleus.AtomicMass( A+3.0, Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[3] = targetNucleus.AtomicMass( A+2.0, Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[4] = targetNucleus.AtomicMass( A+1.0, Z+1.0 )-(Z+1.0)*G4Electron::Electron()->GetPDGMass();
    massVec[5] = theAtomicMass;
    massVec[6] = targetNucleus.AtomicMass( A+2.0, Z+2.0 )-(Z+2.0)*G4Electron::Electron()->GetPDGMass();
    massVec[7] = massVec[3];
    massVec[8] = targetNucleus.AtomicMass( A+2.0, Z     )-(Z)*G4Electron::Electron()->GetPDGMass();
    //
    G4FastVector<G4ReactionProduct,4> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    //
    theReactionDynamics.NuclearReaction( vec, vecLen, originalIncident,
                                         targetNucleus, theAtomicMass, massVec );
    //
    G4double p = vec[0]->GetMomentum().mag();
    theParticleChange.SetMomentumChange( vec[0]->GetMomentum() *(1./p));
    theParticleChange.SetEnergyChange( vec[0]->GetKineticEnergy() );
    delete vec[0];
    //
    theParticleChange.SetNumberOfSecondaries( vecLen-1 );
    G4DynamicParticle *pd;
    for( G4int i=1; i<vecLen; ++i )
    {
      pd = new G4DynamicParticle();
      pd->SetDefinition( vec[i]->GetDefinition() );
      pd->SetMomentum( vec[i]->GetMomentum() );
      theParticleChange.AddSecondary( pd );
      delete vec[i];
    }
    //
    return &theParticleChange;
  }
 
 /* end of file */
 
