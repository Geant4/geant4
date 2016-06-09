//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
//
 // Hadronic Process: Deuteron Inelastic Process
 // J.L. Chuma, TRIUMF, 25-Feb-1997
 // Last modified: 27-Mar-1997
 // J.L. Chuma, 08-May-2001: Update original incident passed back in vec[0]
 //                          from NuclearReaction
 //
#include "G4LEDeuteronInelastic.hh"
#include "Randomize.hh"
#include "G4Electron.hh"

 G4HadFinalState *
  G4LEDeuteronInelastic::ApplyYourself( const G4HadProjectile &aTrack,
                                        G4Nucleus &targetNucleus )
  { 
    theParticleChange.Clear();
    const G4HadProjectile *originalIncident = &aTrack;
    
    if( verboseLevel > 1 )
    {
      const G4Material *targetMaterial = aTrack.GetMaterial();
      G4cout << "G4LEDeuteronInelastic::ApplyYourself called" << G4endl;
      G4cout << "kinetic energy = " << originalIncident->GetKineticEnergy()/MeV << "MeV, ";
      G4cout << "target material = " << targetMaterial->GetName() << ", ";
    }
    
    // Work-around for lack of model above 100 MeV
    if( originalIncident->GetKineticEnergy()/MeV > 100. ||
        originalIncident->GetKineticEnergy() <= 0.1*MeV )
    {
      theParticleChange.SetStatusChange(isAlive);
      theParticleChange.SetEnergyChange(aTrack.GetKineticEnergy());
      theParticleChange.SetMomentumChange(aTrack.Get4Momentum().vect().unit()); 
      return &theParticleChange;      
    }

    G4double A = targetNucleus.GetN();
    G4double Z = targetNucleus.GetZ();
    G4double theAtomicMass = targetNucleus.AtomicMass( A, Z );
    G4double massVec[9];
    massVec[0] = targetNucleus.AtomicMass( A+2.0, Z+1.0 );
    massVec[1] = targetNucleus.AtomicMass( A+1.0, Z+1.0 );
    massVec[2] = targetNucleus.AtomicMass( A+1.0, Z     );
    massVec[3] = theAtomicMass;
    massVec[4] = 0.;
    if (A > 1.0 && A-1.0 > Z) 
        massVec[4] = targetNucleus.AtomicMass( A-1.0, Z     );
    massVec[5] = 0.;
    if (A > 2.0 && Z > 1.0 && A-2.0 > Z-1.0) 
        massVec[5] = targetNucleus.AtomicMass( A-2.0, Z-1.0 );
    massVec[6] = 0.;
    if (A > Z+1.0) 
        massVec[6] = targetNucleus.AtomicMass( A,     Z+1.0 );
    massVec[7] = massVec[3];
    massVec[8] = 0.;
    if (Z > 1.0) massVec[8] = targetNucleus.AtomicMass( A,Z-1.0 );
    
    G4FastVector<G4ReactionProduct,4> vec;  // vec will contain the secondary particles
    G4int vecLen = 0;
    vec.Initialize( 0 );
    
    theReactionDynamics.NuclearReaction( vec, vecLen, originalIncident,
                                         targetNucleus, theAtomicMass, massVec );
    //
    G4double p = vec[0]->GetMomentum().mag();
    theParticleChange.SetMomentumChange( vec[0]->GetMomentum() * (1.0/p)  );
    theParticleChange.SetEnergyChange( vec[0]->GetKineticEnergy() );
    delete vec[0];
    //
    G4DynamicParticle *pd;
    for( G4int i=1; i<vecLen; ++i )
    {
      pd = new G4DynamicParticle();
      pd->SetDefinition( vec[i]->GetDefinition() );
      pd->SetMomentum( vec[i]->GetMomentum() );
      theParticleChange.AddSecondary( pd );
      delete vec[i];
    }
    return &theParticleChange;
  }
 
 /* end of file */
 
