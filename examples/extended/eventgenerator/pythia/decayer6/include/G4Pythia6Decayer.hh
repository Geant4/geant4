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
// $Id: G4Pythia6Decayer.hh 100687 2016-10-31 11:20:33Z gcosmo $
//
/// \file eventgenerator/pythia/decayer6/include/G4Pythia6Decayer.hh
/// \brief Definition of the G4Pythia6Decayer class
//
#ifndef G4_PYTHIA6_DECAYER_H
#define G4_PYTHIA6_DECAYER_H

#include "G4VExtDecayer.hh"
#include "G4Pythia6DecayerMessenger.hh"
#include "Pythia6.hh"
#include "EDecayType.hh"

#include "globals.hh"

struct Pythia6Particle;
class G4Track;
class G4DecayProducts;

/// Pythia6 decayer
///
/// Implements the G4VExtDecayer abstract class using the Pythia6 interface.
/// According to TPythia6Decayer class in Root:
/// http://root.cern.ch/
/// see http://root.cern.ch/root/License.html

class G4Pythia6Decayer : public G4VExtDecayer
{
  public:

    G4Pythia6Decayer();
    virtual ~G4Pythia6Decayer();

    virtual G4DecayProducts* ImportDecayProducts(const G4Track& track);
    
    void ForceDecayType(EDecayType decayType);
    void SetVerboseLevel(G4int verboseLevel) { fVerboseLevel =  verboseLevel; }
    
  private:

    /// Not implemented
    G4Pythia6Decayer(const G4Pythia6Decayer& right);
    /// Not implemented
    G4Pythia6Decayer& operator=(const G4Pythia6Decayer& right);
    
    G4ParticleDefinition*
    GetParticleDefinition(const Pythia6Particle* p,G4bool warn = true) const;
    G4DynamicParticle* CreateDynamicParticle(const Pythia6Particle* p) const;
    G4ThreeVector GetParticlePosition(const Pythia6Particle* particle) const;
    G4ThreeVector GetParticleMomentum(const Pythia6Particle* particle) const; 
                           
    G4int CountProducts(G4int channel, G4int particle);
    void  ForceParticleDecay(G4int particle, G4int product, G4int mult);
    void  ForceParticleDecay(G4int particle, 
                             G4int* products, G4int* mult, G4int npart);
    void  ForceHadronicD();
    void  ForceOmega();
    void  ForceDecay(EDecayType decayType);
    
    
    void  Decay(G4int pdg, const CLHEP::HepLorentzVector& p);
    G4int ImportParticles(ParticleVector* particles);
    
    static const EDecayType fgkDefaultDecayType; ///< default decay type

    G4Pythia6DecayerMessenger fMessenger;  ///< command messenger 
    G4int            fVerboseLevel;        ///< verbose level
    EDecayType       fDecayType;           ///< selected decay type
    ParticleVector*  fDecayProductsArray ; ///< array of decay products
};

// ----------------------------------------------------------------------------

#endif
