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
#include "G4HadProjectile.hh"

enum{ GHADLISTSIZE=256};

 class G4ReactionDynamics 
 {
 public:
    
    G4ReactionDynamics() {}
    
    virtual ~G4ReactionDynamics() {}
    
    virtual G4double FindInelasticity()
    { return 0.0; }
    
    virtual G4double FindTimeDelay()
    { return 0.0; }
    
    G4bool GenerateXandPt(                   // derived from GENXPT
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen,
     G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4HadProjectile *originalIncident,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4DynamicParticle* originalTarget,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged, 
     G4bool &targetHasChanged,
     G4bool leadFlag,
     G4ReactionProduct &leadingStrangeParticle );
    
    void SuppressChargedPions(
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen,
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged,
     G4bool &targetHasChanged );
      
    G4bool TwoCluster(                       // derived from TWOCLU
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen,
     G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4HadProjectile *originalIncident,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     const G4DynamicParticle* originalTarget,
     const G4Nucleus &targetNucleus,
     G4bool &incidentHasChanged, 
     G4bool &targetHasChanged,
     G4bool leadFlag,
     G4ReactionProduct &leadingStrangeParticle );
    
    void TwoBody(                         // derived from TWOB
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
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
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen );
    
    void ProduceStrangeParticlePairs(
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen,
     const G4ReactionProduct &modifiedOriginal,
     const G4DynamicParticle *originalTarget, 
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4bool &incidentHasChanged,
     G4bool &targetHasChanged );
    
    void NuclearReaction(                     // derived from NUCREC
     G4FastVector<G4ReactionProduct,4> &vec,
     G4int &vecLen,
     const G4HadProjectile *originalIncident,
     const G4Nucleus &aNucleus,
     const G4double theAtomicMass,
     const G4double *massVec );
    
 private:
    
    void Rotate(
     const G4double numberofFinalStateNucleons,
     const G4ThreeVector &temp,
     const G4ReactionProduct &modifiedOriginal, // Fermi motion & evap. effect included
     const G4HadProjectile *originalIncident,
     const G4Nucleus &targetNucleus,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen );
    
    void Defs1(
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
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
     G4int PinNucleus,
     G4int NinNucleus,
     const G4Nucleus &aNucleus,
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen );
    
    std::pair<G4int, G4int> GetFinalStateNucleons(
     const G4DynamicParticle* originalTarget,
     const G4FastVector<G4ReactionProduct,GHADLISTSIZE>& vec,
     const G4int& vecLen );

    void MomentumCheck(
     const G4ReactionProduct &modifiedOriginal,
     G4ReactionProduct &currentParticle,
     G4ReactionProduct &targetParticle,
     G4FastVector<G4ReactionProduct,GHADLISTSIZE> &vec,
     G4int &vecLen );
    
    G4double normal();
    
    G4int Poisson( G4double x );

 };
 
#endif
 
