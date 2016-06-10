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
// $Id: G4RPGReaction.hh 79697 2014-03-12 13:10:09Z gcosmo $
//
// Author: D. H. Wright
// Date:   26 May 2007
//

#ifndef G4RPGReaction_h
#define G4RPGReaction_h 1

// Class Description:
//
// Base class providing methods for various stages of the re-parameterized 
// Gheisha model calculation of primary and secondary final state momenta

#include "G4DynamicParticle.hh"
#include "G4ReactionProduct.hh"
#include "G4Nucleus.hh"
#include "G4FastVector.hh"
#include "G4HadProjectile.hh"


class G4RPGReaction 
{
public:  // with description
    
  G4RPGReaction() {}
    
  virtual ~G4RPGReaction() {}

  G4bool ReactionStage(const G4HadProjectile* /*originalIncident*/,
                       G4ReactionProduct& /*modifiedOriginal*/,
                       G4bool& /*incidentHasChanged*/,
                       const G4DynamicParticle* /*originalTarget*/,
                       G4ReactionProduct& /*targetParticle*/,
                       G4bool& /*targetHasChanged*/,
                       const G4Nucleus& /*targetNucleus*/,
                       G4ReactionProduct& /*currentParticle*/,
                       G4FastVector<G4ReactionProduct,256>& /*vec*/,
                       G4int& /*vecLen*/,
                       G4bool /*leadFlag*/,
                       G4ReactionProduct& /*leadingStrangeParticle*/);


  void AddBlackTrackParticles(const G4double /*epnb*/,
                              const G4int /*npnb*/,
                              const G4double /*edta*/,
                              const G4int /*ndta*/,
                              const G4ReactionProduct& /*modifiedOriginal*/,
                              G4int /*PinNucleus*/,
                              G4int /*NinNucleus*/,
                              const G4Nucleus& /*aNucleus*/,
			      G4FastVector<G4ReactionProduct,256>& /*vec*/,
                              G4int& /*vecLen*/ );


  G4double GenerateNBodyEvent(const G4double totalEnergy,
                              const G4bool constantCrossSection,
                              G4FastVector<G4ReactionProduct,256> &vec,
                              G4int& vecLen);
    
  G4double GenerateNBodyEventT(const G4double totalEnergy,
                               const G4bool constantCrossSection,
                               std::vector<G4ReactionProduct*>& list);
    
  void NuclearReaction(G4FastVector<G4ReactionProduct,4> &vec,
                       G4int& vecLen,
                       const G4HadProjectile* originalIncident,
                       const G4Nucleus& aNucleus,
                       const G4double theAtomicMass,
                       const G4double* massVec);
    
protected:  // with description
    
  void Rotate(const G4double numberofFinalStateNucleons,
              const G4ThreeVector& temp,
              const G4ReactionProduct& modifiedOriginal,
              const G4HadProjectile* originalIncident,
              const G4Nucleus& targetNucleus,
              G4ReactionProduct &currentParticle,
              G4ReactionProduct &targetParticle,
              G4FastVector<G4ReactionProduct,256> &vec,
              G4int& vecLen );
    
  void Defs1(const G4ReactionProduct& modifiedOriginal,
             G4ReactionProduct& currentParticle,
             G4ReactionProduct& targetParticle,
             G4FastVector<G4ReactionProduct,256> &vec,
             G4int& vecLen);
    
  std::pair<G4int, G4int> GetFinalStateNucleons(
             const G4DynamicParticle* originalTarget,
             const G4FastVector<G4ReactionProduct,256>& vec,
             const G4int& vecLen );

  void MomentumCheck(const G4ReactionProduct &modifiedOriginal,
                     G4ReactionProduct &currentParticle,
                     G4ReactionProduct &targetParticle,
                     G4FastVector<G4ReactionProduct,256> &vec,
                     G4int& vecLen);
    
  G4double normal();
    
  G4ThreeVector Isotropic(const G4double&);
 
};
 
#endif
 
