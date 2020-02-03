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
////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  File:   G4LightTargetCollider.hh                                          // 
//  Date:   30 September 2019                                                 //
//  Author: Dennis Wright (SLAC)                                              //
//                                                                            //
//  Description: model for collision of elementary particles with light       //
//               targets (H, D, T, 3He)                                       //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#ifndef G4LIGHT_TARGET_COLLIDER_HH
#define G4LIGHT_TARGET_COLLIDER_HH

#include "G4CascadeColliderBase.hh"
#include "G4CascadeFinalStateGenerator.hh"
#include "G4CollisionOutput.hh"

class G4CascadParticle;
class G4ElementaryParticleCollider;
class G4InuclParticle;
class G4KineticTrackVector;

typedef std::pair<G4InuclElementaryParticle, G4InuclElementaryParticle> NucleonPair;
typedef std::vector<G4InuclElementaryParticle> ScatteringProducts;

class G4LightTargetCollider : public G4CascadeColliderBase {
public:
  G4LightTargetCollider();
  virtual ~G4LightTargetCollider();

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& globalOutput);

  void setVerboseLevel(G4int verbose=0);

private: 
  G4ElementaryParticleCollider* theElementaryParticleCollider;

  G4CollisionOutput output;		// Secondaries from main cascade

private:
  // Copying of modules is forbidden
  G4LightTargetCollider(const G4LightTargetCollider&);
  G4LightTargetCollider& operator=(const G4LightTargetCollider&);

  G4double GammaDCrossSection(G4double /*kineticEnergy*/);

  G4CascadeFinalStateGenerator fsGen;

  NucleonPair AbsorptionOnDeuteron(G4InuclParticle* bullet);

  ScatteringProducts SingleNucleonScattering(const G4InuclElementaryParticle& projectile,
                                             const G4InuclElementaryParticle& targetNucleon);

  G4double mP;   // proton mass
  G4double mN;   // neutron mass
  G4double mD;   // deuteron mass
  G4double pFermiD;  // deuteron Fermi momentum (GeV/c)  
};        

#endif /* G4LIGHT_TARGET_COLLIDER_HH */


