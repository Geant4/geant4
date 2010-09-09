#ifndef G4CASCADE_RECOIL_MAKER_HH
#define G4CASCADE_RECOIL_MAKER_HH
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
// $Id: G4CascadeRecoilMaker.hh,v 1.1 2010-09-09 19:11:27 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Collects generated cascade data (using Collider::collide() interface)
// and computes the nuclear recoil kinematics needed to balance the event.
//
// 20100909  M. Kelsey -- Inspired by G4CascadeCheckBalance

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"
#include "G4LorentzVector.hh"
#include <cmath>
#include <vector>

class G4CascadeCheckBalance;
class G4CascadParticle;
class G4InuclElementaryParticle;
class G4InuclNuclei;
class G4InuclParticle;


class G4CascadeRecoilMaker : public G4VCascadeCollider {
public:
  G4CascadeRecoilMaker();
  virtual ~G4CascadeRecoilMaker();

  // Standard Collider interface (non-const output "buffer")
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  // This is for use with G4IntraNucleiCascader
  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       const std::vector<G4InuclElementaryParticle>& particles,
	       const std::vector<G4CascadParticle>& cparticles);

  // Access constructed nucleus (may return null pointer complete event!)
  const G4InuclNuclei* getRecoilFragment() const { 
    return (goodRecoil() ? &theRecoilFragment : 0); }

  G4double getRecoilA() const { return recoilA; }
  G4double getRecoilZ() const { return recoilZ; }
  const G4LorentzVector& getRecoilMomentum() const { return recoilMomentum; }

  // Data quality checks
  G4bool goodRecoil() const;
  G4bool wholeEvent() const;
  G4bool unphysicalRecoil() const { return !wholeEvent() && !goodRecoil(); }

protected:
  void makeRecoilFragment();	// Convert recoil parameters into object

private:
  G4CascadeCheckBalance* balance;	// Used to do kinematics calculations

  G4double recoilA;			// Nuclear parameters of recoil
  G4double recoilZ;
  G4LorentzVector recoilMomentum;
  G4double excitationEnergy;

  G4InuclNuclei theRecoilFragment;	// Buffer! will be reused every time
};

#endif	/* G4CASCADE_RECOIL_MAKER_HH */
