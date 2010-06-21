#ifndef G4CASCADE_CHECK_BALANCE_HH
#define G4CASCADE_CHECK_BALANCE_HH
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
// $Id: G4CascadeCheckBalance.hh,v 1.1 2010-06-21 03:40:00 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Verify and report four-momentum conservation for collision output; uses
// same interface as collision generators.

#include "G4VCascadeCollider.hh"
#include "globals.hh"
#include "G4LorentzVector.hh"

class G4InuclParticle;
class G4CollisionOutput;

class G4CascadeCheckBalance : public G4VCascadeCollider {
public:
  G4CascadeCheckBalance(G4double relative, G4double absolute);
  virtual ~G4CascadeCheckBalance() {};

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

  G4bool energyOkay() const;
  G4bool momentumOkay() const;
  G4bool okay() const { return energyOkay() && momentumOkay(); }

  G4double deltaE() const { return (final.e() - initial.e()); }
  G4double relativeE() const { return (deltaE()==0.)?0.:deltaE()/initial.e(); }

  G4double deltaP() const { return (final.rho() - initial.rho()); }
  G4double relativeP() const { return (deltaP()==0.)?0.:deltaP()/initial.rho(); }

private:
  G4double relativeLimit;
  G4double absoluteLimit;

  G4LorentzVector initial;	// Four-vectors for computing violations
  G4LorentzVector final;
};

#endif	/* G4CASCADE_CHECK_BALANCE_HH */
