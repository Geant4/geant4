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
// $Id: G4CollisionOutput.hh,v 1.20 2010-05-21 18:07:30 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100407  M. Kelsey -- Replace ::resize(0) with ::clear()
// 20100409  M. Kelsey -- Move function code to .cc files, not inlinable
// 20100418  M. Kelsey -- Add function to boost output lists to lab frame
// 20100520  M. Kelsey -- Add function to rotate Z axis, from G4Casc.Interface

#ifndef G4COLLISION_OUTPUT_HH
#define G4COLLISION_OUTPUT_HH

#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzRotation.hh"
#include <algorithm>
#include <vector>

class G4LorentzConvertor;

class G4CollisionOutput {

public:

  G4CollisionOutput();

  G4CollisionOutput& operator=(const G4CollisionOutput& right);

  void reset();

  void addOutgoingParticle(const G4InuclElementaryParticle& particle) {
    outgoingParticles.push_back(particle);
  }

  void addOutgoingParticles(const std::vector<G4InuclElementaryParticle>& particles);

  void addTargetFragment(const G4InuclNuclei& nuclei) {
    nucleiFragments.push_back(nuclei);
  };

  void addTargetFragments(const std::vector<G4InuclNuclei>& nuclea);

  const std::vector<G4InuclElementaryParticle>& getOutgoingParticles() const {
    return outgoingParticles;
  };

  G4int numberOfNucleiFragments() const { 
    return nucleiFragments.size(); 
  };
 
  const std::vector<G4InuclNuclei>& getNucleiFragments() const {
    return nucleiFragments;
  };

  G4LorentzVector getTotalOutputMomentum() const;

  void printCollisionOutput() const;

  void boostToLabFrame(const G4LorentzConvertor& convertor);

  void rotateEvent(const G4LorentzRotation& rotate);

  void trivialise(G4InuclParticle* bullet, 
		  G4InuclParticle* target);

  void setOnShell(G4InuclParticle* bullet, 
		  G4InuclParticle* target);

  void setRemainingExitationEnergy();

  double getRemainingExitationEnergy() const { 
    return eex_rest; 
  };

  G4bool acceptable() const { 
    return on_shell; 
  };

private: 
  G4int verboseLevel;
  std::vector<G4InuclElementaryParticle> outgoingParticles;
  std::vector<G4InuclNuclei> nucleiFragments;
  G4double eex_rest;

  std::pair<std::pair<G4int,G4int>, G4int> selectPairToTune(G4double de) const; 

  G4bool on_shell;

};        

#endif // G4COLLISION_OUTPUT_HH 

