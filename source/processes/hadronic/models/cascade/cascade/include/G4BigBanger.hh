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
// $Id: G4BigBanger.hh,v 1.14 2010-06-25 09:41:48 gunter Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100407  M. Kelsey -- Replace std::vector<> returns with data members.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class
// 20100519  M. Kelsey -- Get rid of proton and neutron masses as arguments!

#ifndef G4BIG_BANGER_HH
#define G4BIG_BANGER_HH

#include "G4VCascadeCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include <vector>

class G4CollisionOutput;


class G4BigBanger : public G4VCascadeCollider {
public:
  G4BigBanger();
  virtual ~G4BigBanger() {};

  void collide(G4InuclParticle* bullet, G4InuclParticle* target,
	       G4CollisionOutput& output);

private: 
  void generateBangInSCM(G4double etot, G4double a, G4double z);

  void generateMomentumModules(G4double etot, G4double a, G4double z); 

  G4double xProbability(G4double x, 
			G4int ia) const; 

  G4double maxProbability(G4double a) const;

  G4double generateX(G4int ia, 
		     G4double a, 
		     G4double promax) const; 

  // Buffers for big-bang results
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4double> momModules;
};        

#endif /* G4BIG_BANGER_HH */











