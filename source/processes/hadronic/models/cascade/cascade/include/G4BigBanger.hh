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
// $Id: G4BigBanger.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100315  M. Kelsey -- Remove "using" directive and unnecessary #includes.
// 20100407  M. Kelsey -- Replace std::vector<> returns with data members.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class
// 20100519  M. Kelsey -- Get rid of proton and neutron masses as arguments!
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100726  M. Kelsey -- Move std::vector<> buffer to .hh file
// 20100928  M. Kelsey -- Migrate to integer A and Z
// 20130620  Address Coverity complaint about missing copy actions
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment

#ifndef G4BIG_BANGER_HH
#define G4BIG_BANGER_HH

#include "G4CascadeDeexciteBase.hh"
#include "G4InuclElementaryParticle.hh"
#include <vector>


class G4BigBanger : public G4CascadeDeexciteBase {
public:
  G4BigBanger();
  virtual ~G4BigBanger() {};

  virtual void deExcite(const G4Fragment& target, G4CollisionOutput& output);

private: 
  void generateBangInSCM(G4double etot, G4int a, G4int z);
  void generateMomentumModules(G4double etot, G4int a, G4int z); 
  G4double xProbability(G4double x, G4int a) const; 
  G4double maxProbability(G4int a) const;
  G4double generateX(G4int ia, G4double promax) const; 

  // Buffers for big-bang results
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4double> momModules;
  std::vector<G4LorentzVector> scm_momentums;

private:
  // Copying of modules is forbidden
  G4BigBanger(const G4BigBanger&);
  G4BigBanger& operator=(const G4BigBanger&);
};        

#endif /* G4BIG_BANGER_HH */











