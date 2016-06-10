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
// $Id: G4EquilibriumEvaporator.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members.  Rename timeToBigBang() to override
//		base explosion().
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class
// 20100923  M. Kelsey -- Migrate to integer A and Z
// 20100925  M. Kelsey -- Remove no longer necessary explosion() interface
// 20110801  M. Kelsey -- Move "parms" buffer to data member, allocate in
//		constructor.
// 20110809  M. Kelsey -- Move "foutput" buffer to data member
// 20130129  M. Kelsey -- Move QF interpolation to global statics for
//		multi-threaded shared memory.
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Follow base class change to explosion() reference interface
// 20130622  Inherit from G4CascadeDeexciteBase, move to deExcite() interface
//		with G4Fragment
// 20130808  Use new object-version of paraMaker, for thread safety
// 20131001  Move QFinterpolation to data member to be thread local (not shared)

#ifndef G4EQUILIBRIUM_EVAPORATOR_HH
#define G4EQUILIBRIUM_EVAPORATOR_HH

#include "G4CascadeDeexciteBase.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CollisionOutput.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclSpecialFunctions.hh"


class G4EquilibriumEvaporator : public G4CascadeDeexciteBase {
public:
  G4EquilibriumEvaporator();
  virtual ~G4EquilibriumEvaporator();

  virtual void setVerboseLevel(G4int verbose);

  virtual void deExcite(const G4Fragment& target, G4CollisionOutput& output);

private: 
  // Replace base class verision with more complex conditions
  virtual G4bool explosion(G4int a, G4int z, G4double e) const;

  // NOTE:  Must redeclare base-class polymorphisms
  virtual G4bool explosion(const G4Fragment& target) const {
    return G4CascadeDeexciteBase::explosion(target);
  }

  G4bool goodRemnant(G4int a, G4int z) const; 

  G4InuclSpecialFunctions::paraMaker theParaMaker;
  G4double getE0(G4int A) const; 
  G4double getPARLEVDEN(G4int A, G4int Z) const; 
  G4double getQF(G4double x, G4double x2, G4int a, G4int z, G4double e) const;
  G4double getAF(G4double x, G4int a, G4int z, G4double e) const; 

  // Buffer for parameter sets
  std::pair<std::vector<G4double>, std::vector<G4double> > parms;
  G4CollisionOutput fission_output;

  // Interpolation object for QF
  G4CascadeInterpolator<72> QFinterp;

  G4Fissioner theFissioner;
  G4BigBanger theBigBanger;

private:
  // Copying of modules is forbidden
  G4EquilibriumEvaporator(const G4EquilibriumEvaporator&);
  G4EquilibriumEvaporator& operator=(const G4EquilibriumEvaporator&);
};        

#endif /* G4EQUILIBRIUM_EVAPORATOR_HH */


