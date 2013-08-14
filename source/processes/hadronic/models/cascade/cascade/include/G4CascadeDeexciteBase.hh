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
#ifndef G4CASCADE_DEEXCITE_BASE_HH
#define G4CASCADE_DEEXCITE_BASE_HH 1
// $Id$
//
// Semi-concrete base class for de-excitation modules, analogous to
// G4CascadeColliderBase.
//
// 20130628  M.Kelsey -- Add makeFragment() with zero momentum

#include "G4VCascadeDeexcitation.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzVector.hh"
#include <vector>

class G4CascadeCheckBalance;
class G4Fragment;


class G4CascadeDeexciteBase : public G4VCascadeDeexcitation {
public:
  G4CascadeDeexciteBase(const char* name);
  virtual ~G4CascadeDeexciteBase();

  virtual void setVerboseLevel(G4int verbose=0);

protected:
  // Decide whether to use G4BigBanger or not
  virtual G4bool explosion(const G4Fragment& target) const;
  virtual G4bool explosion(G4int A, G4int Z, G4double excitation) const;

protected:
  G4CascadeCheckBalance* balance;	// For conservation checking

  // ==> Provide interfaces to G4CascadeCheckBalance
  virtual G4bool validateOutput(const G4Fragment& target,
				G4CollisionOutput& output);

  // This is for use with G4BigBanger
  virtual G4bool validateOutput(const G4Fragment& target,
		const std::vector<G4InuclElementaryParticle>& particles);

  // This is for use with G4Fissioner
  virtual G4bool validateOutput(const G4Fragment& target,
				const std::vector<G4InuclNuclei>& fragments);

  // Interfaces between local content and G4Fragment
  void getTargetData(const G4Fragment& target);
  G4int A;		// Buffers to collect target data for all modules
  G4int Z;
  G4LorentzVector PEX;	// Four momentum of recoil in Bertini units (GeV)
  G4double EEXS;	// Excitation energy in MeV
  
  // NOTE:  Momentum passed by value so that energy/mass can be adjusted
  const G4Fragment& makeFragment(G4LorentzVector mom, G4int A, G4int Z,
				 G4double EX=0.);
  const G4Fragment& makeFragment(G4int A, G4int Z, G4double EX=0.);
  G4Fragment aFragment;	// Reusable buffer to reduce new/delete cycling

private:
  // Copying of modules is forbidden
  G4CascadeDeexciteBase(const G4CascadeDeexciteBase&);
  G4CascadeDeexciteBase& operator=(const G4CascadeDeexciteBase&);
};

#endif	/* G4CASCADE_DEEXCITE_BASE_HH */
