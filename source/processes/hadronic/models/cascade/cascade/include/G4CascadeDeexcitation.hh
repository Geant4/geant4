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
#ifndef G4CASCADE_DEEXCITATION_HH
#define G4CASCADE_DEEXCITATION_HH
// $Id: G4CascadeDeexcitation.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100926  M. Kelsey -- Move to new G4VCascadeDeexcitation base class.
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Drop collide() interface (base class will throw exception);
//		follow base change to deExcite() with reference;
//		add reusable G4InuclNuclei buffer, verbosity handler
// 20130622  Inherit from G4CascadeDeexciteBase

#include "G4CascadeDeexciteBase.hh"
#include "globals.hh"

class G4InuclParticle;
class G4BigBanger;
class G4NonEquilibriumEvaporator;
class G4EquilibriumEvaporator;


class G4CascadeDeexcitation : public G4CascadeDeexciteBase {
public:
  G4CascadeDeexcitation();
  virtual ~G4CascadeDeexcitation();

  virtual void setVerboseLevel(G4int verbose);

  // Interface specific to pre-compound (post-cascade) processing
  virtual void deExcite(const G4Fragment& fragment,
			G4CollisionOutput& globalOutput);

private:
  G4BigBanger* theBigBanger;
  G4NonEquilibriumEvaporator* theNonEquilibriumEvaporator;
  G4EquilibriumEvaporator* theEquilibriumEvaporator;

  G4CollisionOutput tempOutput;	// Local buffer for de-excitation stages

private:
  // Copying of modules is forbidden
  G4CascadeDeexcitation(const G4CascadeDeexcitation&);
  G4CascadeDeexcitation& operator=(const G4CascadeDeexcitation&);
};

#endif	/* G4CASCADE_DEEXCITATION_HH */
