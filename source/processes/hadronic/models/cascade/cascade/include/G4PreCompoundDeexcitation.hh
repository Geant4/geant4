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
#ifndef G4PRECOMPOUND_DEEXCITATION_HH
#define G4PRECOMPOUND_DEEXCITATION_HH
// $Id: G4PreCompoundDeexcitation.hh 71942 2013-06-28 19:08:11Z mkelsey $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.
//
// 20100922  M. Kelsey -- Remove convertFragment() function, pass buffer
//		instead of copying
// 20100926  M. Kelsey -- Move to new G4VCascadeDeexcitation base class.
// 20130620  Address Coverity complaint about missing copy actions
// 20130621  Drop collide() interface (base class will throw exception);
//		follow base change to deExcite() with reference.
// 20130622  Inherit from G4CascadeDeexciteBase, add verbosity interface
//		to pass to PreCompound

#include "G4CascadeDeexciteBase.hh"
#include "globals.hh"

class G4InuclNuclei;
class G4InuclParticle;
class G4ExcitationHandler;
class G4VPreCompoundModel;


class G4PreCompoundDeexcitation : public G4CascadeDeexciteBase {

public:
  G4PreCompoundDeexcitation();
  virtual ~G4PreCompoundDeexcitation();

  virtual void setVerboseLevel(G4int verbose);

  // Interface specific to pre-compound (post-cascade) processing
  virtual void deExcite(const G4Fragment& fragment,
			G4CollisionOutput& globalOutput);

private:
  G4ExcitationHandler* theExcitationHandler;
  G4VPreCompoundModel* theDeExcitation;

private:
  // Copying of modules is forbidden
  G4PreCompoundDeexcitation(const G4PreCompoundDeexcitation&);
  G4PreCompoundDeexcitation& operator=(const G4PreCompoundDeexcitation&);
};

#endif	/* G4PRECOMPOUND_DEEXCITATION_HH */
