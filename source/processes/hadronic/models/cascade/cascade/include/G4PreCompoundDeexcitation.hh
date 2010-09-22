#ifndef G4PRECOMPOUND_DEEXCITATION_HH
#define G4PRECOMPOUND_DEEXCITATION_HH
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
// $Id: G4PreCompoundDeexcitation.hh,v 1.1 2010-09-22 22:17:08 yarba Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// Takes an arbitrary excited or unphysical nuclear state and produces
// a final state with evaporated particles and (possibly) a stable nucleus.

#include "G4CascadeColliderBase.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"

// class G4InuclParticle;
class G4ExcitationHandler;
class G4VPreCompoundModel;

class G4PreCompoundDeexcitation : public G4CascadeColliderBase {

public:
  G4PreCompoundDeexcitation();
  virtual ~G4PreCompoundDeexcitation();

  // Standard Collider interface (bullet is not used in this case)
  //
  void collide(G4InuclParticle* /*bullet*/, 
               G4InuclParticle* target,
	       G4CollisionOutput& globalOutput);


private:

   void convertFragment( G4InuclNuclei* rfrag, std::vector<G4InuclElementaryParticle> particles );
   void getDeExcitedFragments( G4InuclNuclei* rfrag, std::vector<G4InuclElementaryParticle> particles );

   // data members
   //
   G4ExcitationHandler* theExcitationHandler;
   G4VPreCompoundModel* theDeExcitation;

  G4CollisionOutput output;	// Local buffer for de-excitation stages
};

#endif	/* G4PRECOMPOUND_DEEXCITATION_HH */
