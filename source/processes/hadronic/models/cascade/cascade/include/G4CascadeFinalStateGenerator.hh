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
// $Id$
// Author:  Michael Kelsey (SLAC)
// Date:    15 April 2013
//
// Description: Subclass of models/util G4HadDecayGenerator to support
//		production of two-body and three-body final state momenta,
//		using interaction-specific distributions.
//

#ifndef G4CascadeFinalStateGenerator_hh
#define G4CascadeFinalStateGenerator_hh 1

#include "globals.hh"
#include "G4HadDecayGenerator.hh"
#include <vector>

class G4InuclElementaryParticle;
class G4CascadeFinalStateAlgorithm;


class G4CascadeFinalStateGenerator : public G4HadDecayGenerator {
public:
  G4CascadeFinalStateGenerator();
  virtual ~G4CascadeFinalStateGenerator();

  // Configure algorithm (distributions) based on interaction
  void Configure(G4InuclElementaryParticle* bullet,
		 G4InuclElementaryParticle* target,
		 const std::vector<G4int>& particle_kinds);
};

#endif	/* G4CascadeFinalStateGenerator_hh */
