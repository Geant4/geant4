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
// $Id: G4CascadeFunctions.hh,v 1.8 2010-12-15 07:39:38 gunter Exp $
// GEANT4 tag: $Name: not supported by cvs2svn $
//
// 20100407  M. Kelsey -- Return particle types std::vector<> by const ref,
//		using a static variable in the function as a buffer.
// 20100505  M. Kelsey -- Use new interpolator class, drop std::pair<>, move
//		sampleFlat(...) from G4CascadeChannel, move functionality
//		to new base class, to allow data-member buffers.  Move
//		function definitions to .icc file (needed with templating).
// 20100510  M. Kelsey -- Use both summed and inclusive cross-sections for
//		multiplicity, as done in G4{Pion,Nucleon}Sampler.  Support
//		up to 9-body final states.  Add second argument specifying
//		which Sampler is used.  Move implementations to .icc file.
// 20100511  M. Kelsey -- Pass "kinds" buffer as input to getOutputPartTypes
// 20100803  M. Kelsey -- Add printing function for debugging

#ifndef G4_CASCADE_FUNCTIONS_HH
#define G4_CASCADE_FUNCTIONS_HH

#include "globals.hh"
#include "Randomize.hh"
#include <vector>


template <class DATA, class SAMP>
class G4CascadeFunctions : public SAMP {
public:
  static G4double getCrossSection(double ke) {
    return instance.findCrossSection(ke, DATA::data.tot);
  }

  static G4double getCrossSectionSum(double ke) {
    return instance.findCrossSection(ke, DATA::data.sum);
  }

  static G4int getMultiplicity(G4double ke);

  static void
  getOutgoingParticleTypes(std::vector<G4int>& kinds, G4int mult, G4double ke);

  static void printTable();

private:
  G4CascadeFunctions() : SAMP() {}
  static const G4CascadeFunctions<DATA,SAMP> instance;
};

#include "G4CascadeFunctions.icc"

#endif	/* G4_CASCADE_FUNCTIONS_HH */
