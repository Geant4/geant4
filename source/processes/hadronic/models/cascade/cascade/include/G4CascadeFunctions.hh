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
// $Id: G4CascadeFunctions.hh 67796 2013-03-08 06:18:39Z mkelsey $
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
// 20110719  M. Kelsey -- Add inheritance from non-template base for factory,
//		change static's to virtual (no more direct access)
// 20110725  M. Kelsey -- Move ctor to .icc file for registration in lookup
// 20110923  M. Kelsey -- Add optional ostream& argument to printTable()

#ifndef G4_CASCADE_FUNCTIONS_HH
#define G4_CASCADE_FUNCTIONS_HH

#include "G4CascadeChannel.hh"
#include "globals.hh"
#include "Randomize.hh"
#include <vector>


template <class DATA, class SAMP>
class G4CascadeFunctions : public G4CascadeChannel, public SAMP {
public:
  G4CascadeFunctions();
  virtual ~G4CascadeFunctions() {}

  virtual G4double getCrossSection(double ke) const {
    return this->findCrossSection(ke, DATA::data.tot);
  }

  virtual G4double getCrossSectionSum(double ke) const {
    return this->findCrossSection(ke, DATA::data.sum);
  }

  virtual G4int getMultiplicity(G4double ke) const;

  virtual void getOutgoingParticleTypes(std::vector<G4int>& kinds,
				       G4int mult, G4double ke) const;

  virtual void printTable(std::ostream& os=G4cout) const;
};

#include "G4CascadeFunctions.icc"

#endif	/* G4_CASCADE_FUNCTIONS_HH */
