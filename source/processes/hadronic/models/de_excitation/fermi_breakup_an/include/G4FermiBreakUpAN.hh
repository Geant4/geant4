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
//
// G4FermiBreakUpAN alternative de-excitation model
// by A. Novikov (January 2025)
//

#ifndef G4FERMIBREAKUP_HH
#define G4FERMIBREAKUP_HH

#include "G4FermiParticle.hh"
#include "G4FermiSplitter.hh"

#include <G4VFermiBreakUp.hh>

#include <memory>

class G4FermiBreakUpAN : public G4VFermiBreakUp
{
  public:
    using G4FermiSplitCache = G4FermiVCache<G4FermiNucleiData, G4FermiFragmentSplits>;

    G4FermiBreakUpAN();

    G4FermiBreakUpAN(std::unique_ptr<G4FermiSplitCache>&& cache);

    void Initialise() override;

    // check if the Fermi Break Up model can be used
    // mass is an effective mass of a fragment
    G4bool IsApplicable(G4int Z, G4int A, G4double eexc) const override;

    // vector of products is added to the provided vector
    // if no decay channel is found out for the primary fragment
    // then it is added to the results vector
    // if primary decays then it is deleted
    void BreakFragment(G4FragmentVector* results, G4Fragment* theNucleus) override;

    std::vector<G4FermiParticle> BreakItUp(const G4FermiParticle& nucleus) const;

  private:
    std::vector<G4FermiParticle> SelectSplit(const G4FermiParticle& particle,
                                             const G4FermiFragmentSplits& splits) const;

    mutable std::unique_ptr<G4FermiSplitCache> cache_ = nullptr;

    // improve performance, reusing allocated memory
    mutable std::vector<G4FermiFloat> weights_;
    mutable G4FermiFragmentSplits splits_;
};

#endif  // G4FERMIBREAKUP_HH
