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
// Author: Christian Velten (2025)

#ifndef G4MoleculeReactionCounter_hh
#define G4MoleculeReactionCounter_hh 1

#include "G4DNAMolecularReactionTable.hh"
#include "G4VUserMoleculeReactionCounter.hh"

//------------------------------------------------------------------------------

struct G4MoleculeReactionCounterIndex : public G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex
{
    const G4DNAMolecularReactionData* ReactionData;

    G4MoleculeReactionCounterIndex() : ReactionData(nullptr) {}
    explicit G4MoleculeReactionCounterIndex(const G4DNAMolecularReactionData* reaction) : ReactionData(reaction) {}
    ~G4MoleculeReactionCounterIndex() override = default;

    G4bool operator<(G4VMoleculeReactionCounterIndex const& other) const override
    {
      return std::less{}(ReactionData, static_cast<const G4MoleculeReactionCounterIndex&>(other).ReactionData);
    }
    G4bool operator==(G4VMoleculeReactionCounterIndex const& other) const override
    {
      return std::equal_to{}(ReactionData, static_cast<const G4MoleculeReactionCounterIndex&>(other).ReactionData);
    }
    G4String FormattedReactionString(const G4DNAMolecularReactionData* reactionData) const;

    G4String GetInfo() const override
    {
      G4String null = "This reaction data is null";
      if (ReactionData == nullptr)
        return null;
      else
        return FormattedReactionString(ReactionData);
    }
    const G4DNAMolecularReactionData* GetReactionData() const override { return ReactionData; }
};

class G4MoleculeReactionCounter : public G4VUserMoleculeReactionCounter<G4MoleculeReactionCounterIndex>
{
    //----------------------------------------------------------------------------
  public:
    G4MoleculeReactionCounter();
    G4MoleculeReactionCounter(G4String);
    ~G4MoleculeReactionCounter() override = default;

    void InitializeUser() override;

  public:
    std::unique_ptr<G4VMoleculeReactionCounterIndex> BuildSimpleIndex(const G4DNAMolecularReactionData*) const override;
};

#endif  // G4MoleculeReactionCounter_hh
