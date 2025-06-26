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
// Author: Christian Velten (2025)
#ifndef G4VMOLECULEREACTIONCOUNTER_HH
#define G4VMOLECULEREACTIONCOUNTER_HH 1

#include "G4VMoleculeCounterInternalBase.hh"

#include <map>
#include <memory>

class G4Track;
class G4DNAMolecularReactionData;

//------------------------------------------------------------------------------

class G4VMoleculeReactionCounter : public G4VMoleculeCounterInternalBase
{
    template<typename>
    friend class G4VUserMoleculeReactionCounter;

  public:
    enum MoleculeReactionCounterType
    {
      Other,
      Basic,
    };
    struct G4VMoleculeReactionCounterIndex
    {
        virtual ~G4VMoleculeReactionCounterIndex() = default;
        virtual G4bool operator<(G4VMoleculeReactionCounterIndex const&) const = 0;
        virtual G4bool operator==(G4VMoleculeReactionCounterIndex const&) const = 0;
        virtual G4String GetInfo() const = 0;
        virtual const G4DNAMolecularReactionData* GetReactionData() const = 0;
    };

  private:
    G4VMoleculeReactionCounter();
    G4VMoleculeReactionCounter(const G4String&, MoleculeReactionCounterType =
                                                  MoleculeReactionCounterType::Basic);
    G4VMoleculeReactionCounter(G4VMoleculeReactionCounter const&) = delete;
    void operator=(G4VMoleculeReactionCounter const& x) = delete;

  public:
    ~G4VMoleculeReactionCounter() override = default;

  public:
    // virtual std::unique_ptr<G4VMoleculeReactionCounterIndex> BuildIndex(const G4Track*, const G4Track*, const G4DNAMolecularReactionData*) const = 0;
    virtual std::unique_ptr<G4VMoleculeReactionCounterIndex> BuildSimpleIndex(const G4DNAMolecularReactionData*) const = 0;

    virtual void RecordReaction(std::unique_ptr<G4VMoleculeReactionCounterIndex>, G4double, G4int = 1) = 0;

    virtual std::set<const G4DNAMolecularReactionData*> GetRecordedReactions() const = 0;

  protected:
    MoleculeReactionCounterType fType{MoleculeReactionCounterType::Basic};

  public:
	  MoleculeReactionCounterType GetType() const;
};

//------------------------------------------------------------------------------

inline G4VMoleculeReactionCounter::MoleculeReactionCounterType G4VMoleculeReactionCounter::GetType() const
{
  return fType;
}

#endif
