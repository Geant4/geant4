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

// Template base class the histograms/profiles objects managers.

// Author: Ivana Hrivnacova, 23/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4THnManager_h
#define G4THnManager_h 1

#include "G4Threading.hh"
#include "G4Fcn.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>
#include <string_view>

class G4AnalysisManagerState;
class G4HnManager;
class G4HnInformation;

template <typename HT>
class G4THnManager
{
  public:
    G4THnManager(const G4AnalysisManagerState& state);
    G4THnManager() = delete;
    virtual ~G4THnManager();

    G4int RegisterT(const G4String& name, HT* ht, G4HnInformation* info = nullptr);

    // Reset data
    G4bool Reset();

    // Clear data
    void ClearData();

    // Return true if the H1 vector is empty
    G4bool IsEmpty() const;

    // Method for merge (MT)
    void  AddTVector(const std::vector<HT*>& tVector);

    // New method for merge
    void  Merge(G4Mutex& mergeMutex, G4THnManager<HT>* masterInstance);
            // TO DO: change to shared_ptr<G4THnManager<HT>> masterInstance

    // Get method
    // HT*  GetT(G4int id) const;
    HT*  GetT(G4int id, G4bool warn = true, G4bool onlyIfActive = true) const;

    std::vector<HT*>*   GetTVector();
    const std::vector<HT*>&  GetTVectorRef() const;
    std::vector<std::pair<HT*, G4HnInformation*>>*   GetTHnVector();
    const std::vector<std::pair<HT*, G4HnInformation*>>&  GetTHnVectorRef() const;

    // Methods to list/print histograms
    G4bool List(std::ostream& output, G4bool onlyIfActive) const;

    // Iterators
    typename std::vector<HT*>::iterator BeginT();
    typename std::vector<HT*>::iterator EndT();
    typename std::vector<HT*>::const_iterator BeginConstT() const;
    typename std::vector<HT*>::const_iterator EndConstT() const;

  protected:
    std::pair<HT*, G4HnInformation*>  GetTHnInFunction(G4int id,
                       std::string_view functionName,
                       G4bool warn = true,
                       G4bool onlyIfActive = true) const;

    HT*  GetTInFunction(G4int id,
                       std::string_view functionName,
                       G4bool warn = true,
                       G4bool onlyIfActive = true) const;

    G4int GetTId(const G4String& name, G4bool warn = true) const;

    // Methods for verbose
    G4bool IsVerbose(G4int verboseLevel) const;
    void Message(G4int level,
                 const G4String& action,
                 const G4String& objectType,
                 const G4String& objectName = "",
                 G4bool success = true) const;

    // Static data members
    static constexpr std::string_view fkClass { "G4THnManager<T>" };

    // Data members
    const G4AnalysisManagerState& fState;
    std::vector<HT*>  fTVector;
    std::vector<std::pair<HT*, G4HnInformation*>> fTHnVector;
    std::map<G4String, G4int>    fNameIdMap;
    std::shared_ptr<G4HnManager> fHnManager { nullptr };
};

// inline functions

#include "G4THnManager.icc"

#endif
