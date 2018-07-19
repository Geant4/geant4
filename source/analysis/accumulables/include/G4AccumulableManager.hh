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
// $Id: G4AccumulableManager.hh 91012 2015-06-15 10:38:07Z ihrivnac $

// The common implementation of analysis manager classes.

// Author: Ivana Hrivnacova, 07/09/2015  (ivana@ipno.in2p3.fr)

#ifndef G4AccumulableManager_h
#define G4AccumulableManager_h 1

#include "G4Accumulable.hh"
#include "globals.hh"

#include <map>
#include <vector>

class G4AnalysisManagerState;
class G4VAccumulable;

class G4AccumulableManager
{
  public:
    virtual ~G4AccumulableManager();
     
    // static methods
    static G4AccumulableManager* Instance();
   
    // Methods

    // Create accumulables
    //
    template <typename T>
    G4Accumulable<T>* 
    CreateAccumulable(const G4String& name, T value, 
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    template <typename T> 
    G4Accumulable<T>* 
    CreateAccumulable(T value, 
                      G4MergeMode mergeMode = G4MergeMode::kAddition);

    // Register existing accumulables
    // 
    // templated accumulable
    template <typename T> 
    G4bool RegisterAccumulable(G4Accumulable<T>& accumulable);
    // user defined accumulable
    G4bool RegisterAccumulable(G4VAccumulable* accumulable);

    // Access registered accumulables
    // 
    // Via name
    // templated accumulable
    template <typename T> 
    G4Accumulable<T>*  GetAccumulable(const G4String& name, G4bool warn = true) const;
    // user defined accumulable
    G4VAccumulable*  GetAccumulable(const G4String& name, G4bool warn = true) const;

    // Via id (in the order of registering)
    //  templated accumulable
    template <typename T> 
    G4Accumulable<T>*  GetAccumulable(G4int id, G4bool warn = true) const;
    // user defined accumulable
    G4VAccumulable*  GetAccumulable(G4int id, G4bool warn = true) const;
    G4int GetNofAccumulables() const;

    // Via vector iterators
    std::vector<G4VAccumulable*>::iterator Begin();
    std::vector<G4VAccumulable*>::iterator End();
    std::vector<G4VAccumulable*>::const_iterator BeginConst() const;
    std::vector<G4VAccumulable*>::const_iterator EndConst() const;

    // Methods applied to all accumulables
    void Merge();
    void Reset();

  private:
    // hide ctor requiring master/worker specification
    G4AccumulableManager(G4bool isMaster);

    // metods
    // Generate generic accumulable name: accumulableN, where N is the actual number of accumulables
    G4String GenerateName() const;
    // Check if a name is already used in a map and print a warning 
    G4bool CheckName(const G4String& name, const G4String& where) const;

    template <typename T> 
    G4Accumulable<T>*  GetAccumulable(G4VAccumulable* accumulable, G4bool warn) const;

    // constants
    const G4String kBaseName = "accumulable";

    // static data members
    static G4AccumulableManager* fgMasterInstance;
    static G4ThreadLocal G4AccumulableManager* fgInstance;    

    // data members
    std::vector<G4VAccumulable*>        fVector;
    std::map<G4String, G4VAccumulable*> fMap;
    std::vector<G4VAccumulable*>        fAccumulablesToDelete;
 };

#include "G4AccumulableManager.icc"

#endif

