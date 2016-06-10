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
// $Id: G4THnManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Template base class the histograms/profiles objects managers.

// Author: Ivana Hrivnacova, 23/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4THnManager_h
#define G4THnManager_h 1

#include "G4Fcn.hh"
#include "G4BinScheme.hh"
#include "globals.hh"

#include <vector>
#include <map>
#include <memory>

class G4AnalysisManagerState;
class G4HnManager;

template <typename T>
class G4THnManager
{
  public:
    G4THnManager(const G4AnalysisManagerState& state,
                 const G4String& hnType);
    virtual ~G4THnManager();

    // Reset data
    G4bool Reset();
    // Return true if the H1 vector is empty
    G4bool IsEmpty() const;   

  protected:

    // Method for merge (MT)
    void  AddTVector(const std::vector<T*>& tVector);

    // Iterators
    typename std::vector<T*>::iterator BeginT();
    typename std::vector<T*>::iterator EndT();
    typename std::vector<T*>::const_iterator BeginConstT() const;
    typename std::vector<T*>::const_iterator EndConstT() const;

    T*  GetTInFunction(G4int id, 
                       G4String functionName,
                       G4bool warn = true,
                       G4bool onlyIfActive = true) const;

    G4int RegisterT(T* t, const G4String& name);

    G4int GetTId(const G4String& name, G4bool warn = true) const;

    // data members
    const G4AnalysisManagerState& fState;
    std::vector<T*>  fTVector;
    std::map<G4String, G4int>    fNameIdMap;                
    std::shared_ptr<G4HnManager> fHnManager;
};

// inline functions

#include "G4THnManager.icc"

#endif

