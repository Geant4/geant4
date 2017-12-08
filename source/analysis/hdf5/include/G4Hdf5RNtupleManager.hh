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

// Manager class for Hdf5 read ntuples.
// It implements functions specific to Hdf5 read ntuples.
//
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4Hdf5RNtupleManager_h
#define G4Hdf5RNtupleManager_h 1

#include "G4TRNtupleManager.hh"
#include "globals.hh"

#include "tools/hdf5/ntuple"

class G4Hdf5RNtupleManager : public G4TRNtupleManager<tools::hdf5::ntuple>
{
  friend class G4Hdf5AnalysisReader;

  protected:
    explicit G4Hdf5RNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4Hdf5RNtupleManager();

  private:
    // Methods from the templated base class
    //
    virtual G4bool GetTNtupleRow(G4TRNtupleDescription<tools::hdf5::ntuple>* ntupleDescription) final;
};    

// // inline functions

// inline G4int G4Hdf5RNtupleManager::GetNofNtuples() const
// { return fNtupleVector.size(); }

#endif

