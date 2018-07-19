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

// Manager class for Csv read ntuples.
// It implements functions specific to Csv read ntuples.
//
// Author: Ivana Hrivnacova, 25/07/2014 (ivana@ipno.in2p3.fr)

#ifndef G4CsvRNtupleManager_h
#define G4CsvRNtupleManager_h 1

#include "G4TRNtupleManager.hh"
#include "globals.hh"

#include "tools/rcsv_ntuple"

class G4CsvRNtupleManager : public G4TRNtupleManager<tools::rcsv::ntuple>
{
  friend class G4CsvAnalysisReader;

  protected:
    explicit G4CsvRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4CsvRNtupleManager();
  
  private:
    // Methods from the templated base class
    //
    virtual G4bool GetTNtupleRow(G4TRNtupleDescription<tools::rcsv::ntuple>* ntupleDescription) final;
};    

#endif

