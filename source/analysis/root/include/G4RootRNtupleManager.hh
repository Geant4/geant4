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
// $Id: G4RootRNtupleManager.hh 70604 2013-06-03 11:27:06Z ihrivnac $

// Manager class for Root read ntuples.
// It implements functions specific to Root read ntuples.
//
// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4RootRNtupleManager_h
#define G4RootRNtupleManager_h 1

#include "G4TRNtupleManager.hh"
#include "globals.hh"

#include "tools/rroot/ntuple"

#include <vector>

struct G4RootRNtupleDescription;

class G4RootRNtupleManager : public G4TRNtupleManager<tools::rroot::ntuple>
{
  friend class G4RootAnalysisReader;

  protected:
    explicit G4RootRNtupleManager(const G4AnalysisManagerState& state);
    virtual ~G4RootRNtupleManager();
  
  private:
    // Methods from the templated base class
    //
    virtual G4bool GetTNtupleRow(G4TRNtupleDescription<tools::rroot::ntuple>* ntupleDescription) final;
};    

#endif

