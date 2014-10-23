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

// Structure containing the information related to rroot ntuple
//
// Author: Ivana Hrivnacova, 09/04/2014 (ivana@ipno.in2p3.fr)

#ifndef G4CsvRNtupleDescription_h
#define G4CsvRNtupleDescription_h 1

#include "tools/rcsv_ntuple"
#include "tools/ntuple_binding"

namespace tools { 
class ntuple_binding;
}

struct G4CsvRNtupleDescription
{
  G4CsvRNtupleDescription(
    tools::rcsv::ntuple* rntuple) 
    :  fNtuple(rntuple),
       fNtupleBinding(new tools::ntuple_binding()),
       fIsInitialized(false) {}

  ~G4CsvRNtupleDescription()
      { 
        delete fNtupleBinding;
        delete fNtuple;   // CHECK
      }

  tools::rcsv::ntuple* fNtuple; 
  tools::ntuple_binding* fNtupleBinding;
  G4bool fIsInitialized;
  
  private:
    // disabled (not implemented) copy constructor
    G4CsvRNtupleDescription(const G4CsvRNtupleDescription& rhs);
    // disabled (not implemented) assignement operator
    G4CsvRNtupleDescription& operator=(G4CsvRNtupleDescription& rhs);
};

#endif  
