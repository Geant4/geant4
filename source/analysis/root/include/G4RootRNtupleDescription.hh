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

#ifndef G4RootRNtupleDescription_h
#define G4RootRNtupleDescription_h 1

#include "tools/rroot/ntuple"
#include "tools/rroot/fac"
#include "tools/rroot/tree"

namespace tools { 
class ntuple_binding;
}

struct G4RootRNtupleDescription
{
  G4RootRNtupleDescription(
    tools::rroot::ntuple* rntuple,
    tools::rroot::buffer* buffer,
    tools::rroot::fac* fac,
    tools::rroot::tree* tree) 
    :  fNtuple(rntuple),
       fNtupleBinding(new tools::ntuple_binding()),
       fBuffer(buffer),
       fFac(fac),
       fTree(tree),
       fIsInitialized(false) {}

  ~G4RootRNtupleDescription()
      { 
        delete fNtupleBinding;
        delete fBuffer;
        delete fFac;
        delete fTree;
          // fNtuple is owned by the file CHECK
      }

  // deleted copy constructor
  G4RootRNtupleDescription(const G4RootRNtupleDescription& rhs) = delete;
  // deleted assignement operator
  G4RootRNtupleDescription& operator=(G4RootRNtupleDescription& rhs) = delete;

  tools::rroot::ntuple* fNtuple; 
  tools::ntuple_binding* fNtupleBinding;
  tools::rroot::buffer* fBuffer;
  tools::rroot::fac* fFac;
  tools::rroot::tree* fTree;
  G4bool fIsInitialized;
};

#endif  

