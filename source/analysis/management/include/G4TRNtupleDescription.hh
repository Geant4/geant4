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
// Author: Ivana Hrivnacova, 20/07/2017 (ivana@ipno.in2p3.fr)

#ifndef G4TRNtupleDescription_h
#define G4TRNtupleDescription_h 1

#include "tools/ntuple_binding"

#include <map>
#include <vector>

namespace tools { 
class ntuple_binding;
}

template <typename TNTUPLE>
struct G4TRNtupleDescription
{
  G4TRNtupleDescription(TNTUPLE* rntuple) 
    :  fNtuple(rntuple),
       fNtupleBinding(new tools::ntuple_binding()),
       fIsInitialized(false),
       fIVectorBindingMap(),
       fFVectorBindingMap(),
       fDVectorBindingMap() {}

  ~G4TRNtupleDescription()
      { 
        delete fNtupleBinding;
        delete fNtuple;   // CHECK

        {for ( auto mapElement : fIVectorBindingMap ) {
          delete mapElement.first;        
        }}
        {for ( auto mapElement : fFVectorBindingMap ) {
          delete mapElement.first;        
        }}
        {for ( auto mapElement : fDVectorBindingMap ) {
          delete mapElement.first;        
        }}
      }

  // deleted copy constructor
  G4TRNtupleDescription(const G4TRNtupleDescription& rhs) = delete;
  // deleted assignement operator
  G4TRNtupleDescription& operator=(G4TRNtupleDescription& rhs) = delete;

  TNTUPLE* fNtuple; 
  tools::ntuple_binding* fNtupleBinding;
  G4bool fIsInitialized;

  // needed for XML
  std::map<TNTUPLE*, std::vector<int>* >    fIVectorBindingMap;           
  std::map<TNTUPLE*, std::vector<float>* >  fFVectorBindingMap;           
  std::map<TNTUPLE*, std::vector<double>* > fDVectorBindingMap;             
};

#endif
