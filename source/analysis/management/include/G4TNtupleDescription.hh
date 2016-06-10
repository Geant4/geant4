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

// Structure templated containing the information related to an ntuple
// for all output types.
//
// Author: Ivana Hrivnacova, 19/06/2015  (ivana@ipno.in2p3.fr)

#ifndef G4TNtupleDescription_h
#define G4TNtupleDescription_h 1

#include "globals.hh"

#include "tools/ntuple_booking"

#include <fstream>

template <typename TNTUPLE>
struct G4TNtupleDescription
{
  G4TNtupleDescription() 
    :  fFile(nullptr),
       fNtuple(nullptr),
       fNtupleBooking(),
       fActivation(true),
       fIsNtupleOwner(true) {}

  ~G4TNtupleDescription()
      {  delete fFile;
         if ( fIsNtupleOwner ) delete fNtuple;
      }    

  std::ofstream* fFile;    
  TNTUPLE* fNtuple; 
  tools::ntuple_booking fNtupleBooking; 
  G4bool fActivation;
  G4bool fIsNtupleOwner;
};

#endif  
