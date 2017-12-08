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

// Structure containing the information related to Root pntuple.
//
// Author: Ivana Hrivnacova, 04/10/2016  (ivana@ipno.in2p3.fr)

#ifndef G4RootPNtupleDescription_h
#define G4RootPNtupleDescription_h 1

#include "globals.hh"

#include "tools/ntuple_booking"
#include "tools/wroot/mt_ntuple_row_wise"
#include "tools/wroot/mt_ntuple_column_wise"

#include <fstream>

namespace tools {
namespace wroot {
class file;
class branch;
}
}

struct G4RootPNtupleDescription
{
  G4RootPNtupleDescription() 
    :  fFile(nullptr),
       fNtuple(nullptr),
       fBasePNtuple(nullptr),
       fMainBranches(),
       fNtupleBooking(),
       fActivation(true),
       fIsNtupleOwner(true) {}

  ~G4RootPNtupleDescription()
      {
         if ( fIsNtupleOwner ) delete fNtuple;
      }    

  tools::wroot::file* fFile;    
  tools::wroot::imt_ntuple* fNtuple;
  tools::wroot::base_pntuple* fBasePNtuple;
  std::vector<tools::wroot::branch*> fMainBranches; 
  tools::ntuple_booking fNtupleBooking; 
  G4bool fActivation;
  G4bool fIsNtupleOwner;
};

#endif  
