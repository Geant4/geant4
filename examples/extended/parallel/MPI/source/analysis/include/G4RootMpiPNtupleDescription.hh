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

// Structure containing the information related to Root MPI pntuple.
// This class is temporarily provided with g4mpi,
// it will be integrated in Geant4 analysis category in future.
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4RootMpiPNtupleDescription_h
#define G4RootMpiPNtupleDescription_h 1

#include "globals.hh"

#include "G4TNtupleDescription.hh"
#include "G4RootFileDef.hh"

#include "tools/ntuple_booking"
#include "tools/impi"
#include "tools/wroot/impi_ntuple"
#include "tools/wroot/ntuple"

namespace tools {
namespace wroot {
class branch;
class base_pntuple;
}
}

using RootNtupleDescription = G4TNtupleDescription<tools::wroot::ntuple, G4RootFile>;

struct G4RootMpiPNtupleDescription
{
  G4RootMpiPNtupleDescription(G4NtupleBooking* g4NtupleBooking) 
    :  fDescription(g4NtupleBooking),
       fNtuple(nullptr),
       fBasePNtuple(nullptr),
       fMainBranches(),
       fMainNtupleRank(0),
       fImpi(nullptr) {}

  ~G4RootMpiPNtupleDescription()
      {
         if ( fDescription.GetIsNtupleOwner() ) delete fNtuple;
      }    

  RootNtupleDescription fDescription;
  tools::wroot::impi_ntuple* fNtuple;
  tools::wroot::base_pntuple* fBasePNtuple;
  std::vector<tools::wroot::branch*> fMainBranches; 
  G4int  fMainNtupleRank;
  tools::impi* fImpi;
};

#endif  
