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
// Class for configuring G4analysis for merging ntuples via MPI
//
// Author: Ivana Hrivnacova, 21/11/2018 (ivana@ipno.in2p3.fr)

#ifndef G4MPINTUPLEMERGER_HH
#define G4MPINTUPLEMERGER_HH

#include "G4MPImanager.hh"
#include "G4VMPIextraWorker.hh"

namespace toolx {
namespace mpi {
class wrmpi;    
}
}

class G4RootMpiAnalysisManager;

class G4MPIntupleMerger
{
public:
  G4MPIntupleMerger(G4int nofReducedNtupleFiles = 0,
                    G4bool rowWise = false, G4bool rowMode = true);
  ~G4MPIntupleMerger();

private:
  toolx::mpi::wrmpi* fWrmpi;
};

#endif //G4MPINTUPLEMERGER_HH
