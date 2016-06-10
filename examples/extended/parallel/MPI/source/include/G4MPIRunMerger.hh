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
#ifndef G4MPIRUNMERGER_HH_
#define G4MPIRUNMERGER_HH_

#include "G4VUserMPIrunMerger.hh"
#include "G4Run.hh"

//MPI Merger for default G4Run class

class G4MPIrunMerger : public G4VUserMPIrunMerger {
public:
  G4MPIrunMerger() : G4VUserMPIrunMerger() {}
  G4MPIrunMerger(const G4Run* ar,
                              G4int destination = G4MPImanager::kRANK_MASTER,
                              G4int verboose = 0 ) :
                                G4VUserMPIrunMerger(ar,destination,verboose) {}
protected:
  void Pack() {/*nothing do to*/}
  G4Run* UnPack() { return new G4Run; }
};



#endif /* G4MPIRUNMERGER_HH_ */
