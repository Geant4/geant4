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
#ifndef G4TBBTASK_HH
#define G4TBBTASK_HH

#include <tbb/task.h>
#include "G4Types.hh"
#include "G4String.hh"

class G4VtbbJob;

class G4tbbTask : public tbb::task {
public:
  G4tbbTask(G4VtbbJob* job, 
            G4int i_event,  /* G4int n_select, const G4String& msg,*/ 
            G4int seedslenght);
  virtual ~G4tbbTask();
  tbb::task* execute();
private:
  //const G4int select;
  const G4int event;
  //const G4String msg;
  long* seeds;
  G4VtbbJob* job;
}; 
#endif //G4TBBTASK_HH
