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
//
// --------------------------------------------------------------
//                 GEANT 4 - ULTRA experiment example
// --------------------------------------------------------------
//
// Code developed by:
// B. Tome, M.C. Espirito-Santo, A. Trindade, P. Rodrigues 
//
//   **********************************************
//   *        UltraPhysicsList.hh
//   **********************************************
//
//    Ultra Physics List class; Standard and Low Energy EM processes are defined for
//    the relevant particles. Optical processes are declared.
//
#ifndef UltraPhysicsList_H
#define UltraPhysicsList_H 1

#include "globals.hh"
#include "G4VModularPhysicsList.hh"

class UltraPhysicsList : public G4VModularPhysicsList
{
  public:
    UltraPhysicsList();
    ~UltraPhysicsList();
 
  protected:
    // Construct particles and processes
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  private:

  // hide assignment operator
  UltraPhysicsList & operator=(const UltraPhysicsList &right);
  UltraPhysicsList(const UltraPhysicsList&);

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fOpPhysicsList;
  G4VPhysicsConstructor*  fDecayPhysicsList;

  std::vector<G4VPhysicsConstructor*> fHadronPhys;
  G4String fEmName;

  G4int fVerboseLebel;
  G4int fMaxNumPhotonStep;

  G4bool fHelIsRegisted;
  G4bool fBicIsRegisted;
  G4bool fGnucIsRegisted;
  G4bool fStopIsRegisted;
};

#endif

