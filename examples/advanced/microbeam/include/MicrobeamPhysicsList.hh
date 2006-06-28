//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// -------------------------------------------------------------------
// $Id: MicrobeamPhysicsList.hh,v 1.4 2006-06-28 13:43:07 gunter Exp $
// -------------------------------------------------------------------

#ifndef MicrobeamPhysicsList_h
#define MicrobeamPhysicsList_h 1

#include "G4VUserPhysicsList.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class MicrobeamPhysicsList: public G4VUserPhysicsList
{
public:
  MicrobeamPhysicsList();
  ~MicrobeamPhysicsList();

private:
  
  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  
protected:
  
  void ConstructParticle();
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructBaryons();

  void ConstructProcess();
  void ConstructEM();
  void ConstructHad();
  void ConstructGeneral();
  void SetCuts();
  
};
#endif



