//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
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
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// PhysicsList header
// --------------------------------------------------------------

#ifndef DMXPhysicsList_h
#define DMXPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class DMXPhysicsList: public G4VUserPhysicsList
{
public:
  DMXPhysicsList();
  virtual ~DMXPhysicsList();

protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
  
  virtual void SetCuts();
  
  // these methods Construct physics processes and register them
  virtual void ConstructGeneral();
  virtual void ConstructEM();
  virtual void ConstructHad();
  virtual void ConstructOp();


  /*
  // these methods Construct all particles in each category
  virtual void ConstructAllBosons();
  virtual void ConstructAllLeptons();
  virtual void ConstructAllMesons();
  virtual void ConstructAllBaryons();
  virtual void ConstructAllIons();
  virtual void ConstructAllShortLiveds();
  */

  virtual void AddTransportation();

private:
  G4int VerboseLevel;
  G4int OpVerbLevel;

  G4double cutForGamma;
  G4double cutForElectron;
  G4double cutForPositron;
  G4double cutForProton;
  G4double cutForAlpha;
  G4double cutForGenericIon;
  G4double cutForOpticalPhoton;

  // these methods Construct particles 
  void ConstructMyBosons();
  void ConstructMyLeptons();
  void ConstructMyMesons();
  void ConstructMyBaryons();
  void ConstructMyIons();
  void ConstructMyShortLiveds();

};

#endif
