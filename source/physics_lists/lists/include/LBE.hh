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
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//
// This physics list is taken from the underground_physics example with small
// modifications.  It is an example of a "flat" physics list with no dependence
// on builders.  The physics covered would be suitable for a low background
// experiment including the neutron_hp package
//
//
//
// PhysicsList header
// --------------------------------------------------------------

#ifndef TLBE_h
#define TLBE_h 1

#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4VUserPhysicsList.hh"
#include "CompileTimeConstraints.hh"
#include "G4VModularPhysicsList.hh"

class G4StoppingPhysics;	// This builder encapsulate stopping processes

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

template<class T>
class TLBE: public T
{
public:
  TLBE(G4int ver = 1);
  ~TLBE();
  //  virtual ~TLBE();


public:
  virtual void SetCuts();


private:
  enum {ok = CompileTimeConstraints::IsA<T, G4VModularPhysicsList>::ok };

protected:
  // Construct particle and physics
  virtual void ConstructParticle();
  virtual void ConstructProcess();
    
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

  G4StoppingPhysics* stoppingPhysics;

  // these methods Construct particles 
  void ConstructMyBosons();
  void ConstructMyLeptons();
  void ConstructMyMesons();
  void ConstructMyBaryons();
  void ConstructMyIons();
  void ConstructMyShortLiveds();

};
#include "LBE.icc"
typedef TLBE<G4VModularPhysicsList> LBE;

#endif
