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
//    ---------- fake Test29PhysicsList header -------
//    Created by Mikhail Kossov, 7 Dec 2004 
//
// **********************************************************************

#ifndef Test29PhysicsList_h
#define Test29PhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

class Test29PhysicsList: public G4VUserPhysicsList
{
public:
  Test29PhysicsList();
  ~Test29PhysicsList();

protected:
  // Construct particle and physics process
  void ConstructParticle();
  void ConstructProcess();
  void SetCuts();

public:
  // Set/Get cut values 
  void      SetCutForGamma(G4double);
  void      SetCutForElectron(G4double);
  void      SetCutForProton(G4double);           
  G4double  GetCutForGamma() const;
  G4double  GetCutForElectron() const;
  G4double  GetCutForProton() const;
    
protected:
  // these methods Construct particles 
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructAllShortLiveds();

protected:
  // these methods Construct physics processes and register them
  void ConstructGeneral();
  void ConstructEM();

private:
  G4double cutForGamma;
  G4double cutForElectron; 
  G4double cutForProton;
};

#endif







